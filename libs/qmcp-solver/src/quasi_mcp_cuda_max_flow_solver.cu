#include "config.h"

#ifdef CUDA_ENABLED
#include <list>

#include "bam-api/bam_api.hpp"
#include "bam-api/paired_reads.hpp"
#include "qmcp-solver/cuda_helpers.cuh"
#include "qmcp-solver/quasi_mcp_cuda_max_flow_solver.hpp"
#include "qmcp-solver/solver.hpp"

__global__ void push_relabel_kernel(
    uint32_t kernel_cycles, uint32_t nodes_count,
    qmcp::QuasiMcpCudaMaxFlowSolver::Excess* excess_func,
    qmcp::QuasiMcpCudaMaxFlowSolver::Label* label_func,
    qmcp::QuasiMcpCudaMaxFlowSolver::Capacity* residual_capacity,
    const qmcp::QuasiMcpCudaMaxFlowSolver::NeighborInfoIndex* neighbors_start_ind,
    const qmcp::QuasiMcpCudaMaxFlowSolver::NeighborInfoIndex* neighbors_end_ind,
    const qmcp::QuasiMcpCudaMaxFlowSolver::Node* neighbors,
    const qmcp::QuasiMcpCudaMaxFlowSolver::NeighborInfoIndex* inversed_edge_offset) {
    qmcp::QuasiMcpCudaMaxFlowSolver::Node node = blockDim.x * blockIdx.x + threadIdx.x;
    if (node >= nodes_count - 2) {
        return;
    }

    for (uint32_t cycle = 0; cycle < kernel_cycles; ++cycle) {
        // Ensure that current node is active
        if (excess_func[node] <= 0 || label_func[node] >= nodes_count) {
            continue;
        }

        // Find neighbor with minimum label
        qmcp::QuasiMcpCudaMaxFlowSolver::Label min_label = UINT32_MAX;
        qmcp::QuasiMcpCudaMaxFlowSolver::NeighborInfoIndex neighbor_info_ind = UINT32_MAX;
        bool is_forward = false;

        for (qmcp::QuasiMcpCudaMaxFlowSolver::NeighborInfoIndex i = neighbors_start_ind[node];
             i <= neighbors_end_ind[node]; ++i) {
            if (residual_capacity[i] == 0) {
                // Skip edges that are not part of the residual graph
                continue;
            }

            qmcp::QuasiMcpCudaMaxFlowSolver::Node curr_neighbor = neighbors[i];
            qmcp::QuasiMcpCudaMaxFlowSolver::Label neighbor_label = label_func[curr_neighbor];

            // The is_forward heuristic prefers the edges that are forward
            if (min_label > neighbor_label) {
                min_label = neighbor_label;
                neighbor_info_ind = i;
                is_forward = curr_neighbor > node;
            } else if (min_label == neighbor_label && !is_forward && curr_neighbor > node) {
                is_forward = true;
                neighbor_info_ind = i;
            }
        }

        // Do push-relabel
        if (min_label >= label_func[node]) {
            // Preform relabel operation
            label_func[node] = min_label + 1;
        } else {
            // Perform push operation
            uint32_t delta = static_cast<uint32_t>(excess_func[node]);  // excess_func[node] > 0
            if (delta > residual_capacity[neighbor_info_ind]) {
                delta = residual_capacity[neighbor_info_ind];
            }

            atomicAdd(&residual_capacity[neighbors_start_ind[neighbors[neighbor_info_ind]] +
                                         inversed_edge_offset[neighbor_info_ind]],
                      delta);
            atomicSub(&residual_capacity[neighbor_info_ind], delta);
            atomicAdd(&excess_func[neighbors[neighbor_info_ind]],
                      static_cast<qmcp::QuasiMcpCudaMaxFlowSolver::Excess>(delta));
            atomicSub(&excess_func[node],
                      static_cast<qmcp::QuasiMcpCudaMaxFlowSolver::Excess>(delta));
        }
    }
}

void qmcp::QuasiMcpCudaMaxFlowSolver::add_edge(
    std::vector<std::vector<Node>>& neighbors_dict,
    std::vector<std::vector<EdgeDirection>>& edge_dir_dict,
    std::vector<std::vector<Capacity>>& residual_capacity_dict,
    std::vector<std::vector<uint32_t>>& inversed_edge_offset_dict, Node start, Node end,
    Capacity capacity) {
    size_t start_info_size = neighbors_dict[start].size();
    size_t end_info_size = neighbors_dict[end].size();

    neighbors_dict[start].push_back(end);
    edge_dir_dict[start].push_back(EdgeDirection::Forward);
    residual_capacity_dict[start].push_back(capacity);
    inversed_edge_offset_dict[start].push_back(end_info_size);

    neighbors_dict[end].push_back(start);
    edge_dir_dict[end].push_back(EdgeDirection::Backward);
    residual_capacity_dict[end].push_back(0);
    inversed_edge_offset_dict[end].push_back(start_info_size);
}

void qmcp::QuasiMcpCudaMaxFlowSolver::global_relabel() {
    uint32_t nodes_count = label_func_.size();
    // Step 1: Violation-Cancelation, fix all of the invalid residual edges
    for (Node u = 0; u < nodes_count - 2; ++u) {
        for (NeighborInfoIndex i = neighbors_start_ind_[u]; i <= neighbors_end_ind_[u]; ++i) {
            Node v = neighbors_[i];
            if (label_func_[u] <= label_func_[v] + 1 || residual_capacity_[i] == 0) {
                // Skip valid residual arcs
                continue;
            }

            uint32_t delta = static_cast<uint32_t>(excess_func_[u]);  // excess_func[node] > 0
            if (delta > residual_capacity_[i]) {
                delta = residual_capacity_[i];
            }

            excess_func_[u] -= static_cast<Excess>(delta);
            excess_func_[v] += static_cast<Excess>(delta);
            residual_capacity_[neighbors_start_ind_[v] + inversed_edge_offset_[i]] += delta;
            residual_capacity_[i] -= delta;
        }
    }

    Node sink = nodes_count - 1;

    // Step 2: Perform BFS backwards from the sink and assign each node with the BFS tree level
    typedef uint32_t BFSLevel;
    std::list<std::pair<Node, BFSLevel>> nodes;
    is_visited_.resize(nodes_count, false);

    nodes.push_back(std::make_pair(sink, 0));
    is_visited_[sink] = true;
    while (!nodes.empty()) {
        Node curr_node = nodes.front().first;
        BFSLevel curr_level = nodes.front().second;
        nodes.pop_front();

        label_func_[curr_node] = curr_level;

        for (NeighborInfoIndex i = neighbors_start_ind_[curr_node];
             i <= neighbors_end_ind_[curr_node]; ++i) {
            Node neighbor = neighbors_[i];
            // Capacity of neighbor->curr_node edge
            Capacity edge_cap =
                residual_capacity_[neighbors_start_ind_[neighbor] + inversed_edge_offset_[i]];

            if (is_visited_[neighbor] || edge_cap == 0) {
                continue;
            }

            is_visited_[neighbor] = true;
            nodes.push_back(std::make_pair(neighbor, curr_level + 1));
        }
    }
}

void qmcp::QuasiMcpCudaMaxFlowSolver::create_graph(const bam_api::SOAPairedReads& sequence,
                                                   uint32_t required_cover) {
    // Clear the graph data
    clear_graph();

    // Get genome length
    uint32_t n = sequence.ref_genome_length;

    // We first map 0, 1, 2, ... (n - 1) indices from the ref genome
    // to 1, 2, .... n. Then the graph is created in the following way, for read
    // (u, v), where u,v in (1, 2, ... n) create an edge (u - 1, v).
    //
    // The above procedure is equivalent of getting an original read (u, v),
    // where u,v in (0, 1, 2, ... n - 1) and createing an edge (u, v + 1).
    //
    // Additional nodes:
    // - source: n + 1,
    // - sink: n + 2,
    // - artificial node: 0.

    // Create source and sink
    Node source = n + 1;
    Node sink = n + 2;

    // Temporary dictionares with key=node
    std::vector<std::vector<Node>> neighbors_dict(n + 3);
    std::vector<std::vector<EdgeDirection>> edge_dir_dict(n + 3);
    std::vector<std::vector<Capacity>> residual_capacity_dict(n + 3);
    std::vector<std::vector<NeighborInfoIndex>> inversed_edge_offset_dict(n + 3);

    // Save the neighbor index for future export
    read_ind_to_neighbor_offset_.resize(sequence.get_reads_count());

    // Add edges that are corresponding to the reads
    for (bam_api::ReadIndex i = 0; i < sequence.get_reads_count(); ++i) {
        Node u = sequence.start_inds[i];
        Node v = sequence.end_inds[i] + 1;

        read_ind_to_neighbor_offset_[i] = neighbors_dict[u].size();

        // u --> v
        add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict, inversed_edge_offset_dict,
                 u, v, 1);
    }

    // Add returning edges
    for (Node i = 0; i < n; ++i) {
        // i + 1 --> i
        add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict, inversed_edge_offset_dict,
                 i + 1, i, std::numeric_limits<Capacity>::max());
    }

    // Create demand func basing on the required cover
    std::vector<int32_t> demand_func(n + 1, 0);
    for (bam_api::Index i = 0; i < n; ++i) {
        demand_func[i] = static_cast<int32_t>(std::min(max_coverage_[i], required_cover)) -
                         static_cast<int32_t>(std::min(max_coverage_[i + 1], required_cover));
    }
    demand_func[n] = static_cast<int32_t>(std::min(max_coverage_[n], required_cover));

    // Add edges for sink and source in order to simulate a circulation
    for (Node i = 0; i <= n; ++i) {
        if (demand_func[i] > 0) {
            // i --> sink
            add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict,
                     inversed_edge_offset_dict, i, sink, demand_func[i]);
        } else if (demand_func[i] < 0) {
            // source --> i
            add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict,
                     inversed_edge_offset_dict, source, i, -demand_func[i]);
        }
    }

    // Flatten the dictionaries and save them
    // Assumption: there are no isolated nodes
    neighbors_start_ind_.resize(n + 3, 0);
    neighbors_end_ind_.resize(n + 3, 0);

    uint32_t curr_ind = 0;

    for (bam_api::Index i = 0; i <= n + 2; ++i) {
        neighbors_start_ind_[i] = curr_ind;
        neighbors_end_ind_[i] = curr_ind + neighbors_dict[i].size() - 1;

        neighbors_.insert(neighbors_.end(), neighbors_dict[i].begin(), neighbors_dict[i].end());

        residual_capacity_.insert(residual_capacity_.end(), residual_capacity_dict[i].begin(),
                                  residual_capacity_dict[i].end());

        edge_dir_.insert(edge_dir_.end(), edge_dir_dict[i].begin(), edge_dir_dict[i].end());

        inversed_edge_offset_.insert(inversed_edge_offset_.end(),
                                     inversed_edge_offset_dict[i].begin(),
                                     inversed_edge_offset_dict[i].end());

        curr_ind += neighbors_dict[i].size();
    }

    // Prepare excess and label functions
    excess_func_.resize(n + 3, 0);
    label_func_.resize(n + 3, 0);
    label_func_[source] = n + 3;
}

qmcp::QuasiMcpCudaMaxFlowSolver::Excess qmcp::QuasiMcpCudaMaxFlowSolver::create_preflow() {
    // Get graph node count
    uint32_t n = label_func_.size();
    Node source = n - 2;
    Excess total_excess = 0;

    // Create preflow: saturate all edges coming out of the source
    for (uint32_t i = neighbors_start_ind_[source]; i <= neighbors_end_ind_[source]; ++i) {
        // We are in the source so every edge has forward direction
        // and checking the edge_dir is not requred

        Node curr_neighbor = neighbors_[i];
        Capacity curr_edge_capacity = residual_capacity_[i];

        // Find the inversed edge index
        NeighborInfoIndex inversed_i =
            neighbors_start_ind_[curr_neighbor] + inversed_edge_offset_[i];

        // Saturate the inversed edge
        residual_capacity_[inversed_i] = curr_edge_capacity;
        residual_capacity_[i] = 0;

        // Update the excess function
        excess_func_[curr_neighbor] = static_cast<Excess>(curr_edge_capacity);
        excess_func_[source] -= static_cast<Excess>(curr_edge_capacity);

        total_excess += static_cast<Excess>(curr_edge_capacity);
    }

    return total_excess;
}

void qmcp::QuasiMcpCudaMaxFlowSolver::clear_graph() {
    // From https://en.cppreference.com/w/cpp/container/vector/clear,
    // the capacity is unchanged, so it is beneficial to do clears
    // instead of creating new vectors every time, since we
    // may minimize the number of allocations between multiple calls
    // of the same instance of the solver.
    excess_func_.clear();
    label_func_.clear();
    neighbors_.clear();
    neighbors_start_ind_.clear();
    neighbors_end_ind_.clear();
    read_ind_to_neighbor_offset_.clear();
    inversed_edge_offset_.clear();
    edge_dir_.clear();
    residual_capacity_.clear();
    max_coverage_.clear();
}

void qmcp::QuasiMcpCudaMaxFlowSolver::set_block_size(uint32_t block_size) {
    block_size_ = block_size;
}

void qmcp::QuasiMcpCudaMaxFlowSolver::set_kernel_cycles(uint32_t kernel_cycles) {
    kernel_cycles_ = kernel_cycles;
}

std::unique_ptr<qmcp::Solution> qmcp::QuasiMcpCudaMaxFlowSolver::solve(uint32_t required_cover,
                                                                       bam_api::BamApi& bam_api) {
    input_sequence_ = bam_api.get_paired_reads_soa();

    // Create max coverage function
    max_coverage_.resize(input_sequence_.ref_genome_length + 1, 0);
    for (bam_api::ReadIndex i = 0; i < input_sequence_.get_reads_count(); ++i) {
        for (bam_api::Index j = input_sequence_.start_inds[i]; j <= input_sequence_.end_inds[i];
             ++j) {
            ++max_coverage_[j + 1];
        }
    }

    create_graph(input_sequence_, required_cover);

    uint32_t num_blocks = (excess_func_.size() + block_size_ - 1) / block_size_;

    Node sink = excess_func_.size() - 1;
    Node source = excess_func_.size() - 2;

    Excess excess_left = create_preflow();

    // Malloc the CUDA memory
    auto* dev_label_func = cuda::malloc<Label>(label_func_.size());
    auto* dev_excess_func = cuda::malloc<Excess>(excess_func_.size());
    auto* dev_neighbors_start_ind = cuda::malloc<NeighborInfoIndex>(neighbors_start_ind_.size());
    auto* dev_neighbors_end_ind = cuda::malloc<NeighborInfoIndex>(neighbors_end_ind_.size());
    auto* dev_neighbors = cuda::malloc<Node>(neighbors_.size());
    auto* dev_residual_capacity = cuda::malloc<Capacity>(residual_capacity_.size());
    auto* dev_inversed_edge_offset = cuda::malloc<NeighborInfoIndex>(inversed_edge_offset_.size());
    auto* dev_edge_dir = cuda::malloc<EdgeDirection>(edge_dir_.size());

    // Copy the arrays to device
    cuda::memcpy_host_dev<EdgeDirection>(dev_edge_dir, edge_dir_.data(), edge_dir_.size());
    cuda::memcpy_host_dev<NeighborInfoIndex>(dev_inversed_edge_offset, inversed_edge_offset_.data(),
                                             inversed_edge_offset_.size());
    cuda::memcpy_host_dev<NeighborInfoIndex>(dev_neighbors_start_ind, neighbors_start_ind_.data(),
                                             neighbors_start_ind_.size());
    cuda::memcpy_host_dev<NeighborInfoIndex>(dev_neighbors_end_ind, neighbors_end_ind_.data(),
                                             neighbors_end_ind_.size());
    cuda::memcpy_host_dev<Node>(dev_neighbors, neighbors_.data(), neighbors_.size());
    cuda::memcpy_host_dev<Capacity>(dev_residual_capacity, residual_capacity_.data(),
                                    residual_capacity_.size());
    cuda::memcpy_host_dev<Excess>(dev_excess_func, excess_func_.data(), excess_func_.size());
    cuda::memcpy_host_dev<Label>(dev_label_func, label_func_.data(), label_func_.size());

    // 1. Do the push-relabel with global-relabel heuristic
    while (excess_left > kUseGlobalRelabelMin) {
        // Call push-relabel GPU kernel
        push_relabel_kernel<<<num_blocks, block_size_>>>(
            kernel_cycles_, excess_func_.size(), dev_excess_func, dev_label_func,
            dev_residual_capacity, dev_neighbors_start_ind, dev_neighbors_end_ind, dev_neighbors,
            dev_inversed_edge_offset);

        // Copy the residual capacity, label and excess functions from device to
        // the host memory
        cuda::memcpy_dev_host<Capacity>(residual_capacity_.data(), dev_residual_capacity,
                                        residual_capacity_.size());
        cuda::memcpy_dev_host<Excess>(excess_func_.data(), dev_excess_func, excess_func_.size());
        cuda::memcpy_dev_host<Label>(label_func_.data(), dev_label_func, label_func_.size());

        // Call global-relabel on CPU
        global_relabel();

        excess_left = -excess_func_[source] - excess_func_[sink];

        cuda::memcpy_host_dev<Capacity>(dev_residual_capacity, residual_capacity_.data(),
                                        residual_capacity_.size());
        cuda::memcpy_host_dev<Excess>(dev_excess_func, excess_func_.data(), excess_func_.size());
        cuda::memcpy_host_dev<Label>(dev_label_func, label_func_.data(), label_func_.size());
    }

    // 2. Do raw push-relabel without any heuristics
    while (excess_left > 0) {
        // Call push-relabel GPU kernel
        push_relabel_kernel<<<num_blocks, block_size_>>>(
            kernel_cycles_, excess_func_.size(), dev_excess_func, dev_label_func,
            dev_residual_capacity, dev_neighbors_start_ind, dev_neighbors_end_ind, dev_neighbors,
            dev_inversed_edge_offset);

        Excess sink_excess = 0;
        Excess source_excess = 0;
        cuda::memcpy_dev_host<Excess>(&sink_excess, dev_excess_func + excess_func_.size() - 1, 1);
        cuda::memcpy_dev_host<Excess>(&source_excess, dev_excess_func + excess_func_.size() - 2, 1);

        excess_left = -source_excess - sink_excess;
    }

    cuda::memcpy_dev_host<Capacity>(residual_capacity_.data(), dev_residual_capacity,
                                    residual_capacity_.size());
    cuda::memcpy_dev_host<Excess>(excess_func_.data(), dev_excess_func, excess_func_.size());
    cuda::memcpy_dev_host<Label>(label_func_.data(), dev_label_func, label_func_.size());

    cuda::free(dev_label_func);
    cuda::free(dev_excess_func);
    cuda::free(dev_neighbors_start_ind);
    cuda::free(dev_inversed_edge_offset);
    cuda::free(dev_neighbors_end_ind);
    cuda::free(dev_neighbors);
    cuda::free(dev_residual_capacity);
    cuda::free(dev_edge_dir);

    std::unique_ptr<Solution> output = std::make_unique<Solution>();

    // Create output data
    for (bam_api::ReadIndex i = 0; i < read_ind_to_neighbor_offset_.size(); ++i) {
        Node u = input_sequence_.start_inds[i];
        Capacity cap =
            residual_capacity_[neighbors_start_ind_[u] + read_ind_to_neighbor_offset_[i]];

        if (cap == 0) {
            output->push_back(i);
        }
    }

    return output;
}

#endif