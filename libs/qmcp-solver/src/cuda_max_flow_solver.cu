#include <stdint.h>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <list>
#include <optional>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/cuda_helpers.cuh"
#include "qmcp-solver/cuda_max_flow_solver.hpp"

__global__ void push_relabel_kernel(
    uint32_t kernel_cycles, uint32_t nodes_count, qmcp::CudaMaxFlowSolver::Excess* excess_func,
    qmcp::CudaMaxFlowSolver::Label* label_func,
    qmcp::CudaMaxFlowSolver::Capacity* residual_capacity,
    const qmcp::CudaMaxFlowSolver::NeighborInfoIndex* neighbors_start_ind,
    const qmcp::CudaMaxFlowSolver::NeighborInfoIndex* neighbors_end_ind,
    const qmcp::CudaMaxFlowSolver::Node* neighbors,
    const qmcp::CudaMaxFlowSolver::NeighborInfoIndex* inversed_edge_ind) {
    qmcp::CudaMaxFlowSolver::Node node = blockDim.x * blockIdx.x + threadIdx.x;
    if (node > nodes_count) {
        return;
    }

    for (uint32_t cycle = 0; cycle < kernel_cycles; ++cycle) {
        // Ensure that current node is active
        if (excess_func[node] <= 0 || label_func[node] >= nodes_count) {
            continue;
        }

        // Find neighbor with minimum label
        qmcp::CudaMaxFlowSolver::Label min_label = UINT32_MAX;
        qmcp::CudaMaxFlowSolver::NeighborInfoIndex neighbor_info_ind = 0;

        for (qmcp::CudaMaxFlowSolver::NeighborInfoIndex i = neighbors_start_ind[node];
             i <= neighbors_end_ind[node]; ++i) {
            qmcp::CudaMaxFlowSolver::Node curr_neighbor = neighbors[i];
            qmcp::CudaMaxFlowSolver::Label neighbor_label = label_func[curr_neighbor];

            if (min_label > neighbor_label) {
                min_label = neighbor_label;
                neighbor_info_ind = i;
            }
        }

        // Do push-relabel
        if (min_label >= label_func[node]) {
            // Preform relabel operation
            label_func[node] = min_label + 1;
        } else {
            // Perform push operation
            int32_t delta =
                min(static_cast<int32_t>(residual_capacity[neighbor_info_ind]), excess_func[node]);

            atomicAdd(&residual_capacity[inversed_edge_ind[neighbor_info_ind]], delta);
            atomicSub(&residual_capacity[neighbor_info_ind], delta);

            atomicAdd(&excess_func[neighbors[neighbor_info_ind]], delta);
            atomicSub(&excess_func[node], delta);
        }
    }
}

qmcp::CudaMaxFlowSolver::CudaMaxFlowSolver() : is_data_loaded_(false) {}

qmcp::CudaMaxFlowSolver::CudaMaxFlowSolver(const std::filesystem::path& filepath,
                                           uint32_t min_seq_length, uint32_t min_seq_mapq)
    : is_data_loaded_(false) {
    import_reads(filepath, min_seq_length, min_seq_mapq);
}

void qmcp::CudaMaxFlowSolver::import_reads(const std::filesystem::path& filepath,
                                           uint32_t min_seq_length, uint32_t min_seq_mapq) {
    input_filepath_ = filepath;
    input_sequence_ = bam_api::BamApi::read_bam_soa(filepath, min_seq_length, min_seq_mapq);

    // Create max coverage function
    max_coverage_.resize(input_sequence_.ref_genome_length + 1, 0);
    for (bam_api::ReadIndex i = 0; i < input_sequence_.end_inds.size(); ++i) {
        for (bam_api::Index j = input_sequence_.start_inds[i]; j <= input_sequence_.end_inds[i];
             ++i) {
            ++max_coverage_[j + 1];
        }
    }

    is_data_loaded_ = true;
}
void qmcp::CudaMaxFlowSolver::add_edge(std::vector<std::vector<Node>>& neighbors_dict,
                                       std::vector<std::vector<EdgeDirection>>& edge_dir_dict,
                                       std::vector<std::vector<Capacity>>& residual_capacity_dict,
                                       std::vector<std::vector<uint32_t>>& inversed_edge_ind_dict,
                                       Node start, Node end, Capacity capacity) {
    size_t start_info_size = neighbors_dict[start].size();
    size_t end_info_size = neighbors_dict[end].size();

    neighbors_dict[start].push_back(end);
    edge_dir_dict[start].push_back(EdgeDirection::Forward);
    residual_capacity_dict[start].push_back(capacity);
    inversed_edge_ind_dict[start].push_back(end_info_size);

    neighbors_dict[end].push_back(start);
    edge_dir_dict[end].push_back(EdgeDirection::Backward);
    residual_capacity_dict[end].push_back(0);
    inversed_edge_ind_dict[end].push_back(start_info_size);
}

void qmcp::CudaMaxFlowSolver::global_relabel(Excess& excess_total) {
    uint32_t nodes_count = label_func_.size();
    for (Node u = 0; u < nodes_count; ++u) {
        for (NeighborInfoIndex i = neighbors_start_ind_[u]; i <= neighbors_end_ind_[u]; ++i) {
            if (edge_dir_[i] == EdgeDirection::Backward) {
                continue;
            }

            Node v = neighbors_[i];

            excess_func_[u] = excess_func_[u] - static_cast<Excess>(residual_capacity_[i]);
            excess_func_[v] = excess_func_[v] + static_cast<Excess>(residual_capacity_[i]);
            residual_capacity_[inversed_edge_ind_[i]] += residual_capacity_[i];
            residual_capacity_[i] = 0;
        }
    }

    Node sink = nodes_count - 1;

    // Perform BFS backwards from the sink to do the global relabel
    std::list<Node> nodes;
    std::vector<Node> not_relabeled_nodes;
    std::vector<bool> is_visited(nodes_count, false);
    int32_t bfs_level = 0;

    nodes.push_back(sink);
    while (!nodes.empty()) {
        Node curr_node = nodes.front();
        nodes.pop_front();

        if (bfs_level == label_func_[curr_node] && !is_markded_[curr_node]) {
            not_relabeled_nodes.push_back(curr_node);
        }
        label_func_[curr_node] = bfs_level++;
        is_visited[curr_node] = true;

        // Add new vertices
        for (NeighborInfoIndex i = neighbors_start_ind_[curr_node];
             i < neighbors_end_ind_[curr_node]; ++i) {
            Node neighbor = neighbors_[i];
            if (is_visited[neighbor] || edge_dir_[i] == EdgeDirection::Forward) {
                continue;
            }

            nodes.push_back(neighbor);
        }
    }

    // Update the excess total
    for (Node node : not_relabeled_nodes) {
        is_markded_[node] = true;
        excess_total = excess_total - excess_func_[node];
    }
}

void qmcp::CudaMaxFlowSolver::create_graph(const bam_api::SOAPairedReads& sequence,
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
    std::vector<std::vector<NeighborInfoIndex>> inversed_edge_ind_dict(n + 3);

    // Save the neighbor index for future export
    read_ind_to_neighbor_ind_.resize(sequence.end_inds.size());

    // Add edges that are corresponding to the reads
    for (bam_api::ReadIndex i = 0; i < sequence.end_inds.size(); ++i) {
        Node u = sequence.start_inds[i];
        Node v = sequence.end_inds[i] + 1;

        read_ind_to_neighbor_ind_[i] = neighbors_dict[u].size();

        // u --> v
        add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict, inversed_edge_ind_dict, u,
                 v, 1);
    }

    // Add returning edges
    for (Node i = 0; i < n; ++i) {
        // i + 1 --> i
        add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict, inversed_edge_ind_dict,
                 i + 1, i, std::numeric_limits<Capacity>::max());
    }

    // Create demand func basing on the required cover
    std::vector<int32_t> demand_func(n + 1, 0);
    for (bam_api::Index i = 0; i < n; ++i) {
        demand_func[i] = static_cast<int32_t>(std::min(max_coverage_[i + 1], required_cover)) -
                         static_cast<int32_t>(std::min(max_coverage_[i], required_cover));
    }

    // Add edges for sink and source in order to simulate a circulation
    for (Node i = 0; i <= n; ++i) {
        if (demand_func[i] > 0) {
            // i --> sink
            add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict, inversed_edge_ind_dict,
                     i, sink, demand_func[i]);
        } else if (demand_func[i] < 0) {
            // source --> i
            add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict, inversed_edge_ind_dict,
                     source, i, -demand_func[i]);
        }
    }

    // Flatten the dictionaries and save them
    // Assumption: there are no isolated nodes
    neighbors_end_ind_.resize(n + 3, 0);
    neighbors_start_ind_.resize(n + 3, 0);

    uint32_t curr_ind = 0;

    for (bam_api::Index i = 0; i <= n + 2; ++i) {
        neighbors_start_ind_.push_back(curr_ind);
        neighbors_end_ind_.push_back(curr_ind + neighbors_dict[i].size());

        neighbors_.insert(neighbors_.end(), neighbors_dict[i].begin(), neighbors_dict[i].end());

        residual_capacity_.insert(residual_capacity_.end(), residual_capacity_dict[i].begin(),
                                  residual_capacity_dict[i].end());

        edge_dir_.insert(edge_dir_.end(), edge_dir_dict[i].begin(), edge_dir_dict[i].end());

        inversed_edge_ind_.insert(inversed_edge_ind_.end(), inversed_edge_ind_dict[i].begin(),
                                  inversed_edge_ind_dict[i].end());

        curr_ind += neighbors_dict[i].size() + 1;
    }

    // Prepare excess and label functions
    excess_func_.resize(n + 3, 0);
    label_func_.resize(n + 3, 0);
    label_func_[source] = n + 3;
    create_preflow();

    is_markded_.resize(n + 3, false);
}

void qmcp::CudaMaxFlowSolver::create_preflow() {
    // Get graph node count
    uint32_t n = label_func_.size();
    Node source = n - 2;

    // Create preflow: saturate all edges coming out of the source
    for (uint32_t i = neighbors_start_ind_[source]; i <= neighbors_end_ind_[source]; ++i) {
        // We are in the source so every edge has forward direction
        // and checking the edge_dir is not requred

        // Get current neighbor
        Node curr_neighbor = neighbors_[i];
        Capacity curr_edge_capacity = residual_capacity_[i];

        // Get the inversed edge location
        uint32_t inversed_i = neighbors_start_ind_[curr_neighbor] + inversed_edge_ind_[i];

        // Saturate the edge
        residual_capacity_[inversed_i] = curr_edge_capacity;
        residual_capacity_[i] = 0;

        // Update the excess function
        excess_func_[curr_neighbor] = static_cast<Excess>(curr_edge_capacity);
        excess_func_[source] -= static_cast<Excess>(curr_edge_capacity);
    }
}

void qmcp::CudaMaxFlowSolver::clear_graph() {
    excess_func_.clear();
    label_func_.clear();
    neighbors_.clear();
    neighbors_start_ind_.clear();
    neighbors_end_ind_.clear();
    residual_capacity_.clear();
    max_coverage_.clear();
    output_.clear();
}

void qmcp::CudaMaxFlowSolver::export_reads(const std::filesystem::path& filepath) {
    bam_api::BamApi::write_bam(input_filepath_, filepath, output_);
}

void qmcp::CudaMaxFlowSolver::set_block_size(uint32_t block_size) { block_size_ = block_size; }

void qmcp::CudaMaxFlowSolver::set_kernel_cycles(uint32_t kernel_cycles) {
    kernel_cycles_ = kernel_cycles;
}

void qmcp::CudaMaxFlowSolver::solve(uint32_t required_cover) {
    if (!is_data_loaded_) {
        std::cerr << "Couldn't run solver: data has not been loaded.\n";
        std::exit(EXIT_FAILURE);
    }

    create_graph(input_sequence_, required_cover);

    // Malloc the CUDA memory
    auto* dev_label_func = cuda::malloc<Label>(label_func_.size());
    auto* dev_excess_func = cuda::malloc<Excess>(excess_func_.size());
    auto* dev_neighbors_start_ind = cuda::malloc<NeighborInfoIndex>(neighbors_start_ind_.size());
    auto* dev_neighbors_end_ind = cuda::malloc<NeighborInfoIndex>(neighbors_end_ind_.size());
    auto* dev_neighbors = cuda::malloc<Node>(neighbors_.size());
    auto* dev_residual_capacity = cuda::malloc<Capacity>(residual_capacity_.size());
    auto* dev_inversed_edge_ind = cuda::malloc<NeighborInfoIndex>(inversed_edge_ind_.size());
    auto* dev_edge_dir = cuda::malloc<EdgeDirection>(edge_dir_.size());

    // Copy the constant arrays to device
    cuda::memcpy<EdgeDirection>(dev_edge_dir, edge_dir_.data(), edge_dir_.size(),
                                cudaMemcpyHostToDevice);
    cuda::memcpy<NeighborInfoIndex>(dev_inversed_edge_ind, inversed_edge_ind_.data(),
                                    inversed_edge_ind_.size(), cudaMemcpyHostToDevice);
    cuda::memcpy<NeighborInfoIndex>(dev_neighbors_start_ind, neighbors_start_ind_.data(),
                                    neighbors_start_ind_.size(), cudaMemcpyHostToDevice);
    cuda::memcpy<NeighborInfoIndex>(dev_neighbors_end_ind, neighbors_end_ind_.data(),
                                    neighbors_end_ind_.size(), cudaMemcpyHostToDevice);
    cuda::memcpy<Node>(dev_neighbors, neighbors_.data(), neighbors_.size(), cudaMemcpyHostToDevice);

    // Copy the excess function and residual capacity to the device memory
    cuda::memcpy<Capacity>(dev_residual_capacity, residual_capacity_.data(),
                           residual_capacity_.size(), cudaMemcpyHostToDevice);
    cuda::memcpy<Excess>(dev_excess_func, excess_func_.data(), excess_func_.size(),
                         cudaMemcpyHostToDevice);

    uint32_t num_blocks = (excess_func_.size() + block_size_ - 1) / block_size_;

    Node sink = excess_func_.size() - 1;
    Node source = excess_func_.size() - 2;

    Excess total_excess = 0;
    while (excess_func_[source] + excess_func_[sink] < total_excess) {
        // Copy the label function to the device memory
        cuda::memcpy<Label>(dev_label_func, label_func_.data(),
                            label_func_.size() * sizeof(uint32_t), cudaMemcpyHostToDevice);

        // Call push-relabel GPU kernel
        push_relabel_kernel<<<num_blocks, block_size_>>>(
            kernel_cycles_, excess_func_.size(), dev_excess_func, dev_label_func,
            dev_residual_capacity, dev_neighbors_start_ind, dev_neighbors_end_ind, dev_neighbors,
            dev_inversed_edge_ind);

        // Copy the residual capacity, label and excess functions from device to
        // the host memory
        cuda::memcpy<Capacity>(residual_capacity_.data(), dev_residual_capacity,
                               residual_capacity_.size(), cudaMemcpyDeviceToHost);
        cuda::memcpy<Excess>(excess_func_.data(), dev_excess_func, excess_func_.size(),
                             cudaMemcpyDeviceToHost);
        cuda::memcpy<Label>(label_func_.data(), dev_label_func, label_func_.size(),
                            cudaMemcpyHostToDevice);

        // Call global-relabel on CPU
        global_relabel(total_excess);
    }

    cuda::free(dev_label_func);
    cuda::free(dev_excess_func);
    cuda::free(dev_neighbors_start_ind);
    cuda::free(dev_inversed_edge_ind);
    cuda::free(dev_neighbors_end_ind);
    cuda::free(dev_neighbors);
    cuda::free(dev_residual_capacity);
    cuda::free(dev_edge_dir);

    // Create output data
    for (bam_api::ReadIndex i = 0; i < read_ind_to_neighbor_ind_.size(); ++i) {
        Node u = input_sequence_.start_inds[i];
        Capacity cap = residual_capacity_[neighbors_start_ind_[u] + read_ind_to_neighbor_ind_[i]];

        if (cap == 0) {
            output_.push_back(i);

            bam_api::Index pair_ind = input_sequence_.is_first_reads[i] ? (i + 1) : (i - 1);
            Node pair_u = input_sequence_.start_inds[pair_ind];

            // Check if the pair index also must be added
            Capacity pair_cap = residual_capacity_[neighbors_start_ind_[pair_u] +
                                                   read_ind_to_neighbor_ind_[pair_ind]];
            if (pair_cap != 0) {
                output_.push_back(pair_ind);
            }
        }
    }
}