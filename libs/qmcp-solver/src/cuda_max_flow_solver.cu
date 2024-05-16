#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/cuda_helpers.cuh"
#include "qmcp-solver/cuda_max_flow_solver.hpp"

__global__ void push_relabel_kernel(int* data) {
    // TODO(billyk):
}

__host__ void global_relabel(int* data) {
    // TODO(billyk):
}

qmcp::CudaMaxFlowSolver::CudaMaxFlowSolver() : is_data_loaded_(false) {}

qmcp::CudaMaxFlowSolver::CudaMaxFlowSolver(
    const std::filesystem::path& filepath)
    : is_data_loaded_(false) {
    import_data(filepath);
}

void qmcp::CudaMaxFlowSolver::import_data(
    const std::filesystem::path& filepath) {
    input_sequence_ = bam_api::BamApi::read_bam_soa(filepath);

    // Create max coverage function
    max_coverage_.resize(input_sequence_.ref_genome_length + 1, 0);
    for (bam_api::ReadIndex i = 0; i < input_sequence_.end_inds.size(); ++i) {
        for (bam_api::Index j = input_sequence_.start_inds[i];
             j <= input_sequence_.end_inds[i]; ++i) {
            ++max_coverage_[j + 1];
        }
    }

    is_data_loaded_ = true;
}
void qmcp::CudaMaxFlowSolver::add_edge(
    std::vector<std::vector<Node>>& neighbors_dict,
    std::vector<std::vector<EdgeDirection>>& edge_dir_dict,
    std::vector<std::vector<Capacity>>& residual_capacity_dict,
    std::vector<std::vector<uint32_t>>& inversed_edge_ind_dict, Node start,
    Node end, Capacity capacity) {
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

void qmcp::CudaMaxFlowSolver::create_graph(
    const bam_api::SOAPairedReads& sequence, uint32_t required_cover) {
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
    std::vector<std::vector<uint32_t>> inversed_edge_ind_dict(n + 3);

    // Add edges that are corresponding to the reads
    for (bam_api::ReadIndex i = 0; i < sequence.end_inds.size(); ++i) {
        Node u = sequence.start_inds[i];
        Node v = sequence.end_inds[i] + 1;

        // u --> v
        add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict,
                 inversed_edge_ind_dict, u, v, 1);
    }

    // Add returning edges
    for (Node i = 0; i < n; ++i) {
        // i + 1 --> i
        add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict,
                 inversed_edge_ind_dict, i + 1, i,
                 std::numeric_limits<Capacity>::max());
    }

    // Create demand func basing on the required cover
    std::vector<uint32_t> demand_func(n + 1, 0);
    for (bam_api::Index i = 0; i < n; ++i) {
        demand_func[i] = std::min(max_coverage_[i + 1], required_cover) -
                         std::min(max_coverage_[i], required_cover);
    }

    // Add edges for sink and source in order to simulate a circulation
    for (Node i = 0; i <= n; ++i) {
        if (demand_func[i] > 0) {
            // i --> sink
            add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict,
                     inversed_edge_ind_dict, i, sink, demand_func[i]);
        } else if (demand_func[i] < 0) {
            // source --> i
            add_edge(neighbors_dict, edge_dir_dict, residual_capacity_dict,
                     inversed_edge_ind_dict, source, i, -demand_func[i]);
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

        neighbors_.insert(neighbors_.end(), neighbors_dict[i].begin(),
                          neighbors_dict[i].end());

        residual_capacity_.insert(residual_capacity_.end(),
                                  residual_capacity_dict[i].begin(),
                                  residual_capacity_dict[i].end());

        edge_dir_.insert(edge_dir_.end(), edge_dir_dict[i].begin(),
                         edge_dir_dict[i].end());

        inversed_edge_ind_.insert(neighbors_end_ind_.end(),
                                  inversed_edge_ind_dict[i].begin(),
                                  inversed_edge_ind_dict[i].end());

        curr_ind += neighbors_dict[i].size() + 1;
    }

    // Prepare excess and label functions
    excess_func_.resize(n + 3, 0);
    label_func_.resize(n + 3, 0);
    label_func_[source] = n + 3;
    create_preflow();
}

void qmcp::CudaMaxFlowSolver::create_preflow() {
    // Get graph node count
    uint32_t n = label_func_.size();
    Node source = n - 2;

    // Create preflow: saturate all edges coming out of the source
    for (uint32_t i = neighbors_start_ind_[source];
         i <= neighbors_end_ind_[source]; ++i) {
        // We are in the source so every edge has forward direction
        // and checking the edge_dir is not requred

        // Get current neighbor
        Node curr_neighbor = neighbors_[i];
        Capacity curr_edge_capacity = residual_capacity_[i];

        // Get the inversed edge location
        uint32_t inversed_i =
            neighbors_start_ind_[curr_neighbor] + inversed_edge_ind_[i];

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
    auto* dev_neighbors_start_ind =
        cuda::malloc<NeighborInfoIndex>(neighbors_start_ind_.size());
    auto* dev_neighbors_end_ind =
        cuda::malloc<NeighborInfoIndex>(neighbors_end_ind_.size());
    auto* dev_neighbors = cuda::malloc<Node>(neighbors_.size());
    auto* dev_residual_capacity =
        cuda::malloc<Capacity>(residual_capacity_.size());
    auto* dev_inversed_edge_ind =
        cuda::malloc<NeighborInfoIndex>(inversed_edge_ind_.size());
    auto* dev_edge_dir = cuda::malloc<EdgeDirection>(edge_dir_.size());

    // Copy the constant arrays to device
    cuda::memcpy(dev_edge_dir, edge_dir_.data(), edge_dir_.size(),
                 cudaMemcpyHostToDevice);
    cuda::memcpy(dev_inversed_edge_ind, inversed_edge_ind_.data(),
                 inversed_edge_ind_.size(), cudaMemcpyHostToDevice);
    cuda::memcpy(dev_neighbors_start_ind, neighbors_start_ind_.data(),
                 neighbors_start_ind_.size(), cudaMemcpyHostToDevice);
    cuda::memcpy(dev_neighbors_end_ind, neighbors_end_ind_.data(),
                 neighbors_end_ind_.size(), cudaMemcpyHostToDevice);
    cuda::memcpy(dev_neighbors, neighbors_.data(), neighbors_.size(),
                 cudaMemcpyHostToDevice);

    // Copy the excess function and residual capacity
    cuda::memcpy(dev_residual_capacity, residual_capacity_.data(),
                 residual_capacity_.size(), cudaMemcpyHostToDevice);
    cuda::memcpy(dev_excess_func, excess_func_.data(), excess_func_.size(),
                 cudaMemcpyHostToDevice);

    Node sink = excess_func_.size() - 1;
    Node source = excess_func_.size() - 2;

    Excess total_excess = 0;
    while (excess_func_[source] + excess_func_[sink] < total_excess) {
        // Copy the label function to the device memory
        cuda::memcpy(dev_label_func, label_func_.data(),
                     label_func_.size() * sizeof(uint32_t),
                     cudaMemcpyHostToDevice);
        // Call push-relabel GPU kernel
        // TODO(billyk):

        // Copy the residual capacity, label and excess functions from device to
        // the host memory
        cuda::memcpy(residual_capacity_.data(), dev_residual_capacity,
                     residual_capacity_.size(), cudaMemcpyDeviceToHost);
        cuda::memcpy(excess_func_.data(), dev_excess_func, excess_func_.size(),
                     cudaMemcpyDeviceToHost);
        cuda::memcpy(label_func_.data(), dev_label_func, label_func_.size(),
                     cudaMemcpyHostToDevice);

        // Call global-relabel on CPU
        // TODO(billyk):
    }

    cuda::free(dev_label_func);
    cuda::free(dev_excess_func);
    cuda::free(dev_neighbors_start_ind);
    cuda::free(dev_inversed_edge_ind);
    cuda::free(dev_neighbors_end_ind);
    cuda::free(dev_neighbors);
    cuda::free(dev_residual_capacity);
    cuda::free(dev_edge_dir);
}
