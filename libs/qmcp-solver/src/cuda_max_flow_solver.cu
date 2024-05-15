#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/check_cuda_error.hpp"
#include "qmcp-solver/cuda_max_flow_solver.hpp"

__global__ void push_relabel_kernel(int* data) {
    // TODO(billyk):
}

__host__ void global_relabel(int* data) {
    // TODO(billyk):
}

qmcp::CudaMaxFlowSolver::CudaMaxFlowSolver()
    : input_sequence_(), is_data_loaded_(false) {}

qmcp::CudaMaxFlowSolver::CudaMaxFlowSolver(
    const std::filesystem::path& filepath)
    : input_sequence_(), is_data_loaded_(false) {
    import_data(filepath);
}

void qmcp::CudaMaxFlowSolver::clear_data() {
    excess_func_.clear();
    label_func_.clear();
    neighbors_.clear();
    neighbors_start_ind_.clear();
    neighbors_end_ind_.clear();
    residual_capacity_.clear();
}

void qmcp::CudaMaxFlowSolver::import_data(
    const std::filesystem::path& filepath) {
    input_sequence_ = bam_api::BamApi::read_bam_soa(filepath);
    init_data(input_sequence_);
    is_data_loaded_ = true;
}

void qmcp::CudaMaxFlowSolver::init_data(const bam_api::SOAPairedReads& sequence,
                                        uint32_t required_cover) {
    // Clear the graph data
    clear_data();

    uint32_t n = sequence.ref_genome_length;

    // The nodes 0, 1, 2, ... (n - 1) are first mapped to the
    // indices 1, 2, .... n in the ref genome.
    //
    // Then the graph is created in the following way, for read (u, v), where
    // u,v in (1, 2, ... n) create an edge (u - 1, v).
    //
    // So the above procedure is equivalent of getting an original read (u, v),
    // where u,v in (0, 1, 2, ... n - 1) and createing an edge (u, v + 1).
    //
    // Additional nodes:
    // - source: n + 1,
    // - sink: n + 2,
    // - artificial node: 0.

    Node s = n + 1;
    Node t = n + 1;

    // Create coverage function
    std::vector<uint32_t> max_coverage(n + 1, 0);
    for (bam_api::ReadIndex i = 0; i < sequence.end_inds.size(); ++i) {
        for (bam_api::Index j = sequence.start_inds[i];
             j <= sequence.end_inds[i]; ++i) {
            ++max_coverage[j + 1];
        }
    }

    // Temporary dictionares with key=node
    std::vector<std::vector<Node>> neighbors_dict(n + 3);
    std::vector<std::vector<bool>> is_forward_dict(n + 3);
    std::vector<std::vector<Capacity>> residual_capacity_dict(n + 3);

    // Add edges that are corresponding to the reads
    for (bam_api::ReadIndex i = 0; i < sequence.end_inds.size(); ++i) {
        Node u = sequence.start_inds[i];
        Node v = sequence.end_inds[i] + 1;

        // u --> v
        neighbors_dict[u].push_back(v);
        is_forward_dict[u].push_back(true);
        residual_capacity_dict[u].push_back(1);

        neighbors_dict[v].push_back(u);
        is_forward_dict[v].push_back(false);
        residual_capacity_dict[v].push_back(0);
    }

    // Add returning edges
    for (Node i = 0; i < n; ++i) {
        // i <-- i + 1
        neighbors_dict[i + 1].push_back(i);
        is_forward_dict[i + 1].push_back(true);
        residual_capacity_dict[i + 1].push_back(
            std::numeric_limits<Capacity>::max());

        neighbors_dict[i].push_back(i + 1);
        is_forward_dict[i].push_back(false);
        residual_capacity_dict[i].push_back(0);
    }

    // Create demand func basing on the required cover
    std::vector<uint32_t> demand_func(n + 1, 0);
    for (bam_api::Index i = 0; i < n; ++i) {
        demand_func[i] = std::min(max_coverage[i + 1], required_cover) -
                         std::min(max_coverage[i], required_cover);
    }

    // Add edges for sink and source in order to simulate a circulation
    for (Node i = 0; i <= n; ++i) {
        if (demand_func[i] > 0) {
            // Sink: i --> t
            neighbors_dict[i].push_back(t);
            is_forward_dict[i].push_back(true);
            residual_capacity_dict[i].push_back(demand_func[i]);

            neighbors_dict[t].push_back(i);
            is_forward_dict[t].push_back(false);
            residual_capacity_dict[t].push_back(0);
        } else if (demand_func[i] < 0) {
            // Source: s --> i
            neighbors_dict[s].push_back(i);
            is_forward_dict[s].push_back(true);
            residual_capacity_dict[s].push_back(demand_func[i]);

            neighbors_dict[i].push_back(s);
            is_forward_dict[i].push_back(false);
            residual_capacity_dict[i].push_back(0);
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

        is_forward_.insert(is_forward_.end(), is_forward_dict[i].begin(),
                           is_forward_dict[i].end());

        curr_ind += neighbors_dict[i].size() + 1;
    }

    // Prepare excess and label functions sizes
    excess_func_.resize(n + 3, 0);
    label_func_.resize(n + 3, 0);
}

void qmcp::CudaMaxFlowSolver::solve(uint32_t required_cover) {
    if (!is_data_loaded_) {
        std::cerr << "Couldn't run solver: data has not been loaded.\n";
        std::exit(EXIT_FAILURE);
    }

    // Prepare graph
    init_data(input_sequence_, required_cover);

    // Malloc and initialize CUDA memory
    int32_t* dev_label_func = nullptr;
    int32_t* dev_excess_func = nullptr;

    CHECK_CUDA_ERROR(
        cudaMalloc(reinterpret_cast<void**>(&dev_label_func),
                   input_sequence_.ref_genome_length * sizeof(int32_t)));
    CHECK_CUDA_ERROR(cudaMemcpy(dev_excess_func, label_func_.data(),
                                label_func_.size() * sizeof(int32_t),
                                cudaMemcpyKind::cudaMemcpyHostToDevice));
    CHECK_CUDA_ERROR(
        cudaMalloc(reinterpret_cast<void**>(&dev_excess_func),
                   input_sequence_.ref_genome_length * sizeof(int32_t)));
    CHECK_CUDA_ERROR(cudaMemcpy(dev_excess_func, label_func_.data(),
                                label_func_.size() * sizeof(int32_t),
                                cudaMemcpyKind::cudaMemcpyHostToDevice));

    CHECK_CUDA_ERROR(cudaFree(dev_label_func));
    CHECK_CUDA_ERROR(cudaFree(dev_excess_func));
}
