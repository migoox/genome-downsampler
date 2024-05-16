#ifndef QMCP_CUDA_MAX_FLOW_SOLVER_HPP
#define QMCP_CUDA_MAX_FLOW_SOLVER_HPP

#include <cstdint>
#include <filesystem>
#include <memory>

#include "bam-api/bam_paired_reads.hpp"
#include "solver.hpp"

namespace qmcp {

class CudaMaxFlowSolver : public Solver {
   public:
    typedef uint32_t Node;
    typedef uint32_t Capacity;
    typedef int32_t Excess;

    CudaMaxFlowSolver();
    enum class EdgeDirection : uint8_t { Forward, Backward };
    explicit CudaMaxFlowSolver(const std::filesystem::path& filepath);

    void import_data(const std::filesystem::path& filepath);
    void solve(uint32_t required_cover) override;
    void export_data(const std::filesystem::path& filepath);

   private:
    void clear_graph();

    // Calling this function is valid only after the import_data is called
    void create_graph(const bam_api::SOAPairedReads& sequence,
                      uint32_t required_cover);

    // Helper functions
    static void add_edge(
        std::vector<std::vector<Node>>& neighbors_dict,
        std::vector<std::vector<EdgeDirection>>& edge_dir_dict,
        std::vector<std::vector<Capacity>>& residual_capacity_dict,
        std::vector<std::vector<uint32_t>>& inversed_edge_ind_dict, Node start,
        Node end, Capacity capacity);

    // This function is responsible for first step of push-relabel algorithm
    void create_preflow();

    bool is_data_loaded_;

    bam_api::SOAPairedReads input_sequence_;
    std::vector<uint32_t> max_coverage_;

    // Graph data
    std::vector<Excess> excess_func_;
    std::vector<uint32_t> label_func_;
    std::vector<uint32_t> neighbors_start_ind_;
    std::vector<uint32_t> neighbors_end_ind_;
    std::vector<Node> neighbors_;
    std::vector<uint32_t> inversed_edge_ind_;
    std::vector<Capacity> residual_capacity_;
    std::vector<EdgeDirection> edge_dir_;
};

}  // namespace qmcp

#endif