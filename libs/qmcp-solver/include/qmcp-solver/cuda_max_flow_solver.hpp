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

    CudaMaxFlowSolver();
    explicit CudaMaxFlowSolver(const std::filesystem::path& filepath);

    void import_data(const std::filesystem::path& filepath);
    void solve(uint32_t required_cover) override;
    void export_data(const std::filesystem::path& filepath);

   private:
    void clear_data();
    void init_data(const bam_api::SOAPairedReads& sequence,
                   uint32_t required_cover);

    bool is_data_loaded_;

    bam_api::SOAPairedReads input_sequence_;

    // Graph data
    std::vector<int32_t> excess_func_;
    std::vector<uint32_t> label_func_;
    std::vector<uint32_t> neighbors_start_ind_;
    std::vector<uint32_t> neighbors_end_ind_;
    std::vector<Node> neighbors_;
    std::vector<Capacity> residual_capacity_;
    std::vector<bool> is_forward_;
};

}  // namespace qmcp

#endif