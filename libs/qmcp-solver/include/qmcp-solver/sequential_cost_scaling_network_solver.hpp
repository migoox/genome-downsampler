#ifndef QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP
#define QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP
#include <filesystem>
#include <iostream>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
#include "solver.hpp"

namespace qmcp {

class SequentialCostScalingNetworkSolver : public Solver {
   public:
    explicit SequentialCostScalingNetworkSolver(
        const std::filesystem::path& filepath) {
        input_sequence_ = bam_api::BamApi::read_bam_aos(filepath);
    }

    void solve(uint32_t required_cover) override;

   private:
    // static void
    // create_network_flow_graph(operations_research::SimpleMinCostFlow&
    // min_cost_flow, const bam_api::AOSPairedReads& sequence, unsigned int M);
    bam_api::AOSPairedReads input_sequence_;
    bam_api::AOSPairedReads output_sequence_;
};

}  // namespace qmcp

#endif
