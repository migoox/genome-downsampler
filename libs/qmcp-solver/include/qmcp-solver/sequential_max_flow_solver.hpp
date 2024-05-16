#ifndef QMCP_SEQUENTIAL_MAX_FLOW_SOLVER_HPP
#define QMCP_SEQUENTIAL_MAX_FLOW_SOLVER_HPP
#include <filesystem>
#include <iostream>
#include <vector>

#include "solver.hpp"
// #include <ortools/graph/min_cost_flow.h>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"

namespace qmcp {

class SequentialMaxFlowSolver : public Solver {
   public:
    SequentialMaxFlowSolver(unsigned int M,
                            const std::filesystem::path& filepath)
        : M_(M) {
        input_sequence_ = bam_api::BamApi::read_bam_aos(filepath);
        std::cout << "read bam\n";
    }
    void solve() override;
    std::vector<bam_api::ReadIndex> output_sequence();

   private:
    // static void
    // create_network_flow_graph(operations_research::SimpleMinCostFlow&
    // min_cost_flow, const bam_api::AOSPairedReads& sequence, unsigned int M);
    bam_api::AOSPairedReads input_sequence_;
    std::vector<bam_api::ReadIndex> output_sequence_;
    unsigned int M_;
};

}  // namespace qmcp

#endif
