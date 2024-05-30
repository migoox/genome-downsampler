#ifndef QMCP_SEQUENTIAL_MAX_FLOW_SOLVER_HPP
#define QMCP_SEQUENTIAL_MAX_FLOW_SOLVER_HPP
#include <ortools/graph/max_flow.h>

#include <vector>

#include "bam-api/bam_paired_reads.hpp"
#include "solver.hpp"

namespace qmcp {

class SequentialMaxFlowSolver : public Solver {
   public:
    using Solver::Solver;
    std::vector<bam_api::BAMReadId> solve(uint32_t max_coverage) override;

   private:
    static void create_network_flow_graph(operations_research::SimpleMaxFlow& max_flow,
                                          const bam_api::AOSPairedReads& sequence, unsigned int M);
    static std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence,
                                                   unsigned int M);
    static std::vector<bam_api::BAMReadId> obtain_sequence(
        const bam_api::AOSPairedReads& sequence, const operations_research::SimpleMaxFlow& max_flow);
    static std::vector<int> create_b_function(const bam_api::AOSPairedReads& sequence,
                                              unsigned int M);

    bam_api::AOSPairedReads input_sequence_;
};

}  // namespace qmcp

#endif
