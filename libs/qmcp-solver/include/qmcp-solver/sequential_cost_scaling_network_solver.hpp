#ifndef QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP
#define QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP

#include <ortools/graph/min_cost_flow.h>

#include <vector>

#include "bam-api/bam_paired_reads.hpp"
#include "solver.hpp"

namespace qmcp {
class SequentialCostScalingNetworkSolver : public Solver {
   public:
    using Solver::Solver;
    std::vector<bam_api::BAMReadId> solve(uint32_t max_coverage) override;

   private:
    static void create_network_flow_graph(operations_research::SimpleMinCostFlow& min_cost_flow,
                                          const bam_api::AOSPairedReads& sequence, unsigned int M);
    static std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence,
                                                   unsigned int M);
    static std::vector<bam_api::ReadIndex> obtain_sequence(
        const bam_api::AOSPairedReads& sequence,
        const operations_research::SimpleMinCostFlow& min_cost_flow);

    bam_api::AOSPairedReads input_sequence_;
};
}  // namespace qmcp

#endif
