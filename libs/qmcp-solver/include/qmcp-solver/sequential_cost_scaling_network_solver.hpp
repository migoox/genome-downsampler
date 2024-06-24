#ifndef QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP
#define QMCP_SEQUENTIAL_COST_SCALING_NETWORK_SOLVER_HPP

#include <ortools/graph/min_cost_flow.h>

#include <cstdint>
#include <memory>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "solver.hpp"

namespace qmcp {
class SequentialCostScalingNetworkSolver : public Solver {
   public:
    std::unique_ptr<Solution> solve(uint32_t max_coverage, bam_api::BamApi& bam_api) override;
    bool uses_quality_of_reads() override { return true; }

   private:
    static std::vector<int> create_b_function(const bam_api::AOSPairedReads& sequence, uint32_t M);
    static void create_network_flow_graph(operations_research::SimpleMinCostFlow& min_cost_flow,
                                          const bam_api::AOSPairedReads& sequence, uint32_t M);
    static std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence,
                                                   uint32_t M);
    static std::unique_ptr<Solution> obtain_sequence(
        const bam_api::AOSPairedReads& sequence,
        const operations_research::SimpleMinCostFlow& min_cost_flow);

    bam_api::AOSPairedReads input_sequence_;
};
}  // namespace qmcp

#endif
