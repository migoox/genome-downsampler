#ifndef QMCP_QUASI_MCP_CPU_MAX_FLOW_SOLVER_HPP
#define QMCP_QUASI_MCP_CPU_MAX_FLOW_SOLVER_HPP
#include <ortools/graph/max_flow.h>

#include <memory>
#include <vector>

#include "bam-api/bam_api.hpp"

#include "solver.hpp"

namespace qmcp {

class QuasiMcpCpuMaxFlowSolver : public Solver {
   public:
    std::unique_ptr<Solution> solve(uint32_t max_coverage, bam_api::BamApi& bam_api) override;
    bool uses_quality_of_reads() override { return false; }

   private:
    static void create_network_flow_graph(operations_research::SimpleMaxFlow& max_flow,
                                          const bam_api::AOSPairedReads& sequence, unsigned int M);
    static std::vector<int> create_demand_function(const bam_api::AOSPairedReads& sequence,
                                                   unsigned int M);
    static std::unique_ptr<Solution> obtain_sequence(
        const bam_api::AOSPairedReads& sequence, const operations_research::SimpleMaxFlow& max_flow);
    static std::vector<int> create_b_function(const bam_api::AOSPairedReads& sequence,
                                              unsigned int M);

    bam_api::AOSPairedReads input_sequence_;
};

}  // namespace qmcp

#endif
