#ifndef QMCP_CPU_MIXED_INTEGER_SOLVER
#define QMCP_CPU_MIXED_INTEGER_SOLVER

#include <cstdint>
#include <memory>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/soa_paired_reads.hpp"
#include "qmcp-solver/qmcp_cpu_cost_scaling_solver.hpp"
#include "solver.hpp"

namespace qmcp
{
class CpuMixedIntegerSolver : public Solver
{
   public:
    CpuMixedIntegerSolver();
    std::unique_ptr<Solution> solve(uint32_t max_coverage, bam_api::BamApi& bam_api) override;
    bool uses_quality_of_reads() override { return true; }

   private:
    std::unique_ptr<Solver> helper_solver;
    uint32_t get_solution_reads_count(uint32_t max_coverage, bam_api::BamApi& bam_api);
    static std::unique_ptr<Solution> obtain_sequence(const bam_api::AOSPairedReads& sequence,
                                                     std::vector<int64_t>& solution);

    bam_api::AOSPairedReads input_sequence_;
};
}  // namespace qmcp

#endif
