#ifndef QMCP_CPU_CRAWLING_SOLVER
#define QMCP_CPU_CRAWLING_SOLVER

#include <cstdint>
#include <memory>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/soa_paired_reads.hpp"
#include "solver.hpp"

namespace qmcp
{
class CpuCrawlingSolver : public Solver
{
   public:
    std::unique_ptr<Solution> solve(uint32_t max_coverage, bam_api::BamApi& bam_api) override;
    bool uses_quality_of_reads() override { return true; }

   private:
    static std::unique_ptr<Solution> obtain_sequence(const bam_api::AOSPairedReads& sequence,
                                                     std::vector<int64_t>& solution);

    bam_api::AOSPairedReads input_sequence_;
};
}  // namespace qmcp

#endif
