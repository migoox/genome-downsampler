#pragma once

#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/solver.hpp"

namespace qmcp {
class TestSolver : public Solver {
   public:
    void import_reads(const std::filesystem::path& input, uint32_t min_seq_length,
                      uint32_t min_seq_mapq, const std::filesystem::path& bed_amplicon, const std::filesystem::path& tsv_amplicon) override;
    void solve(uint32_t max_coverage) override;
    void export_reads(const std::filesystem::path& output) override;

   private:
    bam_api::SOAPairedReads paired_reads_;
    std::vector<bam_api::BAMReadId> solution_;
    std::filesystem::path input_;
};
}  // namespace qmcp
