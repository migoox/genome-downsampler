#pragma once

#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/solver.hpp"

namespace qmcp {
class TestSolver : public Solver {
   public:
    static constexpr uint32_t kReadsFromEndOffset = 1000;
    static constexpr uint32_t kReadsCount = 998;

    void Import(const std::filesystem::path& input, uint32_t min_seq_length,
                uint32_t min_seq_mapq) override;
    void Solve(uint32_t max_coverage) override;
    void Export(const std::filesystem::path& output) override;

   private:
    bam_api::SOAPairedReads paired_reads_;
    std::vector<bam_api::ReadIndex> solution_;
    std::filesystem::path input_;
};
}  // namespace qmcp
