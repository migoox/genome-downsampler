#pragma once

#include <filesystem>
#include <vector>

#include "bam-api/bam_paired_reads.hpp"
#include "solver.hpp"

namespace qmcp {
class SequentialCostScalingNetworkSolver : public Solver {
   public:
    void import_reads(const std::filesystem::path& input, uint32_t min_seq_length,
                uint32_t min_seq_mapq) override;
    void solve(uint32_t max_coverage) override;
    void export_reads(const std::filesystem::path& output) override;

   private:
    bam_api::AOSPairedReads input_sequence_;
    std::filesystem::path input_;
    std::vector<bam_api::ReadIndex> output_sequence_;
};
}  // namespace qmcp
