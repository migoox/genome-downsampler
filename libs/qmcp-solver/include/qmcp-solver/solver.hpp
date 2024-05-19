#pragma once

#include <filesystem>

namespace qmcp {
class Solver {
   public:
    virtual ~Solver() = default;
    virtual void import_reads(const std::filesystem::path& input,
                        uint32_t min_seq_length, uint32_t min_seq_mapq) = 0;
    virtual void solve(uint32_t max_coverage) = 0;
    virtual void export_reads(const std::filesystem::path& output) = 0;
};
}  // namespace qmcp
