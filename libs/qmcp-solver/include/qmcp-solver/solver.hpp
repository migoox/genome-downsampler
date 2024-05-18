#pragma once

#include <filesystem>

#include "logger.hpp"

namespace qmcp {
class Solver {
   public:
    virtual ~Solver() = default;
    virtual void Import(const std::filesystem::path& input,
                        uint32_t min_seq_length, uint32_t min_seq_mapq) = 0;
    virtual void Solve(uint32_t max_coverage) = 0;
    virtual void Export(const std::filesystem::path& output) = 0;
};
}  // namespace qmcp
