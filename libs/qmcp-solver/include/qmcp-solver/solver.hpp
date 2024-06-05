#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cstdint>
#include <memory>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/read.hpp"

namespace qmcp {

typedef std::vector<bam_api::ReadIndex> Solution;

class Solver {
   public:
    virtual ~Solver() = default;
    virtual std::unique_ptr<Solution> solve(uint32_t max_coverage, bam_api::BamApi& bam_api) = 0;
    virtual bool uses_quality_of_reads() = 0;
};
}  // namespace qmcp

#endif
