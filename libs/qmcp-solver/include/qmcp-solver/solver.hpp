#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cstdint>
#include <memory>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
namespace qmcp {

typedef std::vector<bam_api::BAMReadId> Solution;

class Solver {
   public:
    virtual ~Solver() = default;
    virtual std::unique_ptr<Solution> solve(uint32_t max_coverage, bam_api::BamApi& bam_api) = 0;
};
}  // namespace qmcp

#endif
