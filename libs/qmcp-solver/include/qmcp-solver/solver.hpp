#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cstdint>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
namespace qmcp {
class Solver {
   public:
    explicit Solver(bam_api::BamApi& bam_api) : bam_api_(bam_api) {}
    virtual ~Solver() = default;
    virtual std::vector<bam_api::BAMReadId> solve(uint32_t max_coverage) = 0;

   protected:
    bam_api::BamApi& bam_api_;
};
}  // namespace qmcp

#endif
