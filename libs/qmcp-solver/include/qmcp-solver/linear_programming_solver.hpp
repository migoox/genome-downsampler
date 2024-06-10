#ifndef QMCP_LINEAR_PROGRAMMING_SOLVER_HPP
#define QMCP_LINEAR_PROGRAMMING_SOLVER_HPP

#include <memory>

#include "solver.hpp"

namespace qmcp {
class LinearProgrammingSolver : public Solver {
   public:
    std::unique_ptr<Solution> solve(uint32_t max_coverage, bam_api::BamApi& bam_api) override;
    bool uses_quality_of_reads() override { return true; }
};

}  // namespace qmcp
#endif
