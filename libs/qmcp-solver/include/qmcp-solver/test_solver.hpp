#include <memory>
#include "bam-api/bam_api.hpp"
#ifndef QMCP_TEST_SOLVER_HPP
#define QMCP_TEST_SOLVER_HPP()

#include "qmcp-solver/solver.hpp"

namespace qmcp {
class TestSolver : public Solver {
   public:
    std::unique_ptr<Solution> solve(uint32_t max_coverage, bam_api::BamApi& bam_api) override;
    bool uses_quality_of_reads() override { return false; }
};
}  // namespace qmcp

#endif
