#ifndef QMCP_TEST_SOLVER_HPP
#define QMCP_TEST_SOLVER_HPP()

#include <vector>

#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/solver.hpp"

namespace qmcp {
class TestSolver : public Solver {
   public:
    using Solver::Solver;
    std::vector<bam_api::BAMReadId> solve(uint32_t max_coverage) override;
};
}  // namespace qmcp

#endif
