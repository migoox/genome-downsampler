#ifndef TESTS_SOLVER_TESTER
#define TESTS_SOLVER_TESTER

#include "qmcp-solver/solver.hpp"

namespace test {

class SolverTester {
   public:
    virtual ~SolverTester() = default;
    virtual void test(qmcp::Solver& solver) = 0;
};

}  // namespace test

#endif
