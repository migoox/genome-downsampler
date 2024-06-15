#ifndef TESTS_SOLVER_TESTER
#define TESTS_SOLVER_TESTER

#include <filesystem>
#include "qmcp-solver/solver.hpp"

namespace test {

class SolverTester {
   public:
    virtual ~SolverTester() = default;
    virtual void test(qmcp::Solver& solver, std::filesystem::path& outputs_dir_path_) = 0;
};

}  // namespace test

#endif
