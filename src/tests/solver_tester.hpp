#ifndef TESTS_SOLVER_TESTER
#define TESTS_SOLVER_TESTER

#include <filesystem>
#include <memory>

#include "qmcp-solver/solver.hpp"

namespace test {

class SolverTester {
   public:
    virtual ~SolverTester() = default;
    virtual void test(const std::unique_ptr<qmcp::Solver>& solver,
                      std::filesystem::path& outputs_dir_path_) = 0;
};

}  // namespace test

#endif
