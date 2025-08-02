#ifndef TEST_COMMAND_HPP
#define TEST_COMMAND_HPP

#include <CLI/App.hpp>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "qmcp-solver/solver.hpp"
#include "tests/coverage_tester.hpp"
#include "tests/solver_tester.hpp"

class TestCommand {
   public:
    TestCommand(CLI::App& app,
                const std::map<std::string, std::shared_ptr<qmcp::Solver>>& solvers_map);
    void run(const std::map<std::string, std::shared_ptr<qmcp::Solver>>& solvers_map);

   private:
    std::map<std::string, std::unique_ptr<test::SolverTester>> solver_testers_map_;
    std::vector<std::string> solver_testers_names_;
    std::vector<std::string> algorithms_names_;
    std::vector<std::string> algorithms_to_test_;
    std::vector<std::string> solver_testers_;
    std::filesystem::path test_outputs_dir_;
    CLI::App* test_subcmd_;
};

#endif
