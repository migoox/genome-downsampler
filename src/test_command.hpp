#ifndef TEST_COMMAND_HPP
#define TEST_COMMAND_HPP

#include <CLI/App.hpp>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "qmcp-solver/solver.hpp"
#include "solver_manager.hpp"
#include "tester_manager.hpp"

class TestCommand {
   public:
    TestCommand(CLI::App& app, const SolverManager& solver_manager);
    void run(const SolverManager& solver_manager);

   private:
    TesterManager tester_manager_;
    std::vector<std::string> algorithms_to_test_;
    std::vector<std::string> solver_testers_;
    std::filesystem::path test_outputs_dir_;
    CLI::App* test_subcmd_;
};

#endif
