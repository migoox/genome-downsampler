#include "test_command.hpp"

#include <CLI/App.hpp>
#include <functional>
#include <string>

#include "logging/log.hpp"

TestCommand::TestCommand(CLI::App& app, const SolverManager& solver_manager) {
    // Set test subcmd options
    test_subcmd_ = app.add_subcommand("test", "Run unit tests.");
    test_subcmd_->callback([&]() { run(solver_manager); });

    test_subcmd_->add_option("-a,--algorithms", algorithms_to_test_, "Algorithms to test.")
        ->transform(CLI::IsMember(solver_manager.get_names()));

    test_subcmd_->add_option("-t,--tests", solver_testers_, "Tests to run.")
        ->transform(CLI::IsMember(tester_manager_.get_names()));

    test_subcmd_
        ->add_option(
            "-o,--outputs-dir", test_outputs_dir_,
            "Directory for tests outputs. Each SolverTester should have some optional output data "
            "to save for debugging purposes. This option enables it and specifies its directory")
        ->check(CLI::ExistingDirectory);

    // Configure main app
    app.require_subcommand(0, 1);
}

void TestCommand::run(const SolverManager& solver_manager) {
    const std::vector<std::string>& solvers_to_test =
        algorithms_to_test_.empty() ? solver_manager.get_names() : algorithms_to_test_;
    const std::vector<std::string>& tests_to_run =
        solver_testers_.empty() ? tester_manager_.get_names() : solver_testers_;
    bool with_output = !test_outputs_dir_.empty();
    std::filesystem::path tester_outputs_directory_path;
    std::filesystem::path alg_tester_outputs_directory_path;

    if (with_output && !std::filesystem::exists(test_outputs_dir_)) {
        LOG_WITH_LEVEL(logging::ERROR) << "Directory: " << test_outputs_dir_ << " does not exist!";
        std::exit(EXIT_FAILURE);
    }

    for (const std::string& test : tests_to_run) {
        if (with_output) {
            tester_outputs_directory_path = test_outputs_dir_ / test;
            if (!std::filesystem::exists(tester_outputs_directory_path)) {
                std::filesystem::create_directory(tester_outputs_directory_path);
            }
        }
        LOG_WITH_LEVEL(logging::INFO) << "Running test " << test;
        for (const std::string& solver : solvers_to_test) {
            LOG_WITH_LEVEL(logging::INFO) << "\ton algorithm " << solver;

            if (with_output) {
                alg_tester_outputs_directory_path = tester_outputs_directory_path / solver;
                if (!std::filesystem::exists(alg_tester_outputs_directory_path)) {
                    std::filesystem::create_directory(alg_tester_outputs_directory_path);
                }
            }

            // Run tester
            tester_manager_.get(test).test(solver_manager.get(solver),
                                            alg_tester_outputs_directory_path);

            LOG_WITH_LEVEL(logging::INFO) << "\t\t PASSED";
        }
    }
}
