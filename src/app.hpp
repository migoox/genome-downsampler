#ifndef APP_HPP
#define APP_HPP

#include <CLI/CLI.hpp>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>

#include "bam-api/bam_api.hpp"
#include "qmcp-solver/solver.hpp"
#include "solver_factory_functions.hpp"

class App {
    static constexpr uint32_t kDefaultMinSeqLength = 90;
    static constexpr uint32_t kDefaultMinSeqQmap = 30;

   public:
    App();
    void Parse(int argc, char** argv);
    int Exit(const CLI::ParseError& e);
    void Solve();

   private:
    CLI::App app_;
    bam_api::BamApiConfig bam_api_config_;
    bool verbose_mode_ = false;
    uint32_t max_ref_coverage_;
    std::function<std::unique_ptr<qmcp::Solver>(bam_api::BamApi&)> solver_factory_function_;
    std::filesystem::path input_file_path_;
    std::filesystem::path output_file_path_;
    std::filesystem::path csv_historical_runs_file_path_;
    std::map<std::string, std::function<std::unique_ptr<qmcp::Solver>(bam_api::BamApi&)>>
        solvers_map_{
            {"test", solver_factory_functions::createTestSolver},
            {"sequential-cost-scaling",
             solver_factory_functions::createSequentialCostScalingNetworkSolver},
            {"cuda-max-flow", solver_factory_functions::createCudaMaxFlowSolver},
            {"sequential-max-flow", solver_factory_functions::createSequentialMaxFlowSolver},
        };
};

#endif
