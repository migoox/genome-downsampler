#ifndef APP_HPP
#define APP_HPP

#include <CLI/CLI.hpp>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>

#include "bam-api/bam_api.hpp"
#include "qmcp-solver/cuda_max_flow_solver.hpp"
#include "qmcp-solver/sequential_cost_scaling_network_solver.hpp"
#include "qmcp-solver/sequential_max_flow_solver.hpp"
#include "qmcp-solver/solver.hpp"
#include "qmcp-solver/test_solver.hpp"

class App {
    static constexpr uint32_t kDefaultMinSeqLength = 90;
    static constexpr uint32_t kDefaultMinSeqMAPQ = 30;

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
    std::shared_ptr<qmcp::Solver> solver_;
    std::filesystem::path input_file_path_;
    std::filesystem::path output_file_path_;
    std::filesystem::path csv_historical_runs_file_path_;
    std::map<std::string, std::shared_ptr<qmcp::Solver>> solvers_map_{
        {"test", std::make_shared<qmcp::TestSolver>()},
        {"sequential-cost-scaling", std::make_shared<qmcp::SequentialCostScalingNetworkSolver>()},
        {"cuda-max-flow", std::make_shared<qmcp::CudaMaxFlowSolver>()},
        {"sequential-max-flow", std::make_shared<qmcp::SequentialMaxFlowSolver>()},
    };
};

#endif
