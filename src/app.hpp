#pragma once

#include <CLI/CLI.hpp>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>

#include "qmcp-solver/cuda_max_flow_solver.hpp"
#include "qmcp-solver/sequential_cost_scaling_network_solver.hpp"
#include "qmcp-solver/sequential_max_flow_solver.hpp"
#include "qmcp-solver/solver.hpp"
#include "qmcp-solver/test_solver.hpp"

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
    bool verbose_mode_ = false;
    uint32_t min_seq_length_ = kDefaultMinSeqLength;
    uint32_t min_seq_mapq_ = kDefaultMinSeqQmap;
    uint32_t max_ref_coverage_;
    std::unique_ptr<qmcp::Solver> solver_;
    std::filesystem::path input_file_path_;
    std::filesystem::path output_file_path_;
    std::filesystem::path csv_historical_runs_file_path_;
    std::filesystem::path bed_amplicon_file_path_;
    std::filesystem::path tsv_amplicon_pairings_file_path_;
    std::map<std::string, std::unique_ptr<qmcp::Solver>> solvers_map_;

    void FillSolversMap() {
        // To add an algorithm emplace its unique_ptr to the map with
        // identification name for CLI
        solvers_map_.emplace("test", std::make_unique<qmcp::TestSolver>());
        solvers_map_.emplace("sequential-cost-scaling",
                             std::make_unique<qmcp::SequentialCostScalingNetworkSolver>());
        solvers_map_.emplace("cuda-max-flow", std::make_unique<qmcp::CudaMaxFlowSolver>());
        solvers_map_.emplace("sequential-max-flow",
                             std::make_unique<qmcp::SequentialMaxFlowSolver>());
    }
};
