#ifndef APP_HPP
#define APP_HPP

#include <CLI/CLI.hpp>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>

#include "tests/solver_tester.hpp"
#include "tests/coverage_tester.hpp"
#include "qmcp-solver/cuda_max_flow_solver.hpp"
#include "qmcp-solver/sequential_cost_scaling_network_solver.hpp"
#include "qmcp-solver/sequential_max_flow_solver.hpp"
#include "qmcp-solver/solver.hpp"
#include "qmcp-solver/test_solver.hpp"

class App {
    static constexpr uint32_t kDefaultMinSeqLength = 90;
    static constexpr uint32_t kDefaultMinSeqMAPQ = 30;
    static constexpr uint32_t kDefaultThreadCount = 2;

   public:
    App();
    void Parse(int argc, char** argv);
    int Exit(const CLI::ParseError& e);
    void Solve();

   private:
    CLI::App app_;
    uint32_t hts_thread_count_ = kDefaultThreadCount;
    uint32_t min_mapq_ = kDefaultMinSeqMAPQ;
    uint32_t min_seq_length_ = kDefaultMinSeqLength;
    bool verbose_mode_ = false;
    uint32_t max_ref_coverage_;
    std::shared_ptr<qmcp::Solver> solver_;
    std::filesystem::path input_file_path_;
    std::filesystem::path output_file_path_;
    std::filesystem::path filtered_out_path_;
    std::filesystem::path bed_path_;
    std::filesystem::path tsv_path_;
    std::map<std::string, std::shared_ptr<qmcp::Solver>> solvers_map_{
        {"test", std::make_shared<qmcp::TestSolver>()},
        {"sequential-cost-scaling", std::make_shared<qmcp::SequentialCostScalingNetworkSolver>()},
        {"cuda-max-flow", std::make_shared<qmcp::CudaMaxFlowSolver>()},
        {"sequential-max-flow", std::make_shared<qmcp::SequentialMaxFlowSolver>()},
    };
    std::map<std::string, std::shared_ptr<test::SolverTester>> solver_testers_map_{
        {"coverage", std::make_shared<test::CoverageTester>()},
    };
};

#endif
