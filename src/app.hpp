#ifndef APP_HPP
#define APP_HPP

#include <CLI/CLI.hpp>
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>

#include "qmcp-solver/mcp_cpu_cost_scaling_solver.hpp"
#include "qmcp-solver/qmcp_cpu_cost_scaling_solver.hpp"
#include "qmcp-solver/quasi_mcp_cpu_max_flow_solver.hpp"
#include "qmcp-solver/quasi_mcp_cuda_max_flow_solver.hpp"
#include "qmcp-solver/solver.hpp"
#include "qmcp-solver/test_solver.hpp"

#ifdef TESTING_ENABLED
#include "tests/coverage_tester.hpp"
#include "tests/solver_tester.hpp"
#endif

class App {
    static constexpr uint32_t kDefaultMinSeqLength = 90;
    static constexpr uint32_t kDefaultMinSeqMAPQ = 30;
    static constexpr uint32_t kDefaultThreadCount = 2;

   public:
    App();
    void parse(int argc, char** argv);
    void execute();
    int exit(const CLI::ParseError& e);

   private:
    std::map<std::string, std::shared_ptr<qmcp::Solver>> solvers_map_{
        {"quasi-mcp-cpu", std::make_shared<qmcp::QuasiMcpCpuMaxFlowSolver>()},
        {"mcp-cpu", std::make_shared<qmcp::McpCpuCostScalingSolver>()},
        {"qmcp-cpu", std::make_shared<qmcp::QmcpCpuCostScalingSolver>()},
        {"quasi-mcp-cuda", std::make_shared<qmcp::QuasiMcpCudaMaxFlowSolver>()},
    };
    std::vector<std::string> algorithms_names_;
    CLI::App app_;
    CLI::App* test_subcmd_;
    uint32_t hts_thread_count_ = kDefaultThreadCount;
    uint32_t min_mapq_ = kDefaultMinSeqMAPQ;
    uint32_t min_seq_length_ = kDefaultMinSeqLength;
    bool verbose_mode_ = false;
    uint32_t max_ref_coverage_ = 0;
    std::shared_ptr<qmcp::Solver> solver_;
    std::filesystem::path input_file_path_;
    std::filesystem::path output_file_path_;
    std::filesystem::path filtered_out_path_;
    std::filesystem::path bed_path_;
    std::filesystem::path tsv_path_;

    void add_main_command_options();
    std::vector<std::string> get_algorithms_names() const;

#ifdef TESTING_ENABLED
    std::map<std::string, std::shared_ptr<test::SolverTester>> solver_testers_map_{
        {"coverage", std::make_shared<test::CoverageTester>()},
    };
    std::vector<std::string> solver_testers_names_;
    std::vector<std::string> algorithms_to_test_;
    std::vector<std::string> solver_testers_;
    std::filesystem::path test_outputs_dir_;

    void add_test_subcommand_options();
    void run_tests();
    std::vector<std::string> get_solver_testers_names() const;
#endif
};

#endif
