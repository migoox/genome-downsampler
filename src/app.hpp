#ifndef APP_HPP
#define APP_HPP

#include <CLI/CLI.hpp>
#include <cstdint>
#include <filesystem>
#include <string>

#include "qmcp-solver/quasi_mcp_cuda_max_flow_solver.hpp"
#include "solver_manager.hpp"
#include "test_command.hpp"

class App {
    static constexpr uint32_t kDefaultMinSeqLength = 90;
    static constexpr uint32_t kDefaultAmpOverflow = 0;
    static constexpr uint32_t kDefaultMinSeqMAPQ = 30;
    static constexpr float kDefaultMinAlignment = 0.5;
    static constexpr uint32_t kDefaultThreadCount = 2;
    static constexpr const char* kDefaultSolver = "quasi-mcp-cpu";

   public:
    App();
    void parse(int argc, char** argv);
    void execute();
    int exit(const CLI::ParseError& e);

   private:
    CLI::App app_;
    SolverManager solver_manager_;
    uint32_t hts_thread_count_ = kDefaultThreadCount;
    uint32_t min_mapq_ = kDefaultMinSeqMAPQ;
    uint32_t amp_overflow_ = kDefaultAmpOverflow;
    float min_alignment_ = kDefaultMinAlignment;
    uint32_t min_seq_length_ = kDefaultMinSeqLength;
    bool verbose_mode_ = false;
    uint32_t max_ref_coverage_ = 0;
    std::string solver_name_;
    std::filesystem::path input_file_path_;
    std::filesystem::path output_file_path_;
    std::filesystem::path filtered_out_path_;
    std::filesystem::path bed_path_;
    std::filesystem::path tsv_path_;

#ifdef TESTING_ENABLED
    std::unique_ptr<TestCommand> test_command_;
#endif

    void add_main_command_options();
};

#endif
