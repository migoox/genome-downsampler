#include "app.hpp"

#include <CLI/App.hpp>
#include <CLI/Validators.hpp>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_api_config.hpp"
#include "bam-api/bam_api_config_builder.hpp"
#include "logging/log.hpp"
#include "qmcp-solver/solver.hpp"

App::App() {
    algorithms_names_ = get_algorithms_names();
    solver_ = solvers_map_["quasi-mcp-cpu"];
    add_main_command_options();

#ifdef TESTING_ENABLED
    solver_testers_names_ = get_solver_testers_names();
    test_subcmd_ = app_.add_subcommand("test", "Run unit tests.");
    add_test_subcommand_options();
    app_.require_subcommand(0, 1);
#endif
}

void App::App::add_main_command_options() {
    app_.add_option("INPUT_FILEPATH", input_file_path_, ".bam input file path. Required option.")
        ->check(CLI::ExistingFile);

    app_.add_option("MAX_COVERAGE", max_ref_coverage_,
                    "Maximum coverage per reference genome's base pair index.")
        ->check(CLI::PositiveNumber);

    app_.add_option("-o,--output", output_file_path_,
                    ".bam output file path. Default is \"output.bam\" in "
                    "input's directory.");

    // Logic to make positional arguments required when subcommand is not used
    app_.callback([&]() {
        if (app_.get_subcommand() == nullptr) {
            // If subcommand is not invoked, ensure both input and max-coverage are provided
            if (max_ref_coverage_ == 0) {
                throw CLI::ParseError("MAX_COVERAGE must be specified and integer bigger than 0",
                                      1);
            }

            if (input_file_path_.empty()) {
                throw CLI::ParseError("INPUT_FILEPATH must be specified", 1);
            }

            if (output_file_path_.empty()) {
                output_file_path_ = input_file_path_;
                output_file_path_.replace_filename("output.bam");
            }
        }
    });

    app_.add_option_function<std::string>(
            "-a,--algorithm",
            [this](const std::string& algorithm_name) {
                if (solvers_map_.find(algorithm_name) != solvers_map_.end()) {
                    solver_ = solvers_map_[algorithm_name];
                } else {
                    LOG_WITH_LEVEL(logging::ERROR)
                        << "Algorithm not found: " << algorithm_name << std::endl;
                }
            },
            "Algorithm to use. Default is \"quasi-mcp-cpu\"")
        ->check(CLI::IsMember(algorithms_names_));

    app_.add_option("-b,--bed", bed_path_,
                    ".bed amplicon bounds specification. It would be used to filter out or lower "
                    "the priority "
                    "of pair of sequences from different amplicons. The behaviour depends on "
                    "algorithm used.")
        ->check(CLI::ExistingFile);

    app_.add_option("-t,--tsv", tsv_path_,
                    ".tsv file which describes which of the (must be specified with this option) "
                    ".bed amplicon bounds should be paired together creating amplicon. used in "
                    "filtering or prioritizing pairs of sequences.")
        ->check(CLI::ExistingFile);

    app_.add_option("-p,--preprocessing-out", filtered_out_path_,
                    ".bam output file for reads, which were filtered out during preprocessing. It "
                    "can be useful for debugging.");

    app_.add_option("-l,--min-length", min_seq_length_,
                    "Minimal sequence length. Default is 90. Sequences "
                    "shorter than this "
                    "integer would be filtered before algorithm execution.")
        ->check(CLI::NonNegativeNumber);

    app_.add_option("-q,--min-mapq", min_mapq_,
                    "Minimal MAPQ value of the sequence. Default is 30. "
                    "Sequences with smaller MAPQ than this integer would be "
                    "filtered out before algorithm execution.")
        ->check(CLI::NonNegativeNumber);

    app_.add_option("-@,--threads", hts_thread_count_, "Set thread count for htslib read/write.")
        ->check(CLI::PositiveNumber);

    app_.add_flag("-v,--verbose", verbose_mode_,
                  "If specified app executes with additional logging.");
}

#ifdef TESTING_ENABLED
void App::App::add_test_subcommand_options() {
    test_subcmd_->add_option("-a,--algorithms", algorithms_to_test_, "Algorithms to test.")
        ->transform(CLI::IsMember(algorithms_names_));

    test_subcmd_->add_option("-t,--tests", solver_testers_, "Tests to run.")
        ->transform(CLI::IsMember(solver_testers_names_));

    test_subcmd_
        ->add_option(
            "-o,--outputs-dir", test_outputs_dir_,
            "Directory for tests outputs. Each SolverTester should have some optional output data "
            "to save for debugging purposes. This option enables it and specifies its directory")
        ->check(CLI::ExistingDirectory);
}

void App::run_tests() {
    std::vector<std::string> solvers_to_test =
        algorithms_to_test_.empty() ? algorithms_names_ : algorithms_to_test_;
    std::vector<std::string> tests_to_run =
        solver_testers_.empty() ? solver_testers_names_ : solver_testers_;
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
            solver_testers_map_[test]->test(*solvers_map_[solver],
                                            alg_tester_outputs_directory_path);

            LOG_WITH_LEVEL(logging::INFO) << "\t\t PASSED";
        }
    }
}

std::vector<std::string> App::get_solver_testers_names() const {
    std::vector<std::string> all_tests;
    all_tests.reserve(solver_testers_map_.size());

    for (const auto& mapping : solver_testers_map_) {
        all_tests.push_back(mapping.first);
    }
    return all_tests;
}

#endif

void App::parse(int argc, char** argv) {
    app_.parse(argc, argv);

    SET_LOG_LEVEL(verbose_mode_ ? logging::DEBUG : logging::INFO);
}

int App::exit(const CLI::ParseError& e) { return app_.exit(e); }

void App::execute() {
#ifdef TESTING_ENABLED
    if (*test_subcmd_) {
        run_tests();
        return;
    }
#endif

    bam_api::BamApiConfigBuilder config_buider;

    config_buider.add_hts_thread_count(hts_thread_count_);
    config_buider.add_min_mapq(min_mapq_);
    config_buider.add_min_seq_length(min_seq_length_);

    if (!bed_path_.empty()) {
        if (solver_->uses_quality_of_reads()) {
            config_buider.add_amplicon_filtering(bam_api::AmpliconBehaviour::GRADE, bed_path_,
                                                 tsv_path_);
        } else {
            config_buider.add_amplicon_filtering(bam_api::AmpliconBehaviour::FILTER, bed_path_,
                                                 tsv_path_);
        }
    }

    bam_api::BamApi bam_api(input_file_path_, config_buider.build());

    auto start = std::chrono::high_resolution_clock::now();
    std::unique_ptr<qmcp::Solution> solution = solver_->solve(max_ref_coverage_, bam_api);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    LOG_WITH_LEVEL(logging::DEBUG) << "solve took " << elapsed.count() << " seconds";

    std::vector<bam_api::BAMReadId> paired_solution = bam_api.find_pairs(*solution);
    bam_api.write_paired_reads(output_file_path_, paired_solution);

    if (!filtered_out_path_.empty()) {
        bam_api.write_bam_api_filtered_out_reads(filtered_out_path_);
    }
}

std::vector<std::string> App::get_algorithms_names() const {
    std::vector<std::string> all_algorithms;
    all_algorithms.reserve(solvers_map_.size());

    for (const auto& mapping : solvers_map_) {
        all_algorithms.push_back(mapping.first);
    }
    return all_algorithms;
}
