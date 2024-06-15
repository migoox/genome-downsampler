#include "app.hpp"

#include <CLI/App.hpp>
#include <CLI/Validators.hpp>
#include <chrono>
#include <filesystem>
#include <memory>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_api_config.hpp"
#include "bam-api/bam_api_config_builder.hpp"
#include "logging/log.hpp"
#include "qmcp-solver/solver.hpp"

App::App() {
    std::vector<std::string> all_algorithms = GetAllAlgorithms();
    std::vector<std::string> all_solver_testers = GetAllSolverTesters();
    app_.fallthrough();

    app_.add_option("INPUT_FILE", input_file_path_, ".bam input file path. Required option.")
        ->check(CLI::ExistingFile);

    app_.add_option("MAX_COVERAGE", max_ref_coverage_,
                    "Maximum coverage per reference genome's base pair index.")
        ->check(CLI::PositiveNumber);

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
            "Algorithm to use.")
        ->check(CLI::IsMember(all_algorithms));

    app_.add_option("-o,--output", output_file_path_,
                    ".bam output file path. Default is \"output.bam\" in "
                    "input's directory.");

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
                  "If specified app_ executes with additional logging.");

    // Add test subcommand
    test_subcmd_ = app_.add_subcommand("test", "Run unit tests.");

    test_subcmd_->add_option("-a,--algorithms", algorithms_, "Algorithms to test.")
        ->transform(CLI::IsMember(all_algorithms));

    test_subcmd_->add_option("-t,--tests", solver_testers_, "Tests to run.")
        ->transform(CLI::IsMember(all_solver_testers));

    // Logic to make positional arguments required when subcommand is not used
    app_.callback([&]() {
        if (app_.get_subcommand() == nullptr) {
            // If subcommand is not invoked, ensure both input and max-coverage are provided
            if (max_ref_coverage_ == 0) {
                throw CLI::ParseError("MAX_COVERAGE must be specified and integer bigger than 0", 1);
            }

            if (input_file_path_.empty()) {
                throw CLI::ParseError("INPUT_FILE must be specified", 1);
            }
        }
    });

    app_.require_subcommand(0, 1);
}

void App::Parse(int argc, char** argv) {
    app_.parse(argc, argv);

    if (!test_subcmd_ && output_file_path_.empty()) {
        output_file_path_ = input_file_path_;
        output_file_path_.replace_filename("output.bam");
    }

    SET_LOG_LEVEL(verbose_mode_ ? logging::DEBUG : logging::INFO);
}

int App::Exit(const CLI::ParseError& e) { return app_.exit(e); }

void App::Solve() {
    if (*test_subcmd_) {
        RunTests();
        return;
    }

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

void App::RunTests() {
    std::vector<std::string> solvers_to_test =
        algorithms_.empty() ? GetAllAlgorithms() : algorithms_;
    std::vector<std::string> tests_to_run =
        solver_testers_.empty() ? GetAllSolverTesters() : solver_testers_;

    for (const auto& test : tests_to_run) {
        LOG_WITH_LEVEL(logging::INFO) << "Running test " << test;
        for (const auto& solver : solvers_to_test) {
            LOG_WITH_LEVEL(logging::INFO) << "\ton algorithm " << solver;

            // Run tester
            solver_testers_map_[test]->test(*solvers_map_[solver]);
        }
    }
}

std::vector<std::string> App::GetAllAlgorithms() const {
    std::vector<std::string> all_algorithms;
    all_algorithms.reserve(solvers_map_.size());

    for (const auto& mapping : solvers_map_) {
        all_algorithms.push_back(mapping.first);
    }
    return all_algorithms;
}

std::vector<std::string> App::GetAllSolverTesters() const {
    std::vector<std::string> all_tests;
    all_tests.reserve(solver_testers_map_.size());

    for (const auto& mapping : solver_testers_map_) {
        all_tests.push_back(mapping.first);
    }
    return all_tests;
}
