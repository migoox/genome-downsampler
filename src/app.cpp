#include "app.hpp"

#include <CLI/App.hpp>
#include <CLI/Validators.hpp>
#include <chrono>
#include <filesystem>
#include <memory>
#include <string>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_api_config.hpp"
#include "bam-api/bam_api_config_builder.hpp"
#include "logging/log.hpp"
#include "qmcp-solver/solver.hpp"
#include "test_command.hpp"

App::App() : solver_name_(kDefaultSolver) {
    add_main_command_options();

#ifdef TESTING_ENABLED
    test_command_ = std::make_unique<TestCommand>(app_, solver_manager_);
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
        // Check if any subcommands were parsed
        if (app_.get_subcommands().empty()) {
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

            execute();
        }
    });

    app_.add_option_function<std::string>(
            "-a,--algorithm",
            [this](const std::string& algorithm_name) {
                if (solver_manager_.contains(algorithm_name)) {
                    solver_name_ = algorithm_name;
                } else {
                    LOG_WITH_LEVEL(logging::ERROR)
                        << "Algorithm not found: " << algorithm_name << std::endl;
                }
            },
            "Algorithm to use. Default is \"quasi-mcp-cpu\"")
        ->check(CLI::IsMember(solver_manager_.get_names()));

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

void App::parse(int argc, char** argv) { app_.parse(argc, argv); }

int App::exit(const CLI::ParseError& e) { return app_.exit(e); }

void App::execute() {
    bam_api::BamApiConfigBuilder config_buider;

    config_buider.add_hts_thread_count(hts_thread_count_);
    config_buider.add_min_mapq(min_mapq_);
    config_buider.add_min_seq_length(min_seq_length_);

    if (!bed_path_.empty()) {
        if (solver_manager_.get(solver_name_).uses_quality_of_reads()) {
            config_buider.add_amplicon_filtering(bam_api::AmpliconBehaviour::GRADE, bed_path_,
                                                 tsv_path_);
        } else {
            config_buider.add_amplicon_filtering(bam_api::AmpliconBehaviour::FILTER, bed_path_,
                                                 tsv_path_);
        }
    }

    bam_api::BamApi bam_api(input_file_path_, config_buider.build());

    auto start = std::chrono::high_resolution_clock::now();

    std::unique_ptr<qmcp::Solution> solution =
        solver_manager_.get(solver_name_).solve(max_ref_coverage_, bam_api);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    LOG_WITH_LEVEL(logging::DEBUG) << "solve took " << elapsed.count() << " seconds";

    std::vector<bam_api::BAMReadId> paired_solution = bam_api.find_pairs(*solution);
    bam_api.write_paired_reads(output_file_path_, paired_solution);

    if (!filtered_out_path_.empty()) {
        bam_api.write_bam_api_filtered_out_reads(filtered_out_path_);
    }
}
