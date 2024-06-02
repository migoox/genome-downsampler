#include "app.hpp"

#include <CLI/App.hpp>
#include <CLI/Validators.hpp>
#include <filesystem>
#include <memory>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_api_config.hpp"
#include "bam-api/bam_api_config_builder.hpp"
#include "logging/log.hpp"
#include "qmcp-solver/solver.hpp"

App::App() {
    app_.add_option("-i,--input", input_file_path_, ".bam input file path. Required option.")
        ->required()
        ->check(CLI::ExistingFile);

    std::vector<std::string> possible_algorithms;
    possible_algorithms.reserve(solvers_map_.size());
    for (const auto& mapping : solvers_map_) {
        possible_algorithms.push_back(mapping.first);
    }
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
        ->required()
        ->check(CLI::IsMember(possible_algorithms));

    app_.add_option("-M,--max-coverage", max_ref_coverage_,
                    "Maximum coverage per reference genome's base pair index.")
        ->required()
        ->check(CLI::NonNegativeNumber);

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
                    ".bam output file for reads, which was filtered out during preprocessing. It "
                    "can be useful for debugging.");

    app_.add_option("-l,--min-length", min_seq_length_,
                    "Minimal sequence length. Default is 90. Sequences "
                    "shorter than this "
                    "integer would be filtered before algorithm execution.")
        ->check(CLI::NonNegativeNumber);

    app_.add_option("-q,--min-mapq", min_mapq_,
                    "Minimal MAPQ value of the sequence. Default is 30. "
                    "Sequences with smaller MAPQ than this integer would be "
                    "filtered before algorithm execution.")
        ->check(CLI::NonNegativeNumber);

    app_.add_flag("-v,--verbose", verbose_mode_,
                  "If specified app_ executes with additional logging.");
}

void App::Parse(int argc, char** argv) {
    app_.parse(argc, argv);

    if (output_file_path_.empty()) {
        output_file_path_ = input_file_path_;
        output_file_path_.replace_filename("output.bam");
    }

    SET_LOG_LEVEL(verbose_mode_ ? logging::DEBUG : logging::INFO);
}

int App::Exit(const CLI::ParseError& e) { return app_.exit(e); }

// TODO(mytkom): Add DEBUG chrono logging and other informational logs
void App::Solve() {
    bam_api::BamApiConfigBuilder config_buider;

    config_buider.add_min_mapq(min_mapq_);
    config_buider.add_min_seq_length(min_seq_length_);

    if (!bed_path_.empty()) {
        if (solver_->uses_quality_of_reads()) {
            config_buider.add_amplicon_filtering(bam_api::AmpliconBehaviour::GRADE, bed_path_, tsv_path_);
        } else {
            config_buider.add_amplicon_filtering(bam_api::AmpliconBehaviour::FILTER, bed_path_, tsv_path_);
        }
    }

    bam_api::BamApi bam_api(input_file_path_, config_buider.get_config());
    std::unique_ptr<qmcp::Solution> solution = solver_->solve(max_ref_coverage_, bam_api);
    std::vector<bam_api::BAMReadId> paired_solution = bam_api.find_pairs(*solution);
    bam_api.write_paired_reads(output_file_path_, paired_solution);

    if (!filtered_out_path_.empty()) {
        bam_api.write_bam_api_filtered_out_reads(filtered_out_path_);
    }
}
