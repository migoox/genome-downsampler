#include "app.hpp"

App::App() {
    FillSolversMap();

    app_.add_option("-i,--input", input_file_path_,
                    ".bam input file path. Required option.")
        ->required()
        ->check(CLI::ExistingFile);

    std::vector<std::string> possible_algorithms;
    possible_algorithms.reserve(solvers_map_.size());
    for (const auto& mapping : solvers_map_) {
        possible_algorithms.push_back(mapping.first);
    }
    app_.add_option_function<std::string>(
            "-a,--algorithm",
            [this](const std::string &algorithm_name) {
                if (solvers_map_.find(algorithm_name) != solvers_map_.end()) {
                    solver_ = std::move(solvers_map_[algorithm_name]);
                } else {
                    std::cerr << "Algorithm not found: " << algorithm_name
                              << std::endl;
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

    app_.add_option("-c,--csv-history", csv_historical_runs_file_path_,
                    ".csv historical runs data file path. If it is not "
                    "specified, no historical data would be saved!");

    app_.add_option("-l,--min-length", min_seq_length_,
                    "Minimal sequence length. Default is 90. Sequences "
                    "shorter than this "
                    "integer would be filtered before algorithm execution.")
        ->check(CLI::NonNegativeNumber);

    app_.add_option("-q,--min-mapq", min_seq_mapq_,
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
}

int App::Exit(const CLI::ParseError& e) { return app_.exit(e); }

void App::Solve() {
    solver_->Import(input_file_path_, min_seq_length_, min_seq_mapq_);
    solver_->Solve(max_ref_coverage_);
    solver_->Export(output_file_path_);
}
