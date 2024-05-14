#pragma once

#include <CLI/CLI.hpp>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"

class App {
    static constexpr uint32_t kDefaultMinSeqLength = 90;
    static constexpr uint32_t kDefaultMinSeqQmap = 30;

    CLI::App app_;
    bool verbose_mode_ = false;
    uint32_t min_seq_length_ = kDefaultMinSeqLength;
    uint32_t min_seq_mapq_ = kDefaultMinSeqQmap;
    std::filesystem::path input_file_path_;
    std::filesystem::path output_file_path_;
    std::filesystem::path csv_historical_runs_file_path_ =
        "./historical_runs.csv";

   public
    :App() {
        app_.add_option("-i,--input", input_file_path_,
                       ".bam input file path. Required option.")
            ->required()
            ->check(CLI::ExistingFile);

        app_.add_option("-o,--output", output_file_path_,
                       ".bam output file path. Default is \"output.bam\" in "
                       "input's directory.");

        app_.add_option("-c,--csv-history", csv_historical_runs_file_path_,
                       ".csv historical runs data file path. Default is "
                       "\"historical_runs.csv\" in current directory.");

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

    void Parse(int argc, char** argv) {
        app_.parse(argc, argv);

        if (output_file_path_.empty()) {
            output_file_path_ = input_file_path_;
            output_file_path_.replace_filename("output.bam");
        }
    }

    int Exit(const CLI::ParseError &e) {
      return app_.exit(e);
    }

    void Solve() {
        auto ret_soa = bam_api::BamApi::read_bam_soa(
            input_file_path_, min_seq_length_, min_seq_mapq_);
        std::vector<bam_api::ReadIndex> temp(ret_soa.ids.end() - 1000,
                                             ret_soa.ids.end() - 2);
        bam_api::BamApi::write_sam(input_file_path_, output_file_path_, temp);
    }
};
