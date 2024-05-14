#pragma once

#include <CLI/CLI.hpp>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"

class Solver {
   public:
    virtual void Import(const std::filesystem::path& input,
                        uint32_t min_seq_length, uint32_t min_seq_mapq) = 0;
    virtual void Solve() = 0;
    virtual void Export(const std::filesystem::path& output) = 0;
};

class ExampleSolver : public Solver {
   public:
    void Import(const std::filesystem::path& input, uint32_t min_seq_length,
                uint32_t min_seq_mapq) override {
        input_ = input;
        std::cout << "Import, min_len: " << min_seq_length
                  << ", min_mapq: " << min_seq_mapq << ", input: " << input
                  << std::endl;

        paired_reads_ =
            bam_api::BamApi::read_bam_soa(input, min_seq_length, min_seq_mapq);

        std::cout << paired_reads_.ids.size() << " sequences has been imported!"
                  << std::endl;
    }

    void Solve() override {
        std::cout << "Solve" << std::endl;

        solution_ = std::vector<bam_api::ReadIndex>(
            paired_reads_.ids.end() - 1000, paired_reads_.ids.end() - 2);

        std::cout << "Solution have " << solution_.size() << " sequences!"
                  << std::endl;
    }

    void Export(const std::filesystem::path& output) override {
        std::cout << "Output: " << output << std::endl;
        uint32_t reads_written;
        if (output.extension() == ".bam") {
            reads_written =
                bam_api::BamApi::write_sam(input_, output, solution_, true);
        } else {
            reads_written =
                bam_api::BamApi::write_sam(input_, output, solution_, false);
        }
        std::cout << reads_written << " sequences has been written to file "
                  << output << std::endl;
    }

   private:
    bam_api::SOAPairedReads paired_reads_;
    std::vector<bam_api::ReadIndex> solution_;
    std::filesystem::path input_;
};

class App {
    static constexpr uint32_t kDefaultMinSeqLength = 90;
    static constexpr uint32_t kDefaultMinSeqQmap = 30;

   public:
    App() {
        app_.add_option("-i,--input", input_file_path_,
                        ".bam input file path. Required option.")
            ->required()
            ->check(CLI::ExistingFile);

        std::vector<std::string> possible_algorithms;
        for (const auto& mapping : algorithms_map_) {
            possible_algorithms.push_back(mapping.first);
        }
        app_.add_option_function<std::string>(
                "-a,--algorithm",
                [this](const std::string algorithm_name) {
                    if (algorithms_map_.find(algorithm_name) !=
                        algorithms_map_.end()) {
                        algorithm_ = algorithms_map_[algorithm_name];
                    } else {
                        std::cerr << "Algorithm not found: " << algorithm_name
                                  << std::endl;
                    }
                },
                "Algorithm to use.")
            ->required()
            ->check(CLI::IsMember(possible_algorithms));

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

        app_.add_option(
                "-q,--min-mapq", min_seq_mapq_,
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

    int Exit(const CLI::ParseError& e) { return app_.exit(e); }

    void Solve() {
        algorithm_->Import(input_file_path_, min_seq_length_, min_seq_mapq_);
        algorithm_->Solve();
        algorithm_->Export(output_file_path_);
    }

   private:
    CLI::App app_;
    bool verbose_mode_ = false;
    uint32_t min_seq_length_ = kDefaultMinSeqLength;
    uint32_t min_seq_mapq_ = kDefaultMinSeqQmap;
    Solver* algorithm_;
    std::filesystem::path input_file_path_;
    std::filesystem::path output_file_path_;
    std::filesystem::path csv_historical_runs_file_path_;

    // Add new solver by adding a class and new entry to map below
    std::map<std::string, Solver*> algorithms_map_ = {
        {"example", new ExampleSolver()}};
};
