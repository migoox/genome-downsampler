#include "../include/bam-api/bam_api.hpp"

#include <htslib/hts.h>
#include <htslib/regidx.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "../include/bam-api/bam_paired_reads.hpp"
#include "bam-api/amplicons.hpp"
#include "logging/log.hpp"

void bam_api::BamApi::add_min_length_filter(uint32_t min_length) { min_seq_length_ = min_length; }

void bam_api::BamApi::add_min_mapq_filter(uint32_t min_mapq) { min_mapq_ = min_mapq; }

void bam_api::BamApi::add_amplicon_filter(const std::filesystem::path& bed_filepath,
                                          const std::filesystem::path& tsv_filepath) {
    std::map<std::string, std::pair<bam_api::Index, bam_api::Index>> primer_map;
    primer_map = process_bed_file(bed_filepath);

    if (!tsv_filepath.empty()) {
        std::vector<std::pair<std::string, std::string>> pairs;
        pairs = process_tsv_file(tsv_filepath);

        for (auto& [left, right] : pairs) {
            auto& left_primer = primer_map[left];
            auto& right_primer = primer_map[right];

            if (left_primer.first > right_primer.first) {
                std::swap(left_primer, right_primer);
            }

            amplicon_set_.amplicons.push_back(Amplicon(left_primer.first, right_primer.second));
        }
    } else {
        auto it = primer_map.begin();

        while (it != primer_map.end()) {
            auto& left_primer = it->second;
            it++;
            auto& right_primer = it->second;

            if (left_primer.first > right_primer.first) {
                std::swap(left_primer, right_primer);
            }

            amplicon_set_.amplicons.push_back(Amplicon(left_primer.first, right_primer.second));

            it++;
        }
    }
}

void bam_api::BamApi::set_amplicon_behaviour(AmpliconBehaviour amplicon_behaviour) {
    amplicon_behaviour_ = amplicon_behaviour;
}

std::map<std::string, std::pair<bam_api::Index, bam_api::Index>> bam_api::BamApi::process_bed_file(
    const std::filesystem::path& filepath) {
    std::map<std::string, std::pair<bam_api::Index, bam_api::Index>> ret;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        LOG_WITH_LEVEL(logging::kError) << "Error opening .bed file: " << filepath;
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string chrom;
        std::string start;
        std::string end;
        std::string name;

        std::getline(ss, chrom, '\t');
        std::getline(ss, start, '\t');
        std::getline(ss, end, '\t');
        std::getline(ss, name, '\t');

        if (!chrom.empty() && !start.empty() && !end.empty() && !name.empty()) {
            LOG_WITH_LEVEL(logging::kDebug) << "Chromosome: " << chrom << ", Start: " << start
                                            << ", End: " << end << ", Name: " << name;
            ret.emplace(name, std::pair(start, end));
        } else {
            LOG_WITH_LEVEL(logging::kError) << "Invalid BED line: " << line;
        }
    }

    file.close();

    return ret;
}

std::vector<std::pair<std::string, std::string>> bam_api::BamApi::process_tsv_file(
    const std::filesystem::path& filepath) {
    std::vector<std::pair<std::string, std::string>> ret;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        LOG_WITH_LEVEL(logging::kError) << "Error opening .tsv file: " << filepath;
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string left;
        std::string right;

        std::getline(ss, left, '\t');
        std::getline(ss, right, '\t');

        if (!left.empty() && !right.empty()) {
            LOG_WITH_LEVEL(logging::kDebug)
                << "Left primer: " << left << ", Right primer: " << right;
            ret.push_back(std::pair(left, right));
        } else {
            LOG_WITH_LEVEL(logging::kError) << "Invalid TSV line: " << line;
        }
    }

    file.close();

    return ret;
}

const bam_api::AOSPairedReads& bam_api::BamApi::get_paired_reads_aos() {
    if (stored_paired_reads_ == AOS) {
        return aos_paired_reads_;
    }

    if (stored_paired_reads_ == SOA) {
        LOG_WITH_LEVEL(logging::kInfo)
            << "WARNING: You are reading input file more than once. Consider using only one of "
               "SOA and AOS structures for paired reads.";
    }

    if (input_filepath_.empty()) {
        LOG_WITH_LEVEL(logging::kError)
            << "There is no AOS paired reads set nor input filepath specified. BamApi panicked!";
        exit(EXIT_FAILURE);
    }

    // TODO(mytkom): read_bam execution

    return aos_paired_reads_;
}

const bam_api::SOAPairedReads& bam_api::BamApi::get_paired_reads_soa() {
    if (stored_paired_reads_ == SOA) {
        return soa_paired_reads_;
    }

    if (stored_paired_reads_ == AOS) {
        LOG_WITH_LEVEL(logging::kInfo)
            << "WARNING: You are reading input file more than once. Consider using only one of "
               "SOA and AOS structures for paired reads.";
    }

    if (input_filepath_.empty()) {
        LOG_WITH_LEVEL(logging::kError)
            << "There is no SOA paired reads set nor input filepath specified. BamApi panicked!";
        exit(EXIT_FAILURE);
    }

    // TODO(mytkom): read_bam execution

    return soa_paired_reads_;
}

const std::vector<bam_api::Read>& bam_api::BamApi::get_filtered_out_reads() const {
    return filtered_out_reads_;
}

std::vector<uint32_t> bam_api::BamApi::find_input_cover() {
    const bam_api::PairedReads& paired_reads = get_paired_reads();
    std::vector<uint32_t> result(paired_reads.ref_genome_length, 0);
    for (bam_api::ReadIndex i = 0; i < paired_reads.get_reads_count(); ++i) {
        Read curr_read = paired_reads.get_read_by_index(i);
        for (bam_api::Index j = curr_read.start_ind; j <= curr_read.end_ind; ++j) {
            result[j]++;
        }
    }

    return result;
}

std::vector<uint32_t> bam_api::BamApi::find_filtered_cover(
    const std::vector<bam_api::BAMReadId>& bam_ids) {
    const bam_api::PairedReads& paired_reads = get_paired_reads();
    std::vector<uint32_t> result(paired_reads.ref_genome_length, 0);

    for (const auto bam_id : bam_ids) {
        const auto& read = paired_reads.get_read_by_bam_id(bam_id);
        for (bam_api::Index i = read.start_ind; i <= read.end_ind; ++i) {
            result[i]++;
        }
    }

    return result;
}

const bam_api::PairedReads& bam_api::BamApi::get_paired_reads() {
    if (stored_paired_reads_ == SOA) {
        return soa_paired_reads_;
    }

    return aos_paired_reads_;
}

void bam_api::BamApi::read_bam(const std::filesystem::path& input_filepath,
                               PairedReads& paired_reads, uint32_t min_seq_length,
                               uint32_t min_mapq,
                               const std::filesystem::path& bed_amplicon_bounds_filepath,
                               const std::filesystem::path& tsv_amplicon_pairs_filepath) {
    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    samFile* filtered_out_file = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    bamdata = bam_init1();
    if (!bamdata) {
        LOG_WITH_LEVEL(logging::kError) << "Failed to allocate data memory!";
        std::exit(EXIT_FAILURE);
    }

    // open input file
    infile = sam_open(input_filepath.c_str(), "r");
    if (!infile) {
        LOG_WITH_LEVEL(logging::kError) << "Could not open " << input_filepath;
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }

    auto filtered_out_path = input_filepath;
    filtered_out_path.replace_filename("filtered_out_sequences.bam");
    filtered_out_file = sam_open(filtered_out_path.c_str(), "wb");
    if (!filtered_out_file) {
        LOG_WITH_LEVEL(logging::kError) << "Could not open " << filtered_out_path;
        sam_close(filtered_out_file);
        std::exit(EXIT_FAILURE);
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG_WITH_LEVEL(logging::kError) << "Failed to read header from file!";
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }

    if (sam_hdr_write(filtered_out_file, in_samhdr) < 0) {
        LOG_WITH_LEVEL(logging::kError) << "Can't write header to bam file: " << filtered_out_path;
        std::exit(EXIT_FAILURE);
    }

    // Reserve memory for reads if index is present
    hts_idx_t* idx = sam_index_load(infile, input_filepath.c_str());
    if (idx) {
        uint64_t mapped_seq_c, unmapped_seq_c;
        hts_idx_get_stat(idx, 0, &mapped_seq_c, &unmapped_seq_c);
        paired_reads.reserve(mapped_seq_c);
        hts_idx_destroy(idx);
    }

    paired_reads.ref_genome_length = in_samhdr->target_len[0];

    BAMReadId id = 0;
    std::map<std::string, Read> read_map;
    std::string current_qname;
    hts_pos_t rlen = 0;

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        rlen = bam_cigar2rlen(static_cast<int32_t>(bamdata->core.n_cigar), bam_get_cigar(bamdata));
        Read current_read(id, static_cast<Index>(bamdata->core.pos),
                          static_cast<Index>(bamdata->core.pos + rlen - 1), bamdata->core.qual,
                          static_cast<bool>(bamdata->core.flag & BAM_FREAD1));
        current_qname = bam_get_qname(bamdata);  // NOLINT

        auto read_one_iterator = read_map.find(current_qname);
        if (read_one_iterator != read_map.end()) {
            ReadIndex first_read_index = paired_reads.get_reads_count();

            // if first read of pair is under index i, second is under i+1
            if (current_read.is_first_read) {
                paired_reads.push_back(current_read);
                paired_reads.push_back(read_one_iterator->second);

                paired_reads.bam_id_to_read_index.push_back(first_read_index);
                paired_reads.bam_id_to_read_index[read_one_iterator->second.bam_id] =
                    first_read_index + 1;
            } else {
                paired_reads.push_back(read_one_iterator->second);
                paired_reads.push_back(current_read);

                paired_reads.bam_id_to_read_index[read_one_iterator->second.bam_id] =
                    first_read_index;
                paired_reads.bam_id_to_read_index.push_back(first_read_index + 1);
            }

            paired_reads.read_pair_map[read_one_iterator->second.bam_id] = current_read.bam_id;
            paired_reads.read_pair_map.push_back(read_one_iterator->second.bam_id);

        } else {
            if (bamdata->core.l_qseq >= min_seq_length && bamdata->core.qual >= min_mapq) {
                read_map.insert({current_qname, current_read});
            }
            paired_reads.read_pair_map.push_back(std::nullopt);
            paired_reads.bam_id_to_read_index.push_back(std::nullopt);
        }

        id++;
    }

    assert(paired_reads.bam_id_to_read_index.size() == paired_reads.read_pair_map.size());

    if (ret_r >= 0) {
        LOG_WITH_LEVEL(logging::kError)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";
    }

    sam_hdr_destroy(in_samhdr);
    // Reopen infile to read iterate through it second time
    sam_close(infile);
    infile = sam_open(input_filepath.c_str(), "r");
    if (!infile) {
        LOG_WITH_LEVEL(logging::kError) << "Could not open " << input_filepath;
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }
    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG_WITH_LEVEL(logging::kError) << "Failed to read header from file!";
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }

    id = 0;
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        if (!paired_reads.read_pair_map[id]) {
            if (sam_write1(filtered_out_file, in_samhdr, bamdata) < 0) {
                LOG_WITH_LEVEL(logging::kError)
                    << "Can't write line to bam file: " << filtered_out_path;
                std::exit(EXIT_FAILURE);
            }
        }

        id++;
    }

    if (ret_r >= 0) {
        LOG_WITH_LEVEL(logging::kError)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";
    } else {
        LOG_WITH_LEVEL(logging::kDebug) << "Read bam file have been read correctly";
    }

    // cleanup

    sam_hdr_destroy(in_samhdr);
    sam_close(infile);
    sam_close(filtered_out_file);
    bam_destroy1(bamdata);
}


uint32_t bam_api::BamApi::write_bam(const std::filesystem::path& input_filepath,
                                    const std::filesystem::path& output_filepath,
                                    std::vector<BAMReadId>& bam_ids) {
    LOG_WITH_LEVEL(logging::kDebug) << "bam-api writing" << bam_ids.size() << " reads to "
                                    << output_filepath << " on the basis of " << input_filepath;

    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    samFile* outfile = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    bamdata = bam_init1();
    if (!bamdata) {
        LOG_WITH_LEVEL(logging::kError) << "Failed to allocate data memory!";
        std::exit(EXIT_FAILURE);
    }

    // open input file
    infile = sam_open(input_filepath.c_str(), "r");
    if (!infile) {
        LOG_WITH_LEVEL(logging::kError) << "Could not open " << input_filepath;
        std::exit(EXIT_FAILURE);
    }

    // open output file
    std::string open_mode = output_filepath.extension() == ".bam" ? "wb" : "w";
    outfile = sam_open(output_filepath.c_str(), open_mode.c_str());
    if (!outfile) {
        LOG_WITH_LEVEL(logging::kError) << "Could not open " << output_filepath;
        std::exit(EXIT_FAILURE);
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG_WITH_LEVEL(logging::kError) << "Failed to read header from file!";
        std::exit(EXIT_FAILURE);
    }

    if (sam_hdr_write(outfile, in_samhdr) < 0) {
        LOG_WITH_LEVEL(logging::kError) << "Can't write header to bam file: " << output_filepath;
        std::exit(EXIT_FAILURE);
    }

    BAMReadId id = 0;
    std::sort(bam_ids.begin(), bam_ids.end());
    auto current_read_i = bam_ids.begin();
    uint32_t reads_written = 0;

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0 &&
           current_read_i != bam_ids.end()) {
        if (id == *current_read_i) {
            if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
                LOG_WITH_LEVEL(logging::kError)
                    << "Can't write line to bam file: " << output_filepath;
                std::exit(EXIT_FAILURE);
            }

            reads_written++;
            current_read_i++;
        }

        id++;
    }

    if (current_read_i != bam_ids.end() && ret_r >= 0)
        LOG_WITH_LEVEL(logging::kError)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";

    // cleanup
    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (infile) {
        sam_close(infile);
    }
    if (outfile) {
        sam_close(outfile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }

    LOG_WITH_LEVEL(logging::kDebug) << "bam-api " << reads_written << " reads have been written to "
                                    << output_filepath << " on the basis of " << input_filepath;

    return reads_written;
}
