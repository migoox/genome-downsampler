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

bam_api::BamApi::BamApi(const std::filesystem::path& input_filepath, const BamApiConfig& config)
    : input_filepath_(input_filepath) {
    set_min_length_filter(config.min_seq_length);
    set_min_mapq_filter(config.min_mapq);

    if (!config.bed_filepath.empty()) {
        set_amplicon_filter(config.bed_filepath, config.tsv_filepath);
    }
}
bam_api::BamApi::BamApi(const AOSPairedReads& paired_reads_)
    : aos_paired_reads_{paired_reads_}, stored_paired_reads_(PairedReadsType::AOS) {}
bam_api::BamApi::BamApi(const SOAPairedReads& paired_reads_)
    : soa_paired_reads_{paired_reads_}, stored_paired_reads_(PairedReadsType::SOA) {}

void bam_api::BamApi::set_min_length_filter(uint32_t min_length) { min_seq_length_ = min_length; }

void bam_api::BamApi::set_min_mapq_filter(uint32_t min_mapq) { min_mapq_ = min_mapq; }

void bam_api::BamApi::set_amplicon_filter(const std::filesystem::path& bed_filepath,
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
        LOG_WITH_LEVEL(logging::ERROR) << "Error opening .bed file: " << filepath;
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string chrom;
        std::string start;
        std::string end;
        std::string name;

        Index start_idx = 0;
        Index end_idx = 0;

        try {
            start_idx = std::stoull(start);
            end_idx = std::stoull(end);
        } catch (const std::invalid_argument& e) {
            LOG_WITH_LEVEL(logging::ERROR) << "Invalid argument: " << e.what();
            continue;
        } catch (const std::out_of_range& e) {
            LOG_WITH_LEVEL(logging::ERROR) << "Out of range: " << e.what();
            continue;
        }

        std::getline(ss, chrom, '\t');
        std::getline(ss, start, '\t');
        std::getline(ss, end, '\t');
        std::getline(ss, name, '\t');

        if (!chrom.empty() && !start.empty() && !end.empty() && !name.empty()) {
            LOG_WITH_LEVEL(logging::DEBUG) << "Chromosome: " << chrom << ", Start: " << start
                                           << ", End: " << end << ", Name: " << name;
            ret.emplace(name, std::pair(start_idx, end_idx));
        } else {
            LOG_WITH_LEVEL(logging::ERROR) << "Invalid BED line: " << line;
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
        LOG_WITH_LEVEL(logging::ERROR) << "Error opening .tsv file: " << filepath;
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
            LOG_WITH_LEVEL(logging::DEBUG)
                << "Left primer: " << left << ", Right primer: " << right;
            ret.push_back(std::pair(left, right));
        } else {
            LOG_WITH_LEVEL(logging::ERROR) << "Invalid TSV line: " << line;
        }
    }

    file.close();

    return ret;
}

const bam_api::AOSPairedReads& bam_api::BamApi::get_paired_reads_aos() {
    if (stored_paired_reads_ == PairedReadsType::AOS) {
        return aos_paired_reads_;
    }

    if (stored_paired_reads_ == PairedReadsType::SOA) {
        LOG_WITH_LEVEL(logging::INFO)
            << "WARNING: You are reading input file more than once. Consider using only one of "
               "SOA and AOS structures for paired reads.";
    }

    if (input_filepath_.empty()) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "There is no AOS paired reads set nor input filepath specified. BamApi panicked!";
        exit(EXIT_FAILURE);
    }

    read_bam(input_filepath_, aos_paired_reads_);

    return aos_paired_reads_;
}

const bam_api::SOAPairedReads& bam_api::BamApi::get_paired_reads_soa() {
    if (stored_paired_reads_ == PairedReadsType::SOA) {
        return soa_paired_reads_;
    }

    if (stored_paired_reads_ == PairedReadsType::AOS) {
        LOG_WITH_LEVEL(logging::INFO)
            << "WARNING: You are reading input file more than once. Consider using only one of "
               "SOA and AOS structures for paired reads.";
    }

    if (input_filepath_.empty()) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "There is no SOA paired reads set nor input filepath specified. BamApi panicked!";
        exit(EXIT_FAILURE);
    }

    read_bam(input_filepath_, soa_paired_reads_);

    return soa_paired_reads_;
}

const std::vector<bam_api::BAMReadId>& bam_api::BamApi::get_filtered_out_reads() const {
    return filtered_out_reads_;
}

std::vector<bam_api::BAMReadId> bam_api::BamApi::find_pairs(
    const std::vector<BAMReadId>& bam_ids) const {
    std::vector<BAMReadId> paired_bam_ids;
    paired_bam_ids.reserve(bam_ids.size());

    const PairedReads& paired_reads = get_paired_reads();
    std::vector<bool> read_mapped(paired_reads.bam_id_to_read_index.size(), false);

    for (const BAMReadId& bam_id : bam_ids) {
        if (!read_mapped[bam_id]) {
            paired_bam_ids.push_back(bam_id);
            read_mapped[bam_id] = true;
        }

        std::optional<BAMReadId> pair_bam_id = paired_reads.read_pair_map[bam_id];
        if (pair_bam_id && !read_mapped[pair_bam_id.value()]) {
            paired_bam_ids.push_back(pair_bam_id.value());
            read_mapped[pair_bam_id.value()] = true;
        }
    }

    return paired_bam_ids;
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

const bam_api::PairedReads& bam_api::BamApi::get_paired_reads() const {
    if (stored_paired_reads_ == PairedReadsType::SOA) {
        return soa_paired_reads_;
    }

    return aos_paired_reads_;
}

bool bam_api::BamApi::should_be_filtered_out(const Read& r1, const Read& r2) {
    bool ret = !have_min_mapq(r1, r2, min_mapq_) || !have_min_length(r1, r2, min_seq_length_);

    if (amplicon_behaviour_ == AmpliconBehaviour::FILTER) {
        ret = ret || !are_from_single_amplicon(r1, r2, amplicon_set_);
    }

    return ret;
}

bool bam_api::BamApi::have_min_length(const Read& r1, const Read& r2, uint32_t min_length) {
    return (r1.end_ind - r1.start_ind + 1 >= min_length) &&
           (r2.end_ind - r2.start_ind + 1 >= min_length);
}

bool bam_api::BamApi::have_min_mapq(const Read& r1, const Read& r2, uint32_t min_mapq) {
    return (r1.quality >= min_mapq) && (r2.quality >= min_mapq);
}

bool bam_api::BamApi::are_from_single_amplicon(const Read& r1, const Read& r2,
                                               const AmpliconSet& amplicon_set) {
    return amplicon_set.member_includes_both(r1, r2);
}

void bam_api::BamApi::apply_amplicon_inclusion_grading(Read& r1, Read& r2,
                                                       const AmpliconSet& amplicon_set) {
    // TODO(mytkom): think of better grading
    if (are_from_single_amplicon(r1, r2, amplicon_set)) {
        r1.quality += kMaxMAPQ;
        r2.quality += kMaxMAPQ;
    }
}

void bam_api::BamApi::read_bam(const std::filesystem::path& input_filepath,
                               PairedReads& paired_reads) {
    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    bamdata = bam_init1();
    if (!bamdata) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to allocate data memory!";
        std::exit(EXIT_FAILURE);
    }

    // open input file
    infile = sam_open(input_filepath.c_str(), "r");
    if (!infile) {
        LOG_WITH_LEVEL(logging::ERROR) << "Could not open " << input_filepath;
        sam_close(infile);
        std::exit(EXIT_FAILURE);
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to read header from file!";
        sam_close(infile);
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

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        Read current_read(id, bamdata);
        current_qname = bam_get_qname(bamdata);

        // You cannot resize it before because bam index can be not available for input
        paired_reads.read_pair_map.push_back(std::nullopt);
        paired_reads.bam_id_to_read_index.push_back(std::nullopt);

        auto read_one_iterator = read_map.find(current_qname);
        if (read_one_iterator != read_map.end()) {
            Read& r1 = read_one_iterator->second;
            Read& r2 = current_read;

            if (should_be_filtered_out(r1, r2)) {
                id++;
                continue;
            }

            if (amplicon_behaviour_ == AmpliconBehaviour::GRADE) {
                apply_amplicon_inclusion_grading(r1, r2, amplicon_set_);
            }

            if (r2.is_first_read) {
                std::swap(r1, r2);
            }

            ReadIndex first_read_index = paired_reads.get_reads_count();
            paired_reads.push_back(r1);
            paired_reads.push_back(r2);

            paired_reads.bam_id_to_read_index[r1.bam_id] = first_read_index;
            paired_reads.bam_id_to_read_index[r2.bam_id] = first_read_index + 1;

            paired_reads.read_pair_map[r1.bam_id] = r2.bam_id;
            paired_reads.read_pair_map[r2.bam_id] = r1.bam_id;
        } else {
            read_map.insert({ current_qname, current_read });
        }

        id++;
    }

    assert(paired_reads.bam_id_to_read_index.size() == paired_reads.read_pair_map.size());
    assert(paired_reads.bam_id_to_read_index.size() == id);

    for (BAMReadId bam_id = 0; bam_id < paired_reads.bam_id_to_read_index.size(); ++bam_id) {
        // no mapping between bam_id and read_index - this sequence was filtered out
        if (!paired_reads.bam_id_to_read_index[bam_id]) {
            filtered_out_reads_.push_back(bam_id);
        }
    }

    assert((paired_reads.get_reads_count() + filtered_out_reads_.size()) ==
         paired_reads.read_pair_map.size());

    if (ret_r >= 0) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";
    }

    // cleanup
    sam_hdr_destroy(in_samhdr);
    sam_close(infile);
    bam_destroy1(bamdata);
}

uint32_t bam_api::BamApi::write_paired_reads(const std::filesystem::path& output_filepath,
                                             std::vector<BAMReadId>& active_bam_ids) const {
    return write_bam(input_filepath_, output_filepath, active_bam_ids);
}

uint32_t bam_api::BamApi::write_bam_api_filtered_out_reads(
    const std::filesystem::path& output_filepath) {
    return write_bam(input_filepath_, output_filepath, filtered_out_reads_);
}

uint32_t bam_api::BamApi::write_bam(const std::filesystem::path& input_filepath,
                                    const std::filesystem::path& output_filepath,
                                    std::vector<BAMReadId>& bam_ids) {
    LOG_WITH_LEVEL(logging::DEBUG) << "BamApi: writing " << bam_ids.size() << " reads to "
                                   << output_filepath << " on the basis of " << input_filepath;

    sam_hdr_t* in_samhdr = NULL;
    samFile* infile = NULL;
    samFile* outfile = NULL;
    int ret_r = 0;
    bam1_t* bamdata = NULL;

    bamdata = bam_init1();
    if (!bamdata) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to allocate data memory!";
        std::exit(EXIT_FAILURE);
    }

    // open input file
    infile = sam_open(input_filepath.c_str(), "r");
    if (!infile) {
        LOG_WITH_LEVEL(logging::ERROR) << "Could not open " << input_filepath;
        std::exit(EXIT_FAILURE);
    }

    // open output file
    std::string open_mode = output_filepath.extension() == ".bam" ? "wb" : "w";
    outfile = sam_open(output_filepath.c_str(), open_mode.c_str());
    if (!outfile) {
        LOG_WITH_LEVEL(logging::ERROR) << "Could not open " << output_filepath;
        std::exit(EXIT_FAILURE);
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to read header from file!";
        std::exit(EXIT_FAILURE);
    }

    if (sam_hdr_write(outfile, in_samhdr) < 0) {
        LOG_WITH_LEVEL(logging::ERROR) << "Can't write header to bam file: " << output_filepath;
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
                LOG_WITH_LEVEL(logging::ERROR)
                    << "Can't write line to bam file: " << output_filepath;
                std::exit(EXIT_FAILURE);
            }

            reads_written++;
            current_read_i++;
        }

        id++;
    }

    if (current_read_i != bam_ids.end() && ret_r >= 0)
        LOG_WITH_LEVEL(logging::ERROR)
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

    LOG_WITH_LEVEL(logging::DEBUG) << "BamApi: " << reads_written << " reads have been written to "
                                   << output_filepath << " on the basis of " << input_filepath;

    return reads_written;
}
