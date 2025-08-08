#include "bam-api/bam_api.hpp"

#include <htslib/hts.h>
#include <htslib/regidx.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "bam-api/paired_reads.hpp"
#include "bam-api/read.hpp"
#include "logging/log.hpp"

#define RELEASE_TPOOL(X)                  \
    {                                     \
        hts_tpool* ptr = (hts_tpool*)(X); \
        if (ptr) {                        \
            hts_tpool_destroy(ptr);       \
        }                                 \
    }

bam_api::BamApi::BamApi(const std::filesystem::path& input_filepath, const BamApiConfig& config)
    : input_filepath_(input_filepath),
      min_seq_length_(config.min_seq_length),
      amp_overflow_(config.amp_overflow),
      min_mapq_(config.min_mapq),
      min_alignment_(config.min_alignment),
      hts_thread_count_(config.hts_thread_count) {
    if (!config.bed_filepath.empty()) {
        set_amplicon_filter(config.bed_filepath, config.tsv_filepath);
        amplicon_behaviour_ = config.amplicon_behaviour;
    } else {
        amplicon_behaviour_ = AmpliconBehaviour::IGNORE;
    }
}
bam_api::BamApi::BamApi(const AOSPairedReads& paired_reads)
    : aos_paired_reads_{paired_reads}, is_aos_loaded_(true) {}
bam_api::BamApi::BamApi(const SOAPairedReads& paired_reads)
    : soa_paired_reads_{paired_reads}, is_soa_loaded_(true) {}

void bam_api::BamApi::set_amplicon_filter(const std::filesystem::path& bed_filepath,
                                          const std::filesystem::path& tsv_filepath) {
    auto start = std::chrono::high_resolution_clock::now();

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

            amplicon_set_.amplicons.push_back(
                Amplicon(left_primer.first, right_primer.second, amp_overflow_));
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

            amplicon_set_.amplicons.push_back(
                Amplicon(left_primer.first, right_primer.second, amp_overflow_));

            it++;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    LOG_WITH_LEVEL(logging::DEBUG) << "set_amplicon_filter took " << elapsed.count() << " seconds";
}

void bam_api::BamApi::set_amplicon_behaviour(AmpliconBehaviour amplicon_behaviour) {
    amplicon_behaviour_ = amplicon_behaviour;
}

std::map<std::string, std::pair<bam_api::Index, bam_api::Index>> bam_api::BamApi::process_bed_file(
    const std::filesystem::path& filepath) {
    LOG_WITH_LEVEL(logging::INFO) << "Reading " << filepath.filename() << " file...";

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

        std::getline(ss, chrom, '\t');
        std::getline(ss, start, '\t');
        std::getline(ss, end, '\t');
        std::getline(ss, name, '\t');

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

        if (!chrom.empty() && !start.empty() && !end.empty() && !name.empty()) {
            ret.emplace(name, std::pair(start_idx, end_idx));
        } else {
            LOG_WITH_LEVEL(logging::ERROR) << "Invalid BED line: " << line;
        }
    }

    file.close();

    LOG_WITH_LEVEL(logging::DEBUG) << ret.size() << " primers have been read";

    return ret;
}

std::vector<std::pair<std::string, std::string>> bam_api::BamApi::process_tsv_file(
    const std::filesystem::path& filepath) {
    LOG_WITH_LEVEL(logging::INFO) << "Reading " << filepath.filename() << " file...";

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
            ret.push_back(std::pair(left, right));
        } else {
            LOG_WITH_LEVEL(logging::ERROR) << "Invalid TSV line: " << line;
        }
    }

    file.close();

    LOG_WITH_LEVEL(logging::DEBUG) << ret.size() << " pairs of primers have been read";

    return ret;
}

const bam_api::AOSPairedReads& bam_api::BamApi::get_paired_reads_aos() {
    if (is_aos_loaded_) {
        return aos_paired_reads_;
    }

    if (is_soa_loaded_) {
        aos_paired_reads_.from(soa_paired_reads_);
        is_aos_loaded_ = true;
        return aos_paired_reads_;
    }

    if (input_filepath_.empty()) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "There is no AOS paired reads set nor input filepath specified. BamApi panicked!";
        exit(EXIT_FAILURE);
    }

    read_bam(input_filepath_, aos_paired_reads_);
    is_aos_loaded_ = true;

    return aos_paired_reads_;
}

const bam_api::SOAPairedReads& bam_api::BamApi::get_paired_reads_soa() {
    if (is_soa_loaded_) {
        return soa_paired_reads_;
    }

    if (is_aos_loaded_) {
        soa_paired_reads_.from(aos_paired_reads_);
        is_soa_loaded_ = true;
        return soa_paired_reads_;
    }

    if (input_filepath_.empty()) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "There is no SOA paired reads set nor input filepath specified. BamApi panicked!";
        exit(EXIT_FAILURE);
    }

    read_bam(input_filepath_, soa_paired_reads_);
    is_soa_loaded_ = true;

    return soa_paired_reads_;
}

const std::vector<bam_api::BAMReadId>& bam_api::BamApi::get_filtered_out_reads() const {
    return filtered_out_reads_;
}

std::vector<bam_api::ReadIndex> bam_api::BamApi::find_pairs(
    const std::vector<ReadIndex>& ids) const {
    if (!is_paired_) {
        LOG_WITH_LEVEL(logging::INFO) << "Unpaired data, no need to find pairs.";
        return ids;
    }

    LOG_WITH_LEVEL(logging::INFO) << "Finding paired reads for solution...";
    LOG_WITH_LEVEL(logging::DEBUG) << "Unpaired solution have " << ids.size() << " reads";

    auto start = std::chrono::high_resolution_clock::now();

    const PairedReads& paired_reads = get_paired_reads();

    std::vector<ReadIndex> paired_ids;
    paired_ids.reserve(paired_reads.get_reads_count());

    std::vector<bool> read_mapped(paired_reads.get_reads_count(), false);

    for (const BAMReadId& id : ids) {
        if (!read_mapped[id]) {
            paired_ids.push_back(id);
            read_mapped[id] = true;
        }

        ReadIndex pair_id = paired_reads.get_read_by_index(id).is_first_read ? id + 1 : id - 1;
        if (!read_mapped[pair_id]) {
            paired_ids.push_back(pair_id);
            read_mapped[pair_id] = true;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    LOG_WITH_LEVEL(logging::DEBUG) << "find_pairs took " << elapsed.count() << " seconds";

    LOG_WITH_LEVEL(logging::DEBUG) << "Paired solution have " << paired_ids.size() << " reads";

    return paired_ids;
}

std::vector<uint32_t> bam_api::BamApi::find_input_cover() const {
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
    const std::vector<bam_api::ReadIndex>& ids) const {
    const bam_api::PairedReads& paired_reads = get_paired_reads();
    std::vector<uint32_t> result(paired_reads.ref_genome_length, 0);

    for (const auto id : ids) {
        const auto& read = paired_reads.get_read_by_index(id);
        for (bam_api::Index i = read.start_ind; i <= read.end_ind; ++i) {
            result[i]++;
        }
    }

    return result;
}

const bam_api::PairedReads& bam_api::BamApi::get_paired_reads() const {
    if (is_soa_loaded_) {
        return soa_paired_reads_;
    }

    return aos_paired_reads_;
}

bool bam_api::BamApi::should_be_filtered_out(const Read& r, bool is_paired) const {
    bool ret = false;

    ret = ret || !have_min_mapq(r);
    ret = ret || !have_min_length(r);
    if (!is_paired) ret = ret || !have_min_alignment(r);
    ret = ret || r.has_sa;

    if (amplicon_behaviour_ == AmpliconBehaviour::FILTER) {
        ret = ret || !is_from_single_amplicon(r);
    }

    return ret;
}

bool bam_api::BamApi::should_be_filtered_out(const Read& r1, const Read& r2) const {
    bool ret = should_be_filtered_out(r1, true) || should_be_filtered_out(r2, true);

    if (amplicon_behaviour_ == AmpliconBehaviour::FILTER) {
        ret = ret || !are_from_single_amplicon(r1, r2);
    }

    return ret;
}

bool bam_api::BamApi::have_min_alignment(const Read& r) const {
    return (static_cast<float>(r.as) / static_cast<float>(r.seq_length)) >= min_alignment_;
}

bool bam_api::BamApi::have_min_alignment(const Read& r1, const Read& r2) const {
    return (static_cast<float>(r1.as) / static_cast<float>(r1.seq_length + r2.seq_length)) >=
           min_alignment_;
}

bool bam_api::BamApi::have_min_length(const Read& r) const {
    return r.seq_length >= min_seq_length_;
}

bool bam_api::BamApi::have_min_mapq(const Read& r) const { return r.mapq >= min_mapq_; }

bool bam_api::BamApi::are_from_single_amplicon(const Read& r1, const Read& r2) const {
    return amplicon_set_.member_includes_both(r1, r2);
}

bool bam_api::BamApi::is_from_single_amplicon(const Read& r) const {
    return amplicon_set_.any_includes(r);
}

void bam_api::BamApi::apply_amplicon_inclusion_grading(
    bam_api::PairedReads& paired_reads, std::vector<bool>& is_in_single_amplicon) const {
    LOG_WITH_LEVEL(logging::DEBUG)
        << "Grading: min_mapq: " << min_imported_mapq_ << ", max_mapq: " << max_imported_mapq_;
    if (max_imported_mapq_ > 0 && min_imported_mapq_ < UINT32_MAX) {
        for (ReadIndex i = 0; i < paired_reads.get_reads_count(); ++i) {
            ReadQuality quality = paired_reads.get_quality(i);
            quality -= min_imported_mapq_;
            if (is_in_single_amplicon[i]) {
                quality += max_imported_mapq_ - min_imported_mapq_;
            }
            paired_reads.set_quality(i, quality);
        }
    }
}

void bam_api::BamApi::analyse_mapq(const Read& r) {
    if (r.quality < min_imported_mapq_) {
        min_imported_mapq_ = r.mapq;
    }

    if (r.quality > max_imported_mapq_) {
        max_imported_mapq_ = r.mapq;
    }
}

void bam_api::BamApi::read_bam(const std::filesystem::path& input_filepath,
                               PairedReads& paired_reads) {
    LOG_WITH_LEVEL(logging::INFO) << "Reading " << input_filepath.filename() << " input file...";

    LOG_WITH_LEVEL(logging::DEBUG) << "BamApi: reading " << input_filepath;

    auto start = std::chrono::high_resolution_clock::now();

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

    htsThreadPool tpool = {NULL, 0};
    if (hts_thread_count_ > 1) {
        if (!(tpool.pool = hts_tpool_init(hts_thread_count_))) {
            hts_log_info("Could not initialize thread pool!");
        }

        // set threads
        if (hts_set_thread_pool(infile, &tpool) < 0) {
            LOG_WITH_LEVEL(logging::ERROR) << "Cannot set threads for writing " << input_filepath;
            std::exit(EXIT_FAILURE);
        }
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to read header from file!";
        sam_close(infile);
        RELEASE_TPOOL(tpool.pool);
        std::exit(EXIT_FAILURE);
    }

    std::vector<bool> is_accepted;
    std::vector<bool> is_in_single_amplicon;

    // Reserve memory for reads if index is present
    hts_idx_t* idx = sam_index_load(infile, input_filepath.c_str());
    if (idx) {
        uint64_t mapped_seq_c, unmapped_seq_c;
        hts_idx_get_stat(idx, 0, &mapped_seq_c, &unmapped_seq_c);
        paired_reads.reserve(mapped_seq_c);
        is_accepted.reserve(mapped_seq_c);
        is_in_single_amplicon.reserve(mapped_seq_c);
        hts_idx_destroy(idx);
    }

    paired_reads.ref_genome_length = in_samhdr->target_len[0];

    BAMReadId id = 0;
    std::map<std::string, Read> read_map;
    std::string current_qname;

    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        Read current_read(id, bamdata);
        current_qname = bam_get_qname(bamdata);

        is_accepted.push_back(false);

        auto read_one_iterator = read_map.find(current_qname);

        // Try to find paired-end reads on the fly
        if (read_one_iterator != read_map.end()) {
            Read& r1 = read_one_iterator->second;
            Read& r2 = current_read;

            if (should_be_filtered_out(r1, r2)) {
                id++;
                continue;
            }

            if (amplicon_behaviour_ == AmpliconBehaviour::GRADE) {
                analyse_mapq(r1);
                analyse_mapq(r2);
                if (are_from_single_amplicon(r1, r2)) {
                    is_in_single_amplicon.push_back(true);
                    is_in_single_amplicon.push_back(true);
                } else {
                    is_in_single_amplicon.push_back(false);
                    is_in_single_amplicon.push_back(false);
                }
            }

            if (r2.is_first_read) {
                std::swap(r1, r2);
            }

            paired_reads.push_back(r1);
            paired_reads.push_back(r2);

            is_accepted[r1.bam_id] = true;
            is_accepted[r2.bam_id] = true;
        } else {
            read_map.insert({current_qname, current_read});
        }

        id++;
    }

    // If not paired handle it too
    if (paired_reads.get_reads_count() == 0) {
        is_paired_ = false;

        for (auto& [qname, read] : read_map) {
            if (should_be_filtered_out(read, false)) {
                continue;
            }

            if (amplicon_behaviour_ == AmpliconBehaviour::GRADE) {
                analyse_mapq(read);
                if (is_from_single_amplicon(read)) {
                    is_in_single_amplicon.push_back(true);
                } else {
                    is_in_single_amplicon.push_back(false);
                }
            }

            paired_reads.push_back(read);
            is_accepted[read.bam_id] = true;
        }
    }

    LOG_WITH_LEVEL(logging::INFO) << "After reading: " << read_map.size();

    assert(is_accepted.size() == id);

    for (BAMReadId bam_id = 0; bam_id < is_accepted.size(); ++bam_id) {
        if (!is_accepted[bam_id]) {
            filtered_out_reads_.push_back(bam_id);
        }
    }

    assert((paired_reads.get_reads_count() + filtered_out_reads_.size()) == is_accepted.size());

    if (amplicon_behaviour_ == AmpliconBehaviour::GRADE) {
        apply_amplicon_inclusion_grading(paired_reads, is_in_single_amplicon);
    }

    if (ret_r >= 0) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";
    }

    // cleanup
    sam_hdr_destroy(in_samhdr);
    sam_close(infile);
    bam_destroy1(bamdata);
    RELEASE_TPOOL(tpool.pool);

    auto end = std::chrono::high_resolution_clock::now();

    LOG_WITH_LEVEL(logging::DEBUG) << "BamApi: " << id << " reads have been read";
    LOG_WITH_LEVEL(logging::DEBUG)
        << "BamApi: " << paired_reads.get_reads_count() << " reads have been imported";
    LOG_WITH_LEVEL(logging::DEBUG) << "BamApi: " << filtered_out_reads_.size()
                                   << " reads have been filtered out during preprocessing";

    std::chrono::duration<double> elapsed = end - start;
    LOG_WITH_LEVEL(logging::DEBUG) << "read_bam took " << elapsed.count() << " seconds";
}

uint32_t bam_api::BamApi::write_paired_reads(const std::filesystem::path& output_filepath,
                                             std::vector<ReadIndex>& active_ids) const {
    LOG_WITH_LEVEL(logging::INFO) << "Writing solution of size " << active_ids.size() << " reads "
                                  << output_filepath.filename() << "...";

    const PairedReads& paired_reads = get_paired_reads();

    std::vector<BAMReadId> active_bam_ids;
    active_bam_ids.reserve(active_ids.size());

    for (const auto& id : active_ids) {
        active_bam_ids.push_back(paired_reads.get_read_by_index(id).bam_id);
    }

    return write_bam(input_filepath_, output_filepath, active_bam_ids, hts_thread_count_);
}

uint32_t bam_api::BamApi::write_bam_api_filtered_out_reads(
    const std::filesystem::path& output_filepath) {
    LOG_WITH_LEVEL(logging::INFO) << "Writing " << filtered_out_reads_.size()
                                  << " preprocessing filtered out reads to "
                                  << output_filepath.filename() << "...";
    return write_bam(input_filepath_, output_filepath, filtered_out_reads_, hts_thread_count_);
}

uint32_t bam_api::BamApi::write_bam(const std::filesystem::path& input_filepath,
                                    const std::filesystem::path& output_filepath,
                                    std::vector<BAMReadId>& bam_ids, uint32_t hts_thread_count) {
    LOG_WITH_LEVEL(logging::DEBUG) << "BamApi: writing " << bam_ids.size() << " reads to "
                                   << output_filepath << " on the basis of " << input_filepath;

    auto start = std::chrono::high_resolution_clock::now();

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

    htsThreadPool tpool = {NULL, 0};
    if (hts_thread_count > 1) {
        if (!(tpool.pool = hts_tpool_init(hts_thread_count))) {
            hts_log_info("Could not initialize thread pool!");
        }

        // set threads
        if (hts_set_thread_pool(infile, &tpool) < 0) {
            LOG_WITH_LEVEL(logging::ERROR) << "Cannot set threads for reading " << input_filepath;
            std::exit(EXIT_FAILURE);
        }

        if (hts_set_thread_pool(outfile, &tpool) < 0) {
            LOG_WITH_LEVEL(logging::ERROR) << "Cannot set threads for writing " << output_filepath;
            RELEASE_TPOOL(tpool.pool);
            std::exit(EXIT_FAILURE);
        }
    }

    // read header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to read header from file!";
        RELEASE_TPOOL(tpool.pool);
        std::exit(EXIT_FAILURE);
    }

    if (sam_hdr_write(outfile, in_samhdr) < 0) {
        LOG_WITH_LEVEL(logging::ERROR) << "Can't write header to bam file: " << output_filepath;
        RELEASE_TPOOL(tpool.pool);
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
                RELEASE_TPOOL(tpool.pool);
                std::exit(EXIT_FAILURE);
            }

            reads_written++;
            current_read_i++;
        }

        id++;
    }

    if (current_read_i != bam_ids.end() && ret_r >= 0) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";
        RELEASE_TPOOL(tpool.pool);
        std::exit(EXIT_FAILURE);
    }

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

    RELEASE_TPOOL(tpool.pool);

    auto end = std::chrono::high_resolution_clock::now();

    LOG_WITH_LEVEL(logging::DEBUG) << "BamApi: " << reads_written << " reads have been written to "
                                   << output_filepath << " on the basis of " << input_filepath;

    std::chrono::duration<double> elapsed = end - start;
    LOG_WITH_LEVEL(logging::DEBUG) << "write_bam took " << elapsed.count() << " seconds";

    return reads_written;
}
