#include "bam-api/bam_file_manager.hpp"

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
#include <tuple>
#include <vector>

#include "bam-api/amplicon_set.hpp"
#include "bam-api/paired_reads.hpp"
#include "bam-api/read.hpp"
#include "bam-api/region_api.hpp"
#include "logging/log.hpp"

#define RELEASE_TPOOL(X)                                  \
    {                                                     \
        hts_tpool* ptr = static_cast<hts_tpool*>(X); \
        if (ptr) {                                        \
            hts_tpool_destroy(ptr);                       \
        }                                                 \
    }

bam_api::BamFileManager::BamFileManager(const std::filesystem::path& input_filepath,
                                        const BamFileConfig& config)
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

    read_bam(input_filepath_);
}

void bam_api::BamFileManager::set_amplicon_filter(const std::filesystem::path& bed_filepath,
                                                  const std::filesystem::path& tsv_filepath) {
    auto start = std::chrono::high_resolution_clock::now();

    std::map<std::string, Primer> primer_map;
    primer_map = process_bed_file(bed_filepath);

    if (!tsv_filepath.empty()) {
        std::vector<std::pair<std::string, std::string>> pairs;
        pairs = process_tsv_file(tsv_filepath);

        for (const auto& [left, right] : pairs) {
            add_to_amplicon_sets(primer_map[left], primer_map[right]);
        }
    } else {
        auto it = primer_map.begin();

        while (it != primer_map.end()) {
            auto& left_primer = it->second;
            it++;
            auto& right_primer = it->second;

            add_to_amplicon_sets(left_primer, right_primer);

            it++;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    LOG_WITH_LEVEL(logging::DEBUG)
        << "BamFileManager: set_amplicon_filter took " << elapsed.count() << " seconds";
}

void bam_api::BamFileManager::add_to_amplicon_sets(Primer& left_primer, Primer& right_primer) {
    auto& [left_reg, left_beg, left_end] = left_primer;
    auto& [right_reg, right_beg, right_end] = right_primer;

    if (left_beg > right_beg) {
        std::swap(left_reg, right_reg);
        std::swap(left_beg, right_beg);
        std::swap(left_end, right_end);
    }

    if (left_reg != right_reg) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "Wrong pairings of amplicons in .tsv file, regions (chromosomes do not match)";
        std::exit(EXIT_FAILURE);
    }

    if (amplicon_sets_.find(left_reg) == amplicon_sets_.end()) {
        AmpliconSet new_reg_set;
        new_reg_set.amplicons.push_back({left_beg, right_end, amp_overflow_});
        amplicon_sets_[left_reg] = new_reg_set;
    }

    amplicon_sets_[left_reg].amplicons.push_back({left_beg, right_end, amp_overflow_});
}

std::map<std::string, std::tuple<std::string, bam_api::Index, bam_api::Index>>
bam_api::BamFileManager::process_bed_file(const std::filesystem::path& filepath) {
    LOG_WITH_LEVEL(logging::INFO) << "Reading " << filepath.filename() << " file...";

    std::map<std::string, std::tuple<std::string, bam_api::Index, bam_api::Index>> ret;

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
            ret.emplace(name, std::make_tuple(chrom, start_idx, end_idx));
        } else {
            LOG_WITH_LEVEL(logging::ERROR) << "Invalid BED line: " << line;
        }
    }

    file.close();

    LOG_WITH_LEVEL(logging::DEBUG) << "BamFileManager: " << ret.size() << " primers have been read";

    return ret;
}

std::vector<std::pair<std::string, std::string>> bam_api::BamFileManager::process_tsv_file(
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

    LOG_WITH_LEVEL(logging::DEBUG)
        << "BamFileManager: " << ret.size() << " pairs of primers have been read";

    return ret;
}

const std::vector<bam_api::BAMReadId>& bam_api::BamFileManager::get_filtered_out_reads() const {
    return filtered_out_reads_;
}

std::map<std::string, bam_api::RegionApi>& bam_api::BamFileManager::get_regions() {
    return regions_;
}

bool bam_api::BamFileManager::should_be_filtered_out(std::string& region, const Read& r,
                                                     bool is_paired) const {
    bool ret = false;

    ret = ret || !have_min_mapq(r);
    ret = ret || !have_min_length(r);
    if (!is_paired) ret = ret || !have_min_alignment(r);
    ret = ret || r.has_sa;

    if (amplicon_behaviour_ == AmpliconBehaviour::FILTER) {
        ret = ret || !is_from_single_amplicon(region, r);
    }

    return ret;
}

bool bam_api::BamFileManager::should_be_filtered_out(std::string& region, const Read& r1,
                                                     const Read& r2) const {
    bool ret = should_be_filtered_out(region, r1, true) || should_be_filtered_out(region, r2, true);

    if (amplicon_behaviour_ == AmpliconBehaviour::FILTER) {
        ret = ret || !are_from_single_amplicon(region, r1, r2);
    }

    return ret;
}

bool bam_api::BamFileManager::have_min_alignment(const Read& r) const {
    return (static_cast<float>(r.as) / static_cast<float>(r.seq_length)) >= min_alignment_;
}

bool bam_api::BamFileManager::have_min_alignment(const Read& r1, const Read& r2) const {
    return (static_cast<float>(r1.as) / static_cast<float>(r1.seq_length + r2.seq_length)) >=
           min_alignment_;
}

bool bam_api::BamFileManager::have_min_length(const Read& r) const {
    return r.seq_length >= min_seq_length_;
}

bool bam_api::BamFileManager::have_min_mapq(const Read& r) const { return r.mapq >= min_mapq_; }

bool bam_api::BamFileManager::are_from_single_amplicon(std::string& region, const Read& r1,
                                                       const Read& r2) const {
    return amplicon_sets_.at(region).member_includes_both(r1, r2);
}

bool bam_api::BamFileManager::is_from_single_amplicon(std::string& region, const Read& r) const {
    return amplicon_sets_.at(region).any_includes(r);
}

void bam_api::BamFileManager::apply_amplicon_inclusion_grading(
    bam_api::PairedReads& paired_reads, std::vector<bool>& is_in_single_amplicon) const {
    LOG_WITH_LEVEL(logging::DEBUG) << "BamFileManager: grading min_mapq: " << min_imported_mapq_
                                   << ", max_mapq: " << max_imported_mapq_;
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

void bam_api::BamFileManager::analyse_mapq(const Read& r) {
    if (r.quality < min_imported_mapq_) {
        min_imported_mapq_ = r.mapq;
    }

    if (r.quality > max_imported_mapq_) {
        max_imported_mapq_ = r.mapq;
    }
}

void bam_api::BamFileManager::read_bam(const std::filesystem::path& input_filepath) {
    LOG_WITH_LEVEL(logging::INFO) << "Reading " << input_filepath.filename() << " input file...";
    LOG_WITH_LEVEL(logging::DEBUG) << "BamFileManager: reading " << input_filepath;

    auto start = std::chrono::high_resolution_clock::now();

    sam_hdr_t* in_samhdr = nullptr;
    samFile* infile = nullptr;
    bam1_t* bamdata = nullptr;
    hts_idx_t* idx = nullptr;
    htsThreadPool tpool = {nullptr, 0};
    int ret_r = 0;

    auto cleanup = [&]() {
        if (idx) hts_idx_destroy(idx);
        if (in_samhdr) sam_hdr_destroy(in_samhdr);
        if (infile) sam_close(infile);
        if (bamdata) bam_destroy1(bamdata);
        RELEASE_TPOOL(tpool.pool);
    };

    bamdata = bam_init1();
    if (!bamdata) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to allocate data memory!";
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    infile = sam_open(input_filepath.c_str(), "r");
    if (!infile) {
        LOG_WITH_LEVEL(logging::ERROR) << "Could not open " << input_filepath;
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    if (hts_thread_count_ > 1) {
        if (!(tpool.pool = hts_tpool_init(hts_thread_count_))) {
            hts_log_info("Could not initialize thread pool!");
        }
        if (hts_set_thread_pool(infile, &tpool) < 0) {
            LOG_WITH_LEVEL(logging::ERROR) << "Cannot set threads for reading " << input_filepath;
            cleanup();
            std::exit(EXIT_FAILURE);
        }
    }

    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to read header from file!";
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    idx = sam_index_load(infile, input_filepath.c_str());
    if (!idx) {
        if (sam_index_build(input_filepath.c_str(), 0) != 0) {
            LOG_WITH_LEVEL(logging::ERROR) << "Failed to build index!";
            cleanup();
            std::exit(EXIT_FAILURE);
        }
        idx = sam_index_load(infile, input_filepath.c_str());
        if (!idx) {
            LOG_WITH_LEVEL(logging::ERROR) << "Failed to load index file!";
            cleanup();
            std::exit(EXIT_FAILURE);
        }
    }

    BAMReadId id = 0;
    std::vector<bool> is_accepted;
    for (int tid = 0; tid < in_samhdr->n_targets; ++tid) {
        uint32_t tlen = in_samhdr->target_len[tid];
        std::string rname = in_samhdr->target_name[tid];
        std::vector<bool> is_in_single_amplicon;
        bool is_paired = true;
        SOAPairedReads paired_reads;
        paired_reads.ref_genome_length = tlen;

        uint64_t mapped_seq_c, unmapped_seq_c;
        hts_idx_get_stat(idx, tid, &mapped_seq_c, &unmapped_seq_c);
        paired_reads.reserve(mapped_seq_c);
        is_in_single_amplicon.reserve(mapped_seq_c);

        std::map<std::string, Read> read_map;
        std::string current_qname;

        hts_itr_t* iter = sam_itr_queryi(idx, tid, 0, tlen);
        while ((ret_r = sam_itr_next(infile, iter, bamdata)) >= 0) {
            Read current_read(id, bamdata);
            current_qname = bam_get_qname(bamdata);

            is_accepted.push_back(false);
            auto it = read_map.find(current_qname);

            if (it != read_map.end()) {
                Read& r1 = it->second;
                Read& r2 = current_read;

                if (should_be_filtered_out(rname, r1, r2)) {
                    id++;
                    continue;
                }

                if (amplicon_behaviour_ == AmpliconBehaviour::GRADE) {
                    analyse_mapq(r1);
                    analyse_mapq(r2);
                    bool same_amplicon = are_from_single_amplicon(rname, r1, r2);
                    is_in_single_amplicon.push_back(same_amplicon);
                    is_in_single_amplicon.push_back(same_amplicon);
                }

                if (r2.is_first_read) std::swap(r1, r2);

                paired_reads.push_back(r1);
                paired_reads.push_back(r2);

                is_accepted[r1.bam_id] = true;
                is_accepted[r2.bam_id] = true;
            } else {
                read_map.insert({current_qname, current_read});
            }
            id++;
        }
        hts_itr_destroy(iter);

        if (ret_r < -1) {
            LOG_WITH_LEVEL(logging::ERROR)
                << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";
            cleanup();
            std::exit(EXIT_FAILURE);
        }

        if (paired_reads.get_reads_count() == 0) {
            is_paired = false;
            for (auto& [qname, read] : read_map) {
                if (should_be_filtered_out(rname, read, false)) continue;

                if (amplicon_behaviour_ == AmpliconBehaviour::GRADE) {
                    analyse_mapq(read);
                    is_in_single_amplicon.push_back(is_from_single_amplicon(rname, read));
                }

                paired_reads.push_back(read);
                is_accepted[read.bam_id] = true;
            }
        }

        if (amplicon_behaviour_ == AmpliconBehaviour::GRADE)
            apply_amplicon_inclusion_grading(paired_reads, is_in_single_amplicon);

        regions_.emplace(rname, RegionApi(paired_reads, is_paired));
    }

    assert(is_accepted.size() == id);

    for (BAMReadId bam_id = 0; bam_id < is_accepted.size(); ++bam_id) {
        if (!is_accepted[bam_id]) filtered_out_reads_.push_back(bam_id);
    }

    assert((paired_reads.get_reads_count() + filtered_out_reads_.size()) == is_accepted.size());

    cleanup();

    auto end = std::chrono::high_resolution_clock::now();
    LOG_WITH_LEVEL(logging::DEBUG) << "BamFileManager: " << id << " reads have been read";
    LOG_WITH_LEVEL(logging::DEBUG)
        << "BamFileManager: " << regions_.size() << " regions have been imported";
    LOG_WITH_LEVEL(logging::DEBUG) << "BamFileManager: " << filtered_out_reads_.size()
                                   << " reads have been filtered out during preprocessing";

    LOG_WITH_LEVEL(logging::DEBUG)
        << "BamFileManager: read_bam took " << std::chrono::duration<double>(end - start).count()
        << " seconds";
}

uint32_t bam_api::BamFileManager::write_paired_reads(
    const std::filesystem::path& output_filepath,
    std::map<std::string, std::vector<ReadIndex>>& active_ids) const {
    LOG_WITH_LEVEL(logging::INFO) << "Writing solution of size " << active_ids.size() << " reads "
                                  << output_filepath.filename() << "...";

    std::vector<BAMReadId> active_bam_ids;
    for (const auto& [rname, region_api] : regions_) {
        const PairedReads& paired_reads = region_api.get_paired_reads();

        for (const auto& id : active_ids[rname]) {
            active_bam_ids.push_back(paired_reads.get_read_by_index(id).bam_id);
        }
    }

    return write_bam(input_filepath_, output_filepath, active_bam_ids, hts_thread_count_);
}

uint32_t bam_api::BamFileManager::write_bam_api_filtered_out_reads(
    const std::filesystem::path& output_filepath) {
    LOG_WITH_LEVEL(logging::INFO) << "Writing " << filtered_out_reads_.size()
                                  << " preprocessing filtered out reads to "
                                  << output_filepath.filename() << "...";
    return write_bam(input_filepath_, output_filepath, filtered_out_reads_, hts_thread_count_);
}

uint32_t bam_api::BamFileManager::write_bam(const std::filesystem::path& input_filepath,
                                            const std::filesystem::path& output_filepath,
                                            std::vector<BAMReadId>& bam_ids,
                                            uint32_t hts_thread_count) {
    LOG_WITH_LEVEL(logging::DEBUG) << "BamFileManager: writing " << bam_ids.size() << " reads to "
                                   << output_filepath << " on the basis of " << input_filepath;

    auto start = std::chrono::high_resolution_clock::now();

    sam_hdr_t* in_samhdr = nullptr;
    samFile* infile = nullptr;
    samFile* outfile = nullptr;
    hts_idx_t* idx = nullptr;
    bam1_t* bamdata = nullptr;
    htsThreadPool tpool = {nullptr, 0};
    uint32_t reads_written = 0;
    int ret_r = 0;

    auto cleanup = [&]() {
        if (in_samhdr) sam_hdr_destroy(in_samhdr);
        if (idx) hts_idx_destroy(idx);
        if (infile) sam_close(infile);
        if (outfile) sam_close(outfile);
        if (bamdata) bam_destroy1(bamdata);
        RELEASE_TPOOL(tpool.pool);
    };

    bamdata = bam_init1();
    if (!bamdata) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to allocate data memory!";
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    infile = sam_open(input_filepath.c_str(), "r");
    if (!infile) {
        LOG_WITH_LEVEL(logging::ERROR) << "Could not open " << input_filepath;
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    std::string open_mode = output_filepath.extension() == ".bam" ? "wb" : "w";
    outfile = sam_open(output_filepath.c_str(), open_mode.c_str());
    if (!outfile) {
        LOG_WITH_LEVEL(logging::ERROR) << "Could not open " << output_filepath;
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    if (hts_thread_count > 1) {
        if (!(tpool.pool = hts_tpool_init(hts_thread_count))) {
            hts_log_info("Could not initialize thread pool!");
        }
        if (hts_set_thread_pool(infile, &tpool) < 0) {
            LOG_WITH_LEVEL(logging::ERROR) << "Cannot set threads for reading " << input_filepath;
            cleanup();
            std::exit(EXIT_FAILURE);
        }
        if (hts_set_thread_pool(outfile, &tpool) < 0) {
            LOG_WITH_LEVEL(logging::ERROR) << "Cannot set threads for writing " << output_filepath;
            cleanup();
            std::exit(EXIT_FAILURE);
        }
    }

    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to read header from file!";
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    idx = sam_index_load(infile, input_filepath.c_str());
    if (!idx) {
        LOG_WITH_LEVEL(logging::ERROR) << "Failed to load index file!";
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    if (sam_hdr_write(outfile, in_samhdr) < 0) {
        LOG_WITH_LEVEL(logging::ERROR) << "Can't write header to bam file: " << output_filepath;
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    BAMReadId id = 0;
    std::sort(bam_ids.begin(), bam_ids.end());
    auto current_read_i = bam_ids.begin();

    for (int tid = 0; tid < in_samhdr->n_targets; ++tid) {
        uint32_t tlen = in_samhdr->target_len[tid];
        hts_itr_t* iter = sam_itr_queryi(idx, tid, 0, tlen);

        while ((ret_r = sam_itr_next(infile, iter, bamdata)) >= 0 &&
               current_read_i != bam_ids.end()) {
            if (id == *current_read_i) {
                if (sam_write1(outfile, in_samhdr, bamdata) < 0) {
                    LOG_WITH_LEVEL(logging::ERROR)
                        << "Can't write line to bam file: " << output_filepath;
                    hts_itr_destroy(iter);
                    cleanup();
                    std::exit(EXIT_FAILURE);
                }
                ++reads_written;
                ++current_read_i;
            }
            ++id;
        }
        hts_itr_destroy(iter);
    }

    if (current_read_i != bam_ids.end() && ret_r < -1) {
        LOG_WITH_LEVEL(logging::ERROR)
            << "Failed to read bam file (sam_read1 error code:" << ret_r << ")";
        cleanup();
        std::exit(EXIT_FAILURE);
    }

    cleanup();

    auto end = std::chrono::high_resolution_clock::now();
    LOG_WITH_LEVEL(logging::DEBUG)
        << "BamFileManager: " << reads_written << " reads have been written to " << output_filepath
        << " on the basis of " << input_filepath;
    LOG_WITH_LEVEL(logging::DEBUG)
        << "BamFileManager: write_bam took " << std::chrono::duration<double>(end - start).count()
        << " seconds";

    return reads_written;
}
