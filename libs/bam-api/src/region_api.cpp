#include "bam-api/region_api.hpp"

#include <chrono>
#include <unordered_set>
#include <vector>

#include "bam-api/paired_reads.hpp"
#include "bam-api/read.hpp"
#include "logging/log.hpp"

bam_api::RegionApi::RegionApi(const AOSPairedReads& paired_reads, bool is_paired)
    : is_paired_{is_paired}, aos_paired_reads_{paired_reads}, is_aos_loaded_(true) {}
bam_api::RegionApi::RegionApi(const SOAPairedReads& paired_reads, bool is_paired)
    : is_paired_{is_paired}, soa_paired_reads_{paired_reads}, is_soa_loaded_(true) {}

const bam_api::AOSPairedReads& bam_api::RegionApi::get_paired_reads_aos() {
    if (is_aos_loaded_) {
        return aos_paired_reads_;
    }

    if (is_soa_loaded_) {
        aos_paired_reads_.from(soa_paired_reads_);
        is_aos_loaded_ = true;
        return aos_paired_reads_;
    }

    LOG_WITH_LEVEL(logging::ERROR) << "There is no paired reads set loaded. BamApi panicked!";
    exit(EXIT_FAILURE);
}

const bam_api::SOAPairedReads& bam_api::RegionApi::get_paired_reads_soa() {
    if (is_soa_loaded_) {
        return soa_paired_reads_;
    }

    if (is_aos_loaded_) {
        soa_paired_reads_.from(aos_paired_reads_);
        is_soa_loaded_ = true;
        return soa_paired_reads_;
    }

    LOG_WITH_LEVEL(logging::ERROR) << "There is no paired reads set loaded. BamApi panicked!";
    exit(EXIT_FAILURE);
}

const bam_api::PairedReads& bam_api::RegionApi::get_paired_reads() const {
    if (is_soa_loaded_) {
        return soa_paired_reads_;
    }

    return aos_paired_reads_;
}

std::vector<uint32_t> bam_api::RegionApi::find_input_cover() const {
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

std::vector<bam_api::ReadIndex> bam_api::RegionApi::find_pairs(
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

std::unordered_set<bam_api::ReadIndex> bam_api::RegionApi::add_pseudo_reads() {
    std::vector<uint32_t> coverage = find_input_cover();
    std::vector<std::pair<uint32_t, uint32_t>> zero_ranges;
    bool in_zero = false;
    size_t start = 0;

    for (size_t i = 0; i < coverage.size(); ++i) {
        if (!in_zero && coverage[i] == 0) {
            in_zero = true;
            start = i > 0 ? i - 1 : i;
        } else if (in_zero && coverage[i] != 0) {
            in_zero = false;
            zero_ranges.emplace_back(start, i);
        }
    }

    if (in_zero) zero_ranges.emplace_back(start, coverage.size() - 1);

    std::unordered_set<ReadIndex> ids;
    if (is_aos_loaded_) {
        for (auto& [start_id, end_id] : zero_ranges) {
            ids.emplace(aos_paired_reads_.get_reads_count());
            aos_paired_reads_.push_back({coverage.size(), start_id, end_id, 0, end_id - start_id, 0,
                                         0, false, true, false, false});
        }
        is_soa_loaded_ = false;
    }

    if (is_soa_loaded_) {
        for (auto& [start_id, end_id] : zero_ranges) {
            ids.emplace(soa_paired_reads_.get_reads_count());
            soa_paired_reads_.push_back({coverage.size(), start_id, end_id, 0, end_id - start_id, 0,
                                         0, false, true, false, false});
        }
        is_aos_loaded_ = false;
    }

    return ids;
}

std::vector<uint32_t> bam_api::RegionApi::find_filtered_cover(
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
