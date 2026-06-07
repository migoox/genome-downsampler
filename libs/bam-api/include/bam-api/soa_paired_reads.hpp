#ifndef SOA_PAIRED_READS_HPP
#define SOA_PAIRED_READS_HPP

#include <htslib/sam.h>

#include <cstdint>
#include <vector>

#include "bam-api/paired_reads.hpp"
#include "bam-api/read.hpp"

namespace bam_api {

// Forward declaration
struct AOSPairedReads;

// Structure of Arrays (SoA) implementation
struct SOAPairedReads : PairedReads {
    std::vector<BAMReadId> ids;
    std::vector<Index> start_inds;
    std::vector<Index> end_inds;
    std::vector<ReadQuality> qualities;
    std::vector<uint32_t> seq_lengths;
    std::vector<uint32_t> mapqs;
    std::vector<int64_t> as_scores;
    std::vector<bool> is_first_reads;
    std::vector<bool> is_secondary_reads;
    std::vector<bool> is_supplementary_reads;
    std::vector<bool> has_sa_tags;

    inline void push_back(Read&& read) override {
        ids.emplace_back(read.bam_id);
        start_inds.emplace_back(read.start_ind);
        end_inds.emplace_back(read.end_ind);
        qualities.emplace_back(read.quality);
        seq_lengths.emplace_back(read.seq_length);
        mapqs.emplace_back(read.mapq);
        as_scores.emplace_back(read.as);
        is_first_reads.emplace_back(read.is_first_read);
        is_secondary_reads.emplace_back(read.is_secondary);
        is_supplementary_reads.emplace_back(read.is_supplementary);
        has_sa_tags.emplace_back(read.has_sa);
    }

    inline void push_back(const Read& read) override {
        ids.push_back(read.bam_id);
        start_inds.push_back(read.start_ind);
        end_inds.push_back(read.end_ind);
        qualities.push_back(read.quality);
        seq_lengths.push_back(read.seq_length);
        mapqs.push_back(read.mapq);
        as_scores.push_back(read.as);
        is_first_reads.push_back(read.is_first_read);
        is_secondary_reads.push_back(read.is_secondary);
        is_supplementary_reads.push_back(read.is_supplementary);
        has_sa_tags.push_back(read.has_sa);
    }

    inline void reserve(size_t size) override {
        ids.reserve(size);
        start_inds.reserve(size);
        end_inds.reserve(size);
        qualities.reserve(size);
        seq_lengths.reserve(size);
        mapqs.reserve(size);
        as_scores.reserve(size);
        is_first_reads.reserve(size);
        is_secondary_reads.reserve(size);
        is_supplementary_reads.reserve(size);
        has_sa_tags.reserve(size);
    }

    Read get_read_by_index(ReadIndex index) const override;
    ReadQuality get_quality(ReadIndex index) const override;
    void set_quality(ReadIndex index, ReadQuality quality) override;
    ReadIndex get_reads_count() const override;
    void clear();
    SOAPairedReads& from(const AOSPairedReads& aos);
};

}  // namespace bam_api

#endif  // SOA_PAIRED_READS_HPP
