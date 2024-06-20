#ifndef AOS_PAIRED_READS_HPP
#define AOS_PAIRED_READS_HPP

#include <htslib/sam.h>

#include <vector>

#include "bam-api/paired_reads.hpp"
#include "bam-api/read.hpp"

namespace bam_api {
// Forward declaration
struct SOAPairedReads;

// Array of Structures (AoS) implementation
struct AOSPairedReads : PairedReads {
    std::vector<Read> reads;

    inline void push_back(Read&& read) override { reads.emplace_back(read); }
    inline void push_back(const Read& read) override { reads.push_back(read); }
    inline void reserve(size_t size) override { reads.reserve(size); }
    ReadQuality get_quality(ReadIndex index) const override;
    void set_quality(ReadIndex index, ReadQuality quality) override;
    Read get_read_by_index(ReadIndex index) const override;
    ReadIndex get_reads_count() const override;
    void clear();
    AOSPairedReads& from(const SOAPairedReads& soa);
};

}  // namespace bam_api

#endif  // AOS_PAIRED_READS_HPP
