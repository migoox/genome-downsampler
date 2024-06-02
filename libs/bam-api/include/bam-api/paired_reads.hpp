#ifndef PAIRED_READS_HPP
#define PAIRED_READS_HPP

#include <htslib/sam.h>

#include "bam-api/read.hpp"

namespace bam_api {
// PairedReads base structure
struct PairedReads {
    Index ref_genome_length = 0;

    virtual void push_back(Read&& read) = 0;
    virtual void push_back(const Read& read) = 0;
    virtual Read get_read_by_index(ReadIndex index) const = 0;
    virtual ReadQuality get_quality(ReadIndex index) const = 0;
    virtual void set_quality(ReadIndex index, ReadQuality quality) = 0;
    virtual ReadIndex get_reads_count() const = 0;
    virtual void reserve(size_t size) = 0;
    virtual ~PairedReads() = default;
};

}  // namespace bam_api

#endif  // PAIRED_READS_HPP
