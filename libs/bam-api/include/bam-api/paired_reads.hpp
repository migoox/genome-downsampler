#ifndef PAIRED_READS_HPP
#define PAIRED_READS_HPP

#include <htslib/sam.h>

#include <optional>
#include <vector>

#include "bam-api/read.hpp"

namespace bam_api {
// PairedReads base structure
struct PairedReads {
    Index ref_genome_length = 0;
    std::vector<std::optional<BAMReadId>> read_pair_map;
    std::vector<std::optional<ReadIndex>> bam_id_to_read_index;

    std::optional<Read> get_read_by_bam_id(BAMReadId bam_id) const;

    virtual void push_back(Read&& read) = 0;
    virtual void push_back(const Read& read) = 0;
    virtual Read get_read_by_index(ReadIndex index) const = 0;
    virtual ReadIndex get_reads_count() const = 0;
    virtual void reserve(size_t size) = 0;
    virtual ~PairedReads() = default;
};

}  // namespace bam_api

#endif  // PAIRED_READS_HPP
