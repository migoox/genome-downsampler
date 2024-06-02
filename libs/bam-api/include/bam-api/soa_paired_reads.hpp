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
    std::vector<bool> is_first_reads;

    inline void push_back(Read&& read) override;
    inline void push_back(const Read& read) override;
    Read get_read_by_index(ReadIndex index) const override;
    ReadIndex get_reads_count() const override;
    inline void reserve(size_t size) override;
    void clear();
    SOAPairedReads& from(const AOSPairedReads& aos);
};

}  // namespace bam_api

#endif  // SOA_PAIRED_READS_HPP
