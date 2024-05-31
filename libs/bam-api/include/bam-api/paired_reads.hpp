#ifndef PAIRED_READS_HPP
#define PAIRED_READS_HPP

#include <htslib/sam.h>

#include <cstdint>
#include <optional>
#include <vector>

namespace bam_api {

// Type definitions
typedef size_t BAMReadId;
typedef size_t ReadIndex;
typedef size_t Index;
typedef uint32_t ReadQuality;

// Forward declarations
struct Read;
struct SOAPairedReads;
struct AOSPairedReads;

// Read structure definition
struct Read {
    BAMReadId bam_id;
    Index start_ind;
    Index end_ind;
    ReadQuality quality;
    uint32_t seq_length;
    bool is_first_read;

    Read(BAMReadId id, bam1_t* bamdata);
    Read(BAMReadId id, Index start_ind, Index end_ind, ReadQuality quality, uint32_t seq_length,
         bool is_first);
};

// PairedReads base structure
struct PairedReads {
    Index ref_genome_length = 0;
    std::vector<std::optional<BAMReadId>> read_pair_map;
    std::vector<std::optional<ReadIndex>> bam_id_to_read_index;

    inline std::optional<Read> get_read_by_bam_id(BAMReadId bam_id) const {
        if (bam_id >= bam_id_to_read_index.size() || !bam_id_to_read_index[bam_id]) {
            return std::nullopt;
        }

        ReadIndex index = bam_id_to_read_index[bam_id].value();
        return get_read_by_index(index);
    }

    virtual void push_back(Read&& read) = 0;
    virtual void push_back(const Read& read) = 0;
    virtual Read get_read_by_index(ReadIndex index) const = 0;
    virtual ReadIndex get_reads_count() const = 0;
    virtual void reserve(size_t size) = 0;
    virtual ~PairedReads() = default;
};

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
    SOAPairedReads& operator=(const AOSPairedReads& aos);
};

// Array of Structures (AoS) implementation
struct AOSPairedReads : PairedReads {
    std::vector<Read> reads;

    inline void push_back(Read&& read) override;
    inline void push_back(const Read& read) override;
    Read get_read_by_index(ReadIndex index) const override;
    ReadIndex get_reads_count() const override;
    inline void reserve(size_t size) override;
    void clear();
    AOSPairedReads& operator=(const SOAPairedReads& soa);
};

}  // namespace bam_api

#endif  // PAIRED_READS_HPP
