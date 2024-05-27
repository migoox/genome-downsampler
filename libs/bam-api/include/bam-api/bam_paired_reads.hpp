#ifndef PAIRED_READS_HPP
#define PAIRED_READS_HPP

#include <cstdint>
#include <cstdlib>
#include <optional>
#include <vector>

#include "logging/log.hpp"

namespace bam_api {

// This id type corresponds to the line in bam file where the read was found
typedef size_t BAMReadId;
// This index type corresponds to the read index in the array(s) in the RAM
typedef size_t ReadIndex;
// This index is an index in the reference genome
typedef size_t Index;
typedef uint32_t ReadQuality;

struct Read {
    BAMReadId bam_id;
    Index start_ind;
    Index end_ind;
    ReadQuality quality;
    bool is_first_read;

    Read(BAMReadId id, Index start_ind, Index end_ind, ReadQuality quality, bool is_first)
        : bam_id(id),
          start_ind(start_ind),
          end_ind(end_ind),
          quality(quality),
          is_first_read(is_first) {}
};

struct PairedReads {
    Index ref_genome_length = 0;

    // maps BAMReadId -> BAMReadId of paired read
    // IT IS NOT ReadIndex
    std::vector<std::optional<BAMReadId>> read_pair_map;

    // maps BAMReadId -> ReadIndex
    std::vector<std::optional<ReadIndex>> bam_id_to_read_index;

    inline Read get_read_by_bam_id(BAMReadId bam_id) const {
        if (!bam_id_to_read_index[bam_id]) {
            LOG_WITH_LEVEL(logging::kError) << "No read of BAM id " << bam_id << "!";
            exit(EXIT_FAILURE);
        }

        ReadIndex index = bam_id_to_read_index[bam_id].value();
        return get_read_by_index(index);
    }

    virtual void push_back(const Read& read) = 0;
    virtual Read get_read_by_index(ReadIndex index) const = 0;
    virtual ReadIndex get_reads_count() const = 0;
    virtual void reserve(size_t size) = 0;
    virtual ~PairedReads() = default;
};

struct SOAPairedReads : PairedReads {
    std::vector<BAMReadId> ids;
    std::vector<Index> start_inds;
    std::vector<Index> end_inds;
    std::vector<ReadQuality> qualities;
    std::vector<bool> is_first_reads;

    inline void push_back(const Read& read) override {
        ids.push_back(read.bam_id);
        start_inds.push_back(read.start_ind);
        end_inds.push_back(read.end_ind);
        qualities.push_back(read.quality);
        is_first_reads.push_back(read.is_first_read);
    }

    inline Read get_read_by_index(ReadIndex index) const override {
        return Read(ids[index], start_inds[index], end_inds[index], qualities[index],
                    is_first_reads[index]);
    }

    inline ReadIndex get_reads_count() const override { return ids.size(); }

    inline void reserve(size_t size) override {
        ids.reserve(size);
        start_inds.reserve(size);
        end_inds.reserve(size);
        qualities.reserve(size);
        is_first_reads.reserve(size);
    }
};

struct AOSPairedReads : PairedReads {
    std::vector<Read> reads;

    inline void push_back(const Read& read) override { reads.push_back(read); }
    inline Read get_read_by_index(ReadIndex index) const override { return reads[index]; }
    inline ReadIndex get_reads_count() const override { return reads.size(); }
    inline void reserve(size_t size) override { reads.reserve(size); }
};

}  // namespace bam_api

#endif
