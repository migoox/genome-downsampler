#ifndef PAIRED_READS_HPP
#define PAIRED_READS_HPP

#include <cstdint>
#include <optional>
#include <vector>

namespace bam_api {

typedef size_t ReadIndex;
typedef size_t Index;
typedef uint32_t ReadQuality;

struct Read {
    ReadIndex id;
    Index start_ind;
    Index end_ind;
    ReadQuality quality;
    bool is_first_read;
};

struct PairedReads {
    Index ref_genome_length = 0;
    std::vector<std::optional<ReadIndex>> read_pair_map;

    virtual void push_back(const Read& read) = 0;
    virtual void reserve(size_t size) = 0;
    virtual ~PairedReads() = default;
};

struct SOAPairedReads : PairedReads {
    std::vector<ReadIndex> ids;
    std::vector<Index> start_inds;
    std::vector<Index> end_inds;
    std::vector<ReadQuality> qualities;
    std::vector<bool> is_first_reads;

    void push_back(const Read& read) override {
        ids.push_back(read.id);
        start_inds.push_back(read.start_ind);
        end_inds.push_back(read.end_ind);
        qualities.push_back(read.quality);
        is_first_reads.push_back(read.is_first_read);
    }

    void reserve(size_t size) override {
        ids.reserve(size);
        start_inds.reserve(size);
        end_inds.reserve(size);
        qualities.reserve(size);
        is_first_reads.reserve(size);
    }
};

struct AOSPairedReads : PairedReads {
    std::vector<Read> reads;

    void push_back(const Read& read) override { reads.push_back(read); }
    void reserve(size_t size) override { reads.reserve(size); }
};

}  // namespace bam_api

#endif
