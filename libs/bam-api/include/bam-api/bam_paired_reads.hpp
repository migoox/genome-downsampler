#ifndef PAIRED_READS_HPP
#define PAIRED_READS_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace bam_api {

struct Read {
    int64_t id;
    int64_t start_ind;
    int64_t end_ind;
    uint32_t quality;
    std::string qname;
    bool is_first_read;
};

struct PairedReads {
    std::vector<int64_t> read_pair_map;
    int64_t ref_genome_length;

    virtual void push_back(Read& read) = 0;
    virtual int64_t get_id(int32_t index) = 0;
    virtual ~PairedReads() = default;
};

struct AOSPairedReads : PairedReads {
    std::vector<Read> reads;

    void push_back(Read& read) override { reads.push_back(read); }
    int64_t get_id(int32_t index) override { return reads[index].id; }
};

struct SOAPairedReads : PairedReads {
    std::vector<int64_t> ids;
    std::vector<int64_t> start_inds;
    std::vector<int64_t> end_inds;
    std::vector<uint32_t> qualities;
    std::vector<std::string> qnames;
    std::vector<bool> is_first_reads;

    void push_back(Read& read) override {
        ids.push_back(read.id);
        start_inds.push_back(read.start_ind);
        end_inds.push_back(read.end_ind);
        qualities.push_back(read.quality);
        qnames.push_back(read.qname);
        is_first_reads.push_back(read.is_first_read);
    }
    int64_t get_id(int32_t index) override { return ids[index]; }
};

}  // namespace bam_api

#endif
