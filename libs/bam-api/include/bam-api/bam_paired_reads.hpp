#ifndef PAIRED_READS_HPP
#define PAIRED_READS_HPP

#include <cstdint>
#include <map>
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
    std::map<int64_t, int64_t> read_pair_map;
    int64_t ref_genome_length = 0;
};

struct SOAPairedReads : PairedReads {
    std::vector<int64_t> ids;
    std::vector<int64_t> start_inds;
    std::vector<int64_t> end_inds;
    std::vector<uint32_t> qualities;
    std::vector<std::string> qnames;
    std::vector<bool> is_first_reads;
};

struct AOSPairedReads : PairedReads {
    std::vector<Read> reads;

    SOAPairedReads to_soa() {
        SOAPairedReads soa;
        for (auto& read : reads) {
            soa.ids.push_back(read.id);
            soa.start_inds.push_back(read.start_ind);
            soa.end_inds.push_back(read.end_ind);
            soa.qualities.push_back(read.quality);
            soa.qnames.push_back(read.qname);
            soa.is_first_reads.push_back(read.is_first_read);
        }

        return soa;
    }
};

}  // namespace bam_api

#endif
