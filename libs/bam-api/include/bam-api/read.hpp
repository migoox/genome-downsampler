#ifndef READ_HPP
#define READ_HPP

#include <htslib/sam.h>

#include <cstdint>

namespace bam_api {

// This id type corresponds to the line in bam file where the read was found
typedef size_t BAMReadId;
// This index type corresponds to the read index in the array(s) in the RAM
typedef size_t ReadIndex;
// This index is an index in the reference genome
typedef size_t Index;
typedef int64_t ReadQuality;

// Read structure definition
struct Read {
    BAMReadId bam_id;
    Index start_ind;
    Index end_ind;
    ReadQuality quality;
    uint32_t seq_length;
    uint32_t mapq;
    int64_t as;
    bool is_first_read;
    bool is_secondary;
    bool is_supplementary;
    bool has_sa;

    Read(BAMReadId id, bam1_t* bamdata);
    Read(BAMReadId id, Index start_ind, Index end_ind, ReadQuality quality, uint32_t seq_length,
         uint8_t mapq, int32_t as, bool is_first, bool is_secondary, bool is_supplementary,
         bool has_sa);
};

}  // namespace bam_api

#endif  // READ_HPP
