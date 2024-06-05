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
typedef uint32_t ReadQuality;

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

}  // namespace bam_api

#endif  // READ_HPP
