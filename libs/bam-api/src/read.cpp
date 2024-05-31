#include <htslib/sam.h>
#include "bam-api/paired_reads.hpp"

// Read structure constructors
bam_api::Read::Read(BAMReadId id, bam1_t* bamdata)
    : bam_id(id),
      start_ind(static_cast<Index>(bamdata->core.pos)),
      quality(bamdata->core.qual),
      seq_length(bamdata->core.l_qseq),
      is_first_read(static_cast<bool>(bamdata->core.flag & BAM_FREAD1)) {
    uint64_t rlen =
        bam_cigar2rlen(static_cast<int32_t>(bamdata->core.n_cigar), bam_get_cigar(bamdata));
    end_ind = static_cast<Index>(bamdata->core.pos + rlen - 1);
}

bam_api::Read::Read(BAMReadId id, Index start_ind, Index end_ind, ReadQuality quality, uint32_t seq_length,
           bool is_first)
    : bam_id(id),
      start_ind(start_ind),
      end_ind(end_ind),
      quality(quality),
      seq_length(seq_length),
      is_first_read(is_first) {}
