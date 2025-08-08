#include <htslib/sam.h>

#include "bam-api/paired_reads.hpp"

// Read structure constructors
bam_api::Read::Read(BAMReadId id, bam1_t* bamdata)
    : bam_id(id),
      start_ind(static_cast<Index>(bamdata->core.pos)),
      seq_length(bamdata->core.l_qseq),
      mapq(bamdata->core.qual),
      is_first_read(static_cast<bool>(bamdata->core.flag & BAM_FREAD1)),
      is_secondary(static_cast<bool>(bamdata->core.flag & BAM_FSECONDARY)),
      is_supplementary(static_cast<bool>(bamdata->core.flag & BAM_FSUPPLEMENTARY)),
      has_sa(bam_aux_get(bamdata, "SA") != NULL) {
    uint8_t* as_tag = bam_aux_get(bamdata, "AS");
    as = (as_tag != NULL) ? bam_aux2i(as_tag) : 0;
    quality = as + mapq;
    uint64_t rlen =
        bam_cigar2rlen(static_cast<int32_t>(bamdata->core.n_cigar), bam_get_cigar(bamdata));
    end_ind = static_cast<Index>(bamdata->core.pos + rlen - 1);
}

bam_api::Read::Read(BAMReadId id, Index start_ind, Index end_ind, ReadQuality quality,
                    uint32_t seq_length, uint8_t mapq, int32_t as, bool is_first, bool is_secondary,
                    bool is_supplementary, bool has_sa)
    : bam_id(id),
      start_ind(start_ind),
      end_ind(end_ind),
      quality(quality),
      seq_length(seq_length),
      mapq(mapq),
      as(as),
      is_first_read(is_first),
      is_secondary(is_secondary),
      is_supplementary(is_supplementary),
      has_sa(has_sa) {}
