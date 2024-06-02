#include "bam-api/paired_reads.hpp"

std::optional<bam_api::Read> bam_api::PairedReads::get_read_by_bam_id(BAMReadId bam_id) const {
    if (bam_id >= bam_id_to_read_index.size() || !bam_id_to_read_index[bam_id]) {
        return std::nullopt;
    }

    ReadIndex index = bam_id_to_read_index[bam_id].value();
    return get_read_by_index(index);
}
