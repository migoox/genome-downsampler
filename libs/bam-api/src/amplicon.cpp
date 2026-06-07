#include "bam-api/amplicon.hpp"

bam_api::Amplicon::Amplicon(Index start, Index end, uint32_t overflow)
    : start(start), end(end), overflow(overflow) {}

bool bam_api::Amplicon::includes(const Read& read) const {
    if (start >= overflow)
        return start - overflow <= read.start_ind && read.end_ind < end + overflow;
    return start <= read.start_ind && read.end_ind < end + overflow;
}
