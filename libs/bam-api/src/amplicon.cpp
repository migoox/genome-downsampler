#include "bam-api/amplicon.hpp"

bam_api::Amplicon::Amplicon(Index start, Index end) : start(start), end(end) {}

bool bam_api::Amplicon::includes(const Read& read) const {
    return start <= read.start_ind && read.end_ind <= end;
}
