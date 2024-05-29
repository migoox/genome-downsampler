#ifndef BAM_API_AMPLICONS_HPP
#define BAM_API_AMPLICONS_HPP

#include <algorithm>
#include <vector>

#include "bam_paired_reads.hpp"

namespace bam_api {
struct Amplicon {
    Index start;
    Index end;

    Amplicon(Index s, Index e) : start(s), end(e) {}

    bool includes(const Read& read) const { return start <= read.start_ind && read.end_ind <= end; }
};

struct AmpliconSet {
    std::vector<Amplicon> amplicons;

    bool member_includes_both(const Read& r1, const Read& r2) const {
        return std::any_of(amplicons.cbegin(), amplicons.cend(),
                           [&r1, &r2](const Amplicon& amplicon) {
                               return amplicon.includes(r1) && amplicon.includes(r2);
                           });
    }
};
}  // namespace bam_api

#endif
