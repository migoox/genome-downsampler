#include "bam-api/amplicon_set.hpp"

#include <algorithm>

bool bam_api::AmpliconSet::member_includes_both(const Read& r1, const Read& r2) const {
    return std::any_of(amplicons.cbegin(), amplicons.cend(), [&r1, &r2](const Amplicon& amplicon) {
        return amplicon.includes(r1) && amplicon.includes(r2);
    });
}
