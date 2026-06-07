#ifndef AMPLICON_SET_HPP
#define AMPLICON_SET_HPP

#include <vector>

#include "bam-api/amplicon.hpp"

namespace bam_api {

struct AmpliconSet {
    std::vector<Amplicon> amplicons;

    bool member_includes_both(const Read& r1, const Read& r2) const;
    bool any_includes(const Read& r) const;
};

}  // namespace bam_api

#endif
