#ifndef BAM_API_AMPLICONS_HPP
#define BAM_API_AMPLICONS_HPP

#include "bam-api/read.hpp"

namespace bam_api {
struct Amplicon {
    Index start;
    Index end;

    Amplicon(Index start, Index end);

    bool includes(const Read& read) const;
};

}  // namespace bam_api

#endif
