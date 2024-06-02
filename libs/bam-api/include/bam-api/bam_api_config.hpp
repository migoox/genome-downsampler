#ifndef BAM_API_CONFIG_HPP
#define BAM_API_CONFIG_HPP

#include <cstdint>
#include <filesystem>

namespace bam_api {

enum class AmpliconBehaviour {
    // Ignore amplicons at all
    IGNORE,
    // Filter out all non-single amplicon inclusive pairs of reads
    FILTER,
    // Grade pairs of reads with its `quality` field to prioritize single amplicon pairs
    GRADE,
};

struct BamApiConfig {
    std::filesystem::path bed_filepath;
    std::filesystem::path tsv_filepath;
    uint32_t min_seq_length = 0;
    uint32_t min_mapq = 0;
    AmpliconBehaviour amplicon_behaviour = AmpliconBehaviour::IGNORE;
};

}  // namespace bam_api

#endif

