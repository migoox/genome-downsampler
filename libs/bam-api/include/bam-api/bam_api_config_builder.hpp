#ifndef BAM_API_CONFIG_BUILDER_HPP
#define BAM_API_CONFIG_BUILDER_HPP

#include "bam-api/bam_api_config.hpp"

namespace bam_api {

class BamApiConfigBuilder {
   public:
    void add_min_mapq(uint32_t min_mapq);
    void add_min_seq_length(uint32_t min_seq_length);
    void add_amplicon_filtering(AmpliconBehaviour amplicon_behaviour,
                                const std::filesystem::path& bed_filepath,
                                const std::filesystem::path& tsv_filepath = std::filesystem::path());
    const BamApiConfig& get_config() const;

   private:
    BamApiConfig bam_api_config_;
};

}  // namespace bam_api

#endif
