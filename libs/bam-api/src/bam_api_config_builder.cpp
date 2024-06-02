#include "bam-api/bam_api_config_builder.hpp"

#include "bam-api/bam_api_config.hpp"

void bam_api::BamApiConfigBuilder::add_min_mapq(uint32_t min_mapq) {
    bam_api_config_.min_mapq = min_mapq;
}
void bam_api::BamApiConfigBuilder::add_min_seq_length(uint32_t min_seq_length) {
    bam_api_config_.min_seq_length = min_seq_length;
}
void bam_api::BamApiConfigBuilder::add_amplicon_filtering(
    AmpliconBehaviour amplicon_behaviour, const std::filesystem::path& bed_filepath,
    const std::filesystem::path& tsv_filepath) {
    if (amplicon_behaviour == AmpliconBehaviour::IGNORE) {
        return;
    }

    bam_api_config_.amplicon_behaviour = amplicon_behaviour;
    bam_api_config_.bed_filepath = bed_filepath;
    bam_api_config_.tsv_filepath = tsv_filepath;
}

void bam_api::BamApiConfigBuilder::add_hts_thread_count(uint32_t hts_thread_count) {
    bam_api_config_.hts_thread_count = hts_thread_count;
}

const bam_api::BamApiConfig& bam_api::BamApiConfigBuilder::get_config() const {
    return bam_api_config_;
}
