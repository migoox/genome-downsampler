#include "bam-api/bam_file_config_builder.hpp"

#include "bam-api/bam_file_config.hpp"

void bam_api::BamFileConfigBuilder::add_min_mapq(uint32_t min_mapq) {
    bam_file_config_.min_mapq = min_mapq;
}

void bam_api::BamFileConfigBuilder::add_min_seq_length(uint32_t min_seq_length) {
    bam_file_config_.min_seq_length = min_seq_length;
}

void bam_api::BamFileConfigBuilder::add_amp_overflow(uint32_t amp_overflow) {
    bam_file_config_.amp_overflow = amp_overflow;
}

void bam_api::BamFileConfigBuilder::add_min_alignment(float min_alignment) {
    bam_file_config_.min_alignment = min_alignment;
}

void bam_api::BamFileConfigBuilder::add_amplicon_filtering(
    AmpliconBehaviour amplicon_behaviour, const std::filesystem::path& bed_filepath,
    const std::filesystem::path& tsv_filepath) {
    if (amplicon_behaviour == AmpliconBehaviour::IGNORE) {
        return;
    }

    bam_file_config_.amplicon_behaviour = amplicon_behaviour;
    bam_file_config_.bed_filepath = bed_filepath;
    bam_file_config_.tsv_filepath = tsv_filepath;
}

void bam_api::BamFileConfigBuilder::add_hts_thread_count(uint32_t hts_thread_count) {
    bam_file_config_.hts_thread_count = hts_thread_count;
}

bam_api::BamFileConfig bam_api::BamFileConfigBuilder::build() const { return bam_file_config_; }
