#ifndef BAM_API_HPP
#define BAM_API_HPP

#include <filesystem>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "bam-api/amplicon_set.hpp"
#include "bam-api/aos_paired_reads.hpp"
#include "bam-api/bam_api_config.hpp"
#include "bam-api/paired_reads.hpp"
#include "bam-api/read.hpp"
#include "bam-api/soa_paired_reads.hpp"

namespace bam_api {

constexpr uint32_t kMaxMAPQ = 60;

class BamApi {
   public:
    BamApi(const std::filesystem::path& input_filepath, const BamApiConfig& config);
    explicit BamApi(const AOSPairedReads& paired_reads);
    explicit BamApi(const SOAPairedReads& paired_reads);

    // each solver should set it for itself
    void set_amplicon_behaviour(AmpliconBehaviour amplicon_behaviour);

    // if input_filepath provided in constructor
    // reading and filtering occurs on first execution of one of these functions
    const AOSPairedReads& get_paired_reads_aos();
    const SOAPairedReads& get_paired_reads_soa();
    const PairedReads& get_paired_reads() const;
    const std::vector<BAMReadId>& get_filtered_out_reads() const;
    std::vector<ReadIndex> find_pairs(const std::vector<ReadIndex>& ids) const;

    // returns number of reads written
    uint32_t write_paired_reads(const std::filesystem::path& output_filepath,
                                std::vector<ReadIndex>& active_ids) const;
    uint32_t write_bam_api_filtered_out_reads(const std::filesystem::path& output_filepath);

    // testing purposes
    std::vector<uint32_t> find_input_cover();
    std::vector<uint32_t> find_filtered_cover(const std::vector<ReadIndex>& active_ids);

   private:
    SOAPairedReads soa_paired_reads_;
    bool is_soa_loaded_ = false;
    AOSPairedReads aos_paired_reads_;
    bool is_aos_loaded_ = false;
    AmpliconSet amplicon_set_;
    AmpliconBehaviour amplicon_behaviour_ = AmpliconBehaviour::IGNORE;
    std::vector<BAMReadId> filtered_out_reads_;
    std::filesystem::path input_filepath_;
    uint32_t min_seq_length_ = 0;
    uint32_t min_mapq_ = 0;
    uint32_t hts_thread_count_ = 1;
    uint32_t min_imported_mapq_ = UINT32_MAX;
    uint32_t max_imported_mapq_ = 0;

    static std::map<std::string, std::pair<Index, Index>> process_bed_file(
        const std::filesystem::path& filepath);
    static std::vector<std::pair<std::string, std::string>> process_tsv_file(
        const std::filesystem::path& filepath);

    // Returns number of reads written
    static uint32_t write_bam(const std::filesystem::path& input_filepath,
                              const std::filesystem::path& output_filepath,
                              std::vector<BAMReadId>& bam_ids, uint32_t hts_thread_count);
    void read_bam(const std::filesystem::path& input_filepath, PairedReads& paired_reads);

    void set_min_length_filter(uint32_t min_length);
    void set_min_mapq_filter(uint32_t min_mapq);
    void set_amplicon_filter(const std::filesystem::path& bed_filepath,
                             const std::filesystem::path& tsv_filepath = std::filesystem::path());

    void analyse_mapq(const Read& r);

    // Filtering helpers
    bool should_be_filtered_out(const Read& r1, const Read& r2);
    static bool have_min_length(const Read& r1, const Read& r2, uint32_t min_length);
    static bool have_min_mapq(const Read& r1, const Read& r2, uint32_t min_mapq);
    static bool are_from_single_amplicon(const Read& r1, const Read& r2,
                                         const AmpliconSet& amplicon_set);

    void apply_amplicon_inclusion_grading(PairedReads& paired_reads,
                                          std::vector<bool>& is_in_single_amplicon) const;
};

}  // namespace bam_api

#endif
