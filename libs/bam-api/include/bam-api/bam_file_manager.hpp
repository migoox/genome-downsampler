#ifndef BAM_FILE_MANAGER_HPP
#define BAM_FILE_MANAGER_HPP

#include <filesystem>
#include <map>
#include <string>
#include <utility>

#include "bam-api/amplicon_set.hpp"
#include "bam-api/bam_file_config.hpp"
#include "bam-api/read.hpp"
#include "bam-api/region_api.hpp"
namespace bam_api {

constexpr uint32_t kMaxMAPQ = 60;
using Primer = std::tuple<std::string, bam_api::Index, bam_api::Index>;

class BamFileManager {
   public:
    BamFileManager(const std::filesystem::path& input_filepath, const BamFileConfig& config);

    const std::vector<BAMReadId>& get_filtered_out_reads() const;

    // returns number of reads written
    uint32_t write_paired_reads(const std::filesystem::path& output_filepath,
                                std::map<std::string, std::vector<ReadIndex>>& active_ids) const;
    uint32_t write_bam_api_filtered_out_reads(const std::filesystem::path& output_filepath);
    std::map<std::string, RegionApi>& get_regions();

   private:
    std::map<std::string, RegionApi> regions_;
    std::map<std::string, AmpliconSet> amplicon_sets_;
    AmpliconBehaviour amplicon_behaviour_ = AmpliconBehaviour::IGNORE;
    std::vector<BAMReadId> filtered_out_reads_;
    std::filesystem::path input_filepath_;
    uint32_t min_seq_length_ = 0;
    uint32_t amp_overflow_ = 0;
    uint32_t min_mapq_ = 0;
    float min_alignment_ = 0.;
    uint32_t hts_thread_count_ = 1;
    uint32_t min_imported_mapq_ = UINT32_MAX;
    uint32_t max_imported_mapq_ = 0;

    static std::map<std::string, Primer> process_bed_file(const std::filesystem::path& filepath);
    static std::vector<std::pair<std::string, std::string>> process_tsv_file(
        const std::filesystem::path& filepath);

    // Returns number of reads written
    static uint32_t write_bam(const std::filesystem::path& input_filepath,
                              const std::filesystem::path& output_filepath,
                              std::vector<BAMReadId>& bam_ids, uint32_t hts_thread_count);
    void read_bam(const std::filesystem::path& input_filepath);

    void set_amplicon_filter(const std::filesystem::path& bed_filepath,
                             const std::filesystem::path& tsv_filepath = std::filesystem::path());
    void add_to_amplicon_sets(Primer& left_primer, Primer& right_primer);

    void analyse_mapq(const Read& r);

    // Filtering helpers
    bool should_be_filtered_out(std::string& region, const Read& r, bool is_paired) const;
    bool should_be_filtered_out(std::string& region, const Read& r1, const Read& r2) const;
    bool have_min_alignment(const Read& r) const;
    bool have_min_alignment(const Read& r1, const Read& r2) const;
    bool have_min_length(const Read& r) const;
    bool have_min_mapq(const Read& r) const;
    bool are_from_single_amplicon(std::string& region, const Read& r1, const Read& r2) const;
    bool is_from_single_amplicon(std::string& region, const Read& r) const;

    void apply_amplicon_inclusion_grading(PairedReads& paired_reads,
                                          std::vector<bool>& is_in_single_amplicon) const;
};

}  // namespace bam_api

#endif
