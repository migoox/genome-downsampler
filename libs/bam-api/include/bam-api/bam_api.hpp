#ifndef BAM_API_HPP
#define BAM_API_HPP

#include <filesystem>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "bam-api/amplicons.hpp"
#include "bam_paired_reads.hpp"

namespace bam_api {

typedef enum PairedReadsType {
    NONE,
    // Structure Of Arrays
    SOA,
    // Array Of Structures
    AOS,
} PairedReadsType;

typedef enum AmpliconBehaviour {
    // Filter out all non-single amplicon inclusive pairs of reads
    FILTER,
    // Grade pairs of reads with its `quality` field to prioritize single amplicon pairs
    GRADE,
} AmpliconBehaviour;

class BamApi {
   public:
    explicit BamApi(const std::filesystem::path& input_filepath)
        : input_filepath_(input_filepath) {}
    explicit BamApi(const AOSPairedReads& paired_reads_)
        : aos_paired_reads_{paired_reads_}, stored_paired_reads_(AOS) {}
    explicit BamApi(const SOAPairedReads& paired_reads_)
        : soa_paired_reads_{paired_reads_}, stored_paired_reads_(SOA) {}

    void add_min_length_filter(uint32_t min_length);
    void add_min_mapq_filter(uint32_t min_mapq);
    void add_amplicon_filter(const std::filesystem::path& bed_filepath,
                             const std::filesystem::path& tsv_filepath = std::filesystem::path());

    void set_amplicon_behaviour(AmpliconBehaviour amplicon_behaviour);

    // if input_filepath provided in constructor
    // reading and filtering occurs on first execution of one of these functions
    const AOSPairedReads& get_paired_reads_aos();
    const SOAPairedReads& get_paired_reads_soa();
    const std::vector<Read>& get_filtered_out_reads() const;

    // Returns number of reads written
    uint32_t write_paired_reads(const std::filesystem::path& output_filepath,
                                const std::vector<BAMReadId>& active_bam_ids) const;
    uint32_t write_bam_api_filtered_out_reads(const std::filesystem::path& output_filepath) const;

    // Testing purposes
    std::vector<uint32_t> find_input_cover();
    std::vector<uint32_t> find_filtered_cover(const std::vector<BAMReadId>& active_bam_ids);

   private:
    SOAPairedReads soa_paired_reads_;
    AOSPairedReads aos_paired_reads_;
    PairedReadsType stored_paired_reads_ = NONE;
    AmpliconSet amplicon_set_;
    AmpliconBehaviour amplicon_behaviour_ = FILTER;
    std::vector<Read> filtered_out_reads_;
    std::filesystem::path input_filepath_;
    uint32_t min_seq_length_ = 0;
    uint32_t min_mapq_ = 0;

    static std::map<std::string, std::pair<Index, Index>> process_bed_file(const std::filesystem::path& filepath);
    static std::vector<std::pair<std::string, std::string>> process_tsv_file(const std::filesystem::path& filepath);
    // Returns number of reads written
    static uint32_t write_bam(const std::filesystem::path& input_filepath,
                              const std::filesystem::path& output_filepath,
                              std::vector<BAMReadId>& bam_ids);
    static void read_bam(const std::filesystem::path& input_filepath, PairedReads& paired_reads,
                         std::vector<Read> filtered_out_reads, uint32_t min_seq_length,
                         uint32_t min_mapq,
                         const std::filesystem::path& bed_amplicon_bounds_filepath,
                         const std::filesystem::path& tsv_amplicon_pairs_filepath);
    const PairedReads& get_paired_reads();
};

}  // namespace bam_api

#endif
