#ifndef BAM_API_HPP
#define BAM_API_HPP

#include <filesystem>
#include <vector>

#include "bam_paired_reads.hpp"

namespace bam_api {
class BamApi {
   public:
    static AOSPairedReads read_bam_aos(const std::filesystem::path& filepath,
                                       uint32_t min_seq_length, uint32_t min_mapq);
    static SOAPairedReads read_bam_soa(const std::filesystem::path& filepath,
                                       uint32_t min_seq_length, uint32_t min_mapq);
    static void read_bam(const std::filesystem::path& filepath, PairedReads& paired_reads,
                         uint32_t min_seq_length, uint32_t min_mapq);

    static std::vector<uint32_t> find_cover(const PairedReads& paired_reads);

    static std::vector<uint32_t> find_cover_filtered(const PairedReads& paired_reads,
                                                     const std::vector<ReadIndex>& reads_indices);

    // Returns number of reads written
    static uint32_t write_bam(const std::filesystem::path& input_filepath,
                              const std::filesystem::path& output_filepath,
                              std::vector<ReadIndex>& read_ids);
};

}  // namespace bam_api

#endif
