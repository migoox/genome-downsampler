#ifndef BAM_API_HPP
#define BAM_API_HPP

#include <filesystem>

#include "bam_paired_reads.hpp"

namespace bam_api {
class BamApi {
   public:
    static AOSPairedReads read_bam_aos(const std::filesystem::path& filepath);
    static SOAPairedReads read_bam_soa(const std::filesystem::path& filepath);
    static void read_bam(const std::filesystem::path& filepath,
                         PairedReads& paired_reads);
    // Returns number of reads written
    static void write_sam(const std::filesystem::path& input_filepath,
                              const std::filesystem::path& output_filepath,
                              std::vector<ReadIndex>& read_ids,
                              bool use_bam = false);
};

}  // namespace bam_api

#endif
