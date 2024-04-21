#ifndef BAM_API_HPP
#define BAM_API_HPP

#include <cstdint>
#include <filesystem>

#include "bam_paired_reads.hpp"

namespace bam_api {

class BamApi {
   public:
    static AOSPairedReads read_bam_aos(const std::filesystem::path& filepath);
    static SOAPairedReads read_bam_soa(const std::filesystem::path& filepath);

    // Returns number of reads written
    static int32_t write_bam_aos(const std::filesystem::path& input_filepath,
                                 const std::filesystem::path& output_filepath,
                                 AOSPairedReads& pairedReads);
    static int32_t write_bam_soa(const std::filesystem::path& input_filepath,
                                 const std::filesystem::path& output_filepath,
                                 SOAPairedReads& pairedReads);

   private:
    static void read_bam(PairedReads& paired_reads,
                         const std::filesystem::path& filepath);
    static int32_t write_bam(const std::filesystem::path& input_filepath,
                             const std::filesystem::path& output_filepath,
                             PairedReads& pairedReads);
};

}  // namespace bam_api

#endif
