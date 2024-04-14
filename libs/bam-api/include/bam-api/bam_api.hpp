#ifndef BAM_API_HPP
#define BAM_API_HPP

#include <string>

#include "bam_paired_reads.hpp"

namespace bam_api {

class BamApi {
   public:
    static AOSPairedReads read_bam_aos(const std::string& filepath);
    static SOAPairedReads read_bam_soa(const std::string& filepath);

   private:
    static void read_bam(PairedReads& paired_reads,
                         const std::string& filepath);
};

}  // namespace bam_api

#endif
