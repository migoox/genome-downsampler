#ifndef BAM_API_HPP
#define BAM_API_HPP

#include <bam-api/bam_data.hpp>
#include <string>

namespace bam_api {

class BamApi {
   public:
    // API calls here
    static void test_func();

    static AOSPairedReads read_bam_aos(std::string filepath);
    static SOAPairedReads read_bam_soa(std::string filepath);

   private:
    // API helper functions here
};

}  // namespace bam_api

#endif
