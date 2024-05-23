#ifndef SAMPLE_DATA_GENERATOR_HPP
#define SAMPLE_DATA_GENERATOR_HPP

#include "../../../bam-api/include/bam-api/bam_paired_reads.hpp"

namespace sample_data_generator {
class SampleDataGenerator {
   public:
    bam_api::AOSPairedReads static SOA_to_AOS_converter(bam_api::SOAPairedReads paired_reads);
    bam_api::SOAPairedReads static AOS_to_SOA_converter(bam_api::AOSPairedReads paired_reads);

    void static create_mock_1(bam_api::PairedReads& paired_reads, unsigned int max_coverage,
                              unsigned int genom_length);
};

}  // namespace sample_data_generator

#endif