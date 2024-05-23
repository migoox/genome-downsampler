#include "../include/sample-data-generator/sample_data_generator.hpp"

namespace sample_data_generator {
void SampleDataGenerator::create_mock_1(bam_api::PairedReads& paired_reads,
                                        unsigned int max_coverage, unsigned int genom_length) {
    // create sample data with covers all nucleotides exactly max_coverage times there will be no
    // pairs

    paired_reads.ref_genome_length = genom_length;

    for (unsigned int nucleotide_index = 0; nucleotide_index < paired_reads.ref_genome_length - 1;
         nucleotide_index++) {
        for (unsigned int i = 0; i < max_coverage; i++) {
            unsigned int total_index = max_coverage * nucleotide_index + i;

            auto read =
                bam_api::Read{total_index, nucleotide_index, nucleotide_index + 1, 60, true};

            paired_reads.push_back(read);
        }
    }
}
}  // namespace sample_data_generator