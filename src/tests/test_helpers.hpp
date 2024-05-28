#ifndef TEST_HELPERS_HPP
#define TEST_HELPERS_HPP
#include <cstdint>
#include <vector>

#include "bam-api/bam_paired_reads.hpp"

namespace test {
namespace test_helpers {

std::vector<bam_api::Read> small_reads_example_data = {
        {0, 0, 2, 0, true},  {1, 6, 9, 0, false},  {2, 2, 4, 0, true},  {3, 6, 8, 0, false},
        {4, 1, 3, 0, true},  {5, 7, 10, 0, false}, {6, 3, 6, 0, true},  {7, 9, 10, 0, false},
        {8, 0, 4, 0, true},  {9, 7, 9, 0, false},  {10, 4, 6, 0, true}, {11, 9, 10, 0, false},
        {12, 1, 4, 0, true}, {13, 6, 8, 0, false}, {14, 0, 2, 0, true}, {15, 4, 6, 0, false},
    };

bam_api::AOSPairedReads small_aos_reads_example() {
    bam_api::AOSPairedReads result;
    
    for (auto& read : small_reads_example_data) {
        result.push_back(read);
    }

    for (bam_api::ReadIndex i = 0; i < result.get_reads_count(); ++i) {
        result.bam_id_to_read_index.push_back(i);
    }

    // NOLINTNEXTLINE
    result.ref_genome_length = 11;

    return result;
}

bam_api::SOAPairedReads small_soa_reads_example() {
    bam_api::SOAPairedReads result;

    for (auto& read : small_reads_example_data) {
      result.push_back(read);
    }

    for (bam_api::ReadIndex i = 0; i < result.get_reads_count(); ++i) {
        result.bam_id_to_read_index.push_back(i);
    }

    // NOLINTNEXTLINE
    result.ref_genome_length = 11;

    return result;
}

void print_vectors(const std::vector<uint32_t>& vec1, const std::vector<uint32_t>& vec2) {
    for (int i = 0; i < vec1.size(); ++i) {
        std::cout << i << "\t" << vec1[i] << "\t" << vec2[i] << std::endl;
    }
}

void cap_cover(std::vector<uint32_t>& cover, uint32_t cap) {
    for (uint32_t i = 0; i < cover.size(); ++i) {
        cover[i] = cover[i] > cap ? cap : cover[i];
    }
}

bool is_out_cover_valid(std::vector<uint32_t>& in_cover, const std::vector<uint32_t>& out_cover,
                        uint32_t m) {
    test_helpers::cap_cover(in_cover, m);
    for (uint32_t i = 0; i < out_cover.size(); ++i) {
        if (in_cover[i] > out_cover[i]) {
            return false;
        }
    }

    return true;
}
}  // namespace test_helpers
}  // namespace test

#endif
