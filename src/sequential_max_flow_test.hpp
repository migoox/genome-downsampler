// TODO(billyk): introduce testing framework
#ifndef SEQUENTIAL_MAX_FLOW_TESTS_HPP
#define SEQUENTIAL_MAX_FLOW_TESTS_HPP
#include <cstdint>
#include <random>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/sequential_max_flow_solver.hpp"
#include "reads_gen.hpp"

namespace test {
namespace test_helpers {

bam_api::AOSPairedReads small_aos_reads_example() {
    bam_api::AOSPairedReads result;
    bam_api::ReadIndex id = 0;

    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 0, 2, 0, true));
    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 6, 9, 0, false));

    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 2, 4, 0, true));
    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 6, 8, 0, false));

    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 1, 3, 0, true));
    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 7, 10, 0, false));

    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 3, 6, 0, true));
    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 9, 10, 0, false));

    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 0, 4, 0, true));
    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 7, 9, 0, false));

    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 4, 6, 0, true));
    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 9, 10, 0, false));

    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 1, 4, 0, true));
    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 6, 8, 0, false));

    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 0, 2, 0, true));
    // NOLINTNEXTLINE
    result.reads.emplace_back(bam_api::Read(id++, 4, 6, 0, false));

    // NOLINTNEXTLINE
    result.ref_genome_length = 11;

    return result;
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

void small_example_test() {
    // GIVEN
    const uint32_t m = 4;

    auto input = test_helpers::small_aos_reads_example();
    auto input_cover = bam_api::BamApi::find_cover(input);

    qmcp::SequentialMaxFlowSolver solver;
    solver.set_reads(input);

    // WHEN
    solver.solve(m);
    auto output_indices = solver.get_output();
    auto output_cover = bam_api::BamApi::find_cover_filtered(input, output_indices);

    bool valid = test_helpers::is_out_cover_valid(input_cover, output_cover, m);

    // THEN
    assert(valid == true);
}

void random_uniform_dist_test() {
    // GIVEN
    const uint32_t seed = 12345;
    const uint32_t pairs_count = 1'000'000;
    const uint32_t genome_length = 30'000;
    const uint32_t read_length = 150;
    const uint32_t m = 1000;

    std::mt19937 mt(seed);
    auto input = reads_gen::rand_reads_uniform(mt, pairs_count, genome_length, read_length);
    auto input_cover = bam_api::BamApi::find_cover(input);

    qmcp::SequentialMaxFlowSolver solver;
    solver.set_reads(input);

    // WHEN
    solver.solve(m);
    auto output_indices = solver.get_output();
    auto output_cover = bam_api::BamApi::find_cover_filtered(input, output_indices);

    bool valid = test_helpers::is_out_cover_valid(input_cover, output_cover, m);

    // THEN
    assert(valid == true);
}

void data_with_hole_test() {
    // TODO(billyk)
}

}  // namespace test

#endif