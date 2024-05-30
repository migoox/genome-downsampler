// TODO(billyk): introduce testing framework
#ifndef SEQUENTIAL_MAX_FLOW_TESTS_HPP
#define SEQUENTIAL_MAX_FLOW_TESTS_HPP
#include <cstdint>
#include <filesystem>
#include <random>
#include <vector>

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

    for (bam_api::ReadIndex i; i < result.get_reads_count(); ++i) {
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

void small_example_test() {
    // GIVEN
    const uint32_t m = 4;

    auto input = test_helpers::small_aos_reads_example();

    bam_api::BamApi bam_api(input);
    auto input_cover = bam_api.find_input_cover();

    qmcp::SequentialMaxFlowSolver solver(bam_api);

    // WHEN
    auto output_indices = solver.solve(m);
    auto output_cover = bam_api.find_filtered_cover(output_indices);
#ifdef TESTS_VERBOSE_DATA
    test_helpers::print_vectors(input_cover, output_cover);
#endif

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

    bam_api::BamApi bam_api(input);
    auto input_cover = bam_api.find_input_cover();

    qmcp::SequentialMaxFlowSolver solver(bam_api);

    // WHEN
    auto output_indices = solver.solve(m);
    auto output_cover = bam_api.find_filtered_cover(output_indices);
#ifdef TESTS_VERBOSE_DATA
    test_helpers::print_vectors(input_cover, output_cover);
#endif

    bool valid = test_helpers::is_out_cover_valid(input_cover, output_cover, m);

    // THEN
    assert(valid == true);
}

void random_with_func_dist_test(const std::function<double(double)>& dist_func) {
    // GIVEN
    const uint32_t seed = 12345;
    const uint32_t pairs_count = 1'000'000;
    const uint32_t genome_length = 30'000;
    const uint32_t read_length = 150;
    const uint32_t m = 8000;

    std::mt19937 mt(seed);
    auto input = reads_gen::rand_reads(mt, pairs_count, genome_length, read_length, dist_func);

    bam_api::BamApi bam_api(input);
    auto input_cover = bam_api.find_input_cover();

    qmcp::SequentialMaxFlowSolver solver(bam_api);

    // WHEN
    auto output_indices = solver.solve(m);
    auto output_cover = bam_api.find_filtered_cover(output_indices);
#ifdef TESTS_VERBOSE_DATA
    test_helpers::print_vectors(input_cover, output_cover);
#endif

    bool valid = test_helpers::is_out_cover_valid(input_cover, output_cover, m);

    // THEN
    assert(valid == true);
}

void random_low_coverage_on_both_sides_test() {
    auto func = [](double x) { return x - x * x; };
    random_with_func_dist_test(func);
}

void random_with_hole_test() {
    auto func = [](double x) {
        if (x > 0.3684 && x < 0.6316) {                                     // NOLINT
            return 1000.0 * (x * x - x + 0.25) * (x * x - x + 0.25) + 0.2;  // NOLINT
        }
        return 0.5;  // NOLINT
    };
    random_with_func_dist_test(func);
}

void random_zero_coverage_on_both_sides_test() {
    auto func = [](double x) {
        return -10.0 * (x - 0.5) * (x - 0.5) + 1.0;  // NOLINT
    };
    random_with_func_dist_test(func);
}

void bam_file_test(const std::filesystem::path& path) {
    // GIVEN
    const uint32_t m = 1000;
    bam_api::BamApiConfig config{
        .min_seq_length = 100,
        .min_mapq = 30,
    };
    bam_api::BamApi bam_api(path, config);

    qmcp::SequentialMaxFlowSolver solver(bam_api);

    // WHEN
    auto output_indices = solver.solve(m);

    auto input_cover = bam_api.find_input_cover();
    auto output_cover = bam_api.find_filtered_cover(output_indices);
#ifdef TESTS_VERBOSE_DATA
    test_helpers::print_vectors(input_cover, output_cover);
#endif

    bool valid = test_helpers::is_out_cover_valid(input_cover, output_cover, m);

    // THEN
    assert(valid == true);
}
}  // namespace test

#endif
