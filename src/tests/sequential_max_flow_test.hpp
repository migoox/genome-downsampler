// TODO(billyk): introduce testing framework
#ifndef SEQUENTIAL_MAX_FLOW_TESTS_HPP
#define SEQUENTIAL_MAX_FLOW_TESTS_HPP
#include <cstdint>
#include <filesystem>
#include <random>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "test_helpers.hpp"
#include "qmcp-solver/sequential_max_flow_solver.hpp"
#include "reads_gen.hpp"

namespace test {
void small_example_test() {
    // GIVEN
    const uint32_t m = 4;

    auto input = test_helpers::small_aos_reads_example();
    auto input_cover = bam_api::BamApi::find_cover(input);

    qmcp::SequentialMaxFlowSolver solver;
    solver.find_pairs(false);
    solver.set_reads(input);

    // WHEN
    solver.solve(m);
    auto output_indices = solver.get_output();
    auto output_cover = bam_api::BamApi::find_cover_filtered(input, output_indices);
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
    auto input_cover = bam_api::BamApi::find_cover(input);

    qmcp::SequentialMaxFlowSolver solver;
    solver.find_pairs(false);
    solver.set_reads(input);

    // WHEN
    solver.solve(m);
    auto output_indices = solver.get_output();
    auto output_cover = bam_api::BamApi::find_cover_filtered(input, output_indices);
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
    auto input_cover = bam_api::BamApi::find_cover(input);

    qmcp::SequentialMaxFlowSolver solver;
    solver.find_pairs(false);
    solver.set_reads(input);

    // WHEN
    solver.solve(m);
    auto output_indices = solver.get_output();
    auto output_cover = bam_api::BamApi::find_cover_filtered(input, output_indices);
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
    auto input = bam_api::BamApi::read_bam_aos(path, 100, 30);

    auto input_cover = bam_api::BamApi::find_cover(input);

    qmcp::SequentialMaxFlowSolver solver;
    solver.find_pairs(true);
    solver.set_reads(input);

    // WHEN
    solver.solve(m);
    auto output_indices = solver.get_output();
    auto output_cover = bam_api::BamApi::find_cover_filtered(input, output_indices);
#ifdef TESTS_VERBOSE_DATA
    test_helpers::print_vectors(input_cover, output_cover);
#endif
    bool valid = test_helpers::is_out_cover_valid(input_cover, output_cover, m);

    // THEN
    assert(valid == true);
}
}  // namespace test

#endif
