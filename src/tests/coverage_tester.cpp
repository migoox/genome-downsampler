#include "coverage_tester.hpp"
#include <cassert>
#include <random>

#include "bam-api/read.hpp"
#include "qmcp-solver/solver.hpp"
#include "reads_gen.hpp"

void test::CoverageTester::test(qmcp::Solver& solver) {
   small_example_test(solver);
   random_uniform_dist_test(solver);
   random_low_coverage_on_both_sides_test(solver);
   random_with_hole_test(solver);
   random_zero_coverage_on_both_sides_test(solver);
}

bam_api::AOSPairedReads test::CoverageTester::get_small_aos_example() {
    bam_api::AOSPairedReads result;
    bam_api::ReadIndex id = 0;

    std::vector<bam_api::Read> reads = {
        {id++, 0, 2, 0, 3, true},  {id++, 6, 9, 0, 4, false},  {id++, 2, 4, 0, 3, true},
        {id++, 6, 8, 0, 3, false}, {id++, 1, 3, 0, 3, true},   {id++, 7, 10, 0, 4, false},
        {id++, 3, 6, 0, 4, true},  {id++, 9, 10, 0, 2, false}, {id++, 0, 4, 0, 5, true},
        {id++, 7, 9, 0, 3, false}, {id++, 4, 6, 0, 3, true},   {id++, 9, 10, 0, 2, false},
        {id++, 1, 4, 0, 4, true},  {id++, 6, 8, 0, 3, false},  {id++, 0, 2, 0, 3, true},
        {id++, 4, 6, 0, 3, false},
    };

    result.ref_genome_length = 11;
    for(bam_api::Read& read : reads) {
      result.push_back(read);
    }

    return result;
}

void test::CoverageTester::cap_cover(std::vector<uint32_t>& cover, uint32_t cap) {
    for (uint32_t i = 0; i < cover.size(); ++i) {
        cover[i] = cover[i] > cap ? cap : cover[i];
    }
}

bool test::CoverageTester::is_out_cover_valid(std::vector<uint32_t>& in_cover, const std::vector<uint32_t>& out_cover,
                        uint32_t m) {
    cap_cover(in_cover, m);
    for (uint32_t i = 0; i < out_cover.size(); ++i) {
        if (in_cover[i] > out_cover[i]) {
            return false;
        }
    }

    return true;
}

void test::CoverageTester::small_example_test(qmcp::Solver& solver) {
    // GIVEN
    const uint32_t m = 4;

    auto input = get_small_aos_example();

    bam_api::BamApi bam_api(input);
    auto input_cover = bam_api.find_input_cover();

    // WHEN
    auto output_indices = solver.solve(m, bam_api);
    auto output_cover = bam_api.find_filtered_cover(*output_indices);
    bool valid = is_out_cover_valid(input_cover, output_cover, m);

    // THEN
    assert(valid == true);
}

void test::CoverageTester::random_uniform_dist_test(qmcp::Solver& solver) {
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

    // WHEN
    auto output_indices = solver.solve(m, bam_api);
    auto output_cover = bam_api.find_filtered_cover(*output_indices);
    bool valid = is_out_cover_valid(input_cover, output_cover, m);

    // THEN
    assert(valid == true);
}

void test::CoverageTester::random_with_func_dist_test(const std::function<double(double)>& dist_func, qmcp::Solver& solver) {
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

    // WHEN
    auto output_indices = solver.solve(m, bam_api);
    auto output_cover = bam_api.find_filtered_cover(*output_indices);
    bool valid = is_out_cover_valid(input_cover, output_cover, m);

    // THEN
    assert(valid == true);
}

void test::CoverageTester::random_low_coverage_on_both_sides_test(qmcp::Solver& solver) {
    auto func = [](double x) { return x - x * x; };
    random_with_func_dist_test(func, solver);
}

void test::CoverageTester::random_with_hole_test(qmcp::Solver& solver) {
    auto func = [](double x) {
        if (x > 0.3684 && x < 0.6316) {                                     // NOLINT
            return 1000.0 * (x * x - x + 0.25) * (x * x - x + 0.25) + 0.2;  // NOLINT
        }
        return 0.5;  // NOLINT
    };
    random_with_func_dist_test(func, solver);
}

void test::CoverageTester::random_zero_coverage_on_both_sides_test(qmcp::Solver& solver) {
    auto func = [](double x) {
        return -10.0 * (x - 0.5) * (x - 0.5) + 1.0;  // NOLINT
    };
    random_with_func_dist_test(func, solver);
}

