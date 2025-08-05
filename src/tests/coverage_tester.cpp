#include "coverage_tester.hpp"

#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>

#include "bam-api/read.hpp"
#include "logging/log.hpp"
#include "qmcp-solver/solver.hpp"
#include "reads_gen.hpp"
#include "scoped_timer.hpp"

#define RUN_TEST_FUNCTION(SOLVER, OUTPUT, FUNCTION)                                            \
    LOG_WITH_LEVEL(logging::INFO) << "Running " << #FUNCTION << "...";                         \
    {                                                                                          \
        ScopedTimer timer;                                                                     \
        run_test_and_write_output(SOLVER, OUTPUT, std::string(#FUNCTION) + ".cov", &FUNCTION); \
    }                                                                                          \
    LOG_WITH_LEVEL(logging::INFO) << "PASSED!";

namespace fs = std::filesystem;

namespace test {

void CoverageTester::test(const std::unique_ptr<qmcp::Solver>& solver,
                          fs::path& outputs_dir_path_) {
    if (outputs_dir_path_.empty()) {
        small_example_test(solver);
        random_uniform_dist_test(solver);
        random_low_coverage_on_both_sides_test(solver);
        random_with_hole_test(solver);
        random_zero_coverage_on_both_sides_test(solver);
        return;
    }

    RUN_TEST_FUNCTION(solver, outputs_dir_path_, small_example_test);
    RUN_TEST_FUNCTION(solver, outputs_dir_path_, random_uniform_dist_test);
    RUN_TEST_FUNCTION(solver, outputs_dir_path_, random_low_coverage_on_both_sides_test);
    RUN_TEST_FUNCTION(solver, outputs_dir_path_, random_with_hole_test);
    RUN_TEST_FUNCTION(solver, outputs_dir_path_, random_zero_coverage_on_both_sides_test);
}

void CoverageTester::run_test_and_write_output(const std::unique_ptr<qmcp::Solver>& solver,
                                               fs::path& outputs_dir_path_,
                                               const std::string& output_filename,
                                               CoverageTestFunction test_func) {
    CoverageTestResult result;
    fs::path output_path = outputs_dir_path_ / output_filename;
    result = (*test_func)(solver);
    write_covers(result, output_path);
}

void CoverageTester::write_covers(CoverageTestResult& result, fs::path& output_filepath) {
    const auto& input_cov = result.first;
    const auto& output_cov = result.second;

    assert(input_cov.size() == output_cov.size());

    std::ofstream file(output_filepath);
    if (file.is_open()) {
        for (size_t i = 0; i < input_cov.size(); ++i) {
            file << i << "\t" << input_cov[i] << "\t" << output_cov[i] << '\n';
        }
        file.close();
    } else {
        LOG_WITH_LEVEL(logging::ERROR) << "Cannot open or create file: " << output_filepath;
        exit(EXIT_FAILURE);
    }
}

bam_api::AOSPairedReads CoverageTester::get_small_aos_example() {
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

    assert(id == 16);

    result.ref_genome_length = 11;
    for (auto& read : reads) {
        result.push_back(read);
    }

    return result;
}

void CoverageTester::cap_cover(std::vector<uint32_t>& cover, uint32_t cap) {
    for (auto& c : cover) {
        c = std::min(c, cap);
    }
}

bool CoverageTester::is_out_cover_valid(std::vector<uint32_t>& in_cover,
                                        const std::vector<uint32_t>& out_cover, uint32_t m) {
    std::vector<uint32_t> capped_input_cover(in_cover);
    cap_cover(capped_input_cover, m);
    return std::equal(capped_input_cover.begin(), capped_input_cover.end(), out_cover.begin(),
                      std::less_equal<>());
}

CoverageTestResult CoverageTester::small_example_test(const std::unique_ptr<qmcp::Solver>& solver) {
    const uint32_t m = 4;
    auto input = get_small_aos_example();
    bam_api::BamApi bam_api(input);
    auto input_cover = bam_api.find_input_cover();
    auto output_indices = solver->solve(m, bam_api);
    auto output_cover = bam_api.find_filtered_cover(*output_indices);
    assert(is_out_cover_valid(input_cover, output_cover, m));
    return {input_cover, output_cover};
}

CoverageTestResult CoverageTester::random_uniform_dist_test(
    const std::unique_ptr<qmcp::Solver>& solver) {
    const uint32_t seed = 12345;
    const uint32_t pairs_count = 1'000'000;
    const uint32_t genome_length = 30'000;
    const uint32_t read_length = 150;
    const uint32_t m = 1000;

    std::mt19937 mt(seed);
    auto input = reads_gen::rand_reads_uniform(mt, pairs_count, genome_length, read_length);

    bam_api::BamApi bam_api(input);
    auto input_cover = bam_api.find_input_cover();
    auto output_indices = solver->solve(m, bam_api);
    auto output_cover = bam_api.find_filtered_cover(*output_indices);
    assert(is_out_cover_valid(input_cover, output_cover, m));
    return {input_cover, output_cover};
}

CoverageTestResult CoverageTester::random_with_func_dist_test(
    const std::function<double(double)>& dist_func, const std::unique_ptr<qmcp::Solver>& solver) {
    const uint32_t seed = 12345;
    const uint32_t pairs_count = 1'000'000;
    const uint32_t genome_length = 30'000;
    const uint32_t read_length = 150;
    const uint32_t m = 8000;

    std::mt19937 mt(seed);
    auto input = reads_gen::rand_reads(mt, pairs_count, genome_length, read_length, dist_func);

    bam_api::BamApi bam_api(input);
    auto input_cover = bam_api.find_input_cover();
    auto output_indices = solver->solve(m, bam_api);
    auto output_cover = bam_api.find_filtered_cover(*output_indices);
    assert(is_out_cover_valid(input_cover, output_cover, m));
    return {input_cover, output_cover};
}

CoverageTestResult CoverageTester::random_low_coverage_on_both_sides_test(
    const std::unique_ptr<qmcp::Solver>& solver) {
    auto func = [](double x) { return x - x * x; };
    return random_with_func_dist_test(func, solver);
}

CoverageTestResult CoverageTester::random_with_hole_test(
    const std::unique_ptr<qmcp::Solver>& solver) {
    auto func = [](double x) {
        if (x > 0.3684 && x < 0.6316) {
            return 1000.0 * (x * x - x + 0.25) * (x * x - x + 0.25) + 0.2;
        }
        return 0.5;
    };
    return random_with_func_dist_test(func, solver);
}

CoverageTestResult CoverageTester::random_zero_coverage_on_both_sides_test(
    const std::unique_ptr<qmcp::Solver>& solver) {
    auto func = [](double x) { return -10.0 * (x - 0.5) * (x - 0.5) + 1.0; };
    return random_with_func_dist_test(func, solver);
}

}  // namespace test
