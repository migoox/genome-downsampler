#include "efficiency_tester.hpp"

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>

#include "logging/log.hpp"
#include "qmcp-solver/solver.hpp"
#include "reads_gen.hpp"
#include "scoped_timer.hpp"

#define RUN_EFF_TEST_FUNCTION(SOLVER, OUTPUT, FUNCTION)                               \
    LOG_WITH_LEVEL(logging::INFO) << "Running " << #FUNCTION << "...";                \
    {                                                                                 \
        ScopedTimer timer;                                                            \
        run_test_and_write_output(SOLVER, OUTPUT, &FUNCTION, std::string(#FUNCTION)); \
    }                                                                                 \
    LOG_WITH_LEVEL(logging::INFO) << "PASSED!";

namespace fs = std::filesystem;

namespace test {

void EfficiencyTester::test(qmcp::Solver& solver, fs::path& outputs_dir_path_) {
    if (outputs_dir_path_.empty()) return;

    fs::path output_path = outputs_dir_path_ / "results.csv";
    std::ofstream results_file;

    if (std::filesystem::exists(output_path)) {
        results_file.open(output_path, std::ios_base::app);
    } else {
        results_file.open(output_path);
        results_file << "name;time_ms;m;genome_length;reads_count" << std::endl;
    }

    if (!results_file.is_open()) {
        LOG_WITH_LEVEL(logging::ERROR) << "Cannot open or create file: " << output_path;
        exit(EXIT_FAILURE);
    }

    RUN_EFF_TEST_FUNCTION(solver, results_file, random_uniform_dist_test_small);
    RUN_EFF_TEST_FUNCTION(solver, results_file, random_low_coverage_on_both_sides_test_small);
    RUN_EFF_TEST_FUNCTION(solver, results_file, random_with_hole_test_small);
    RUN_EFF_TEST_FUNCTION(solver, results_file, random_zero_coverage_on_both_sides_test_small);

    RUN_EFF_TEST_FUNCTION(solver, results_file, random_uniform_dist_test_medium);
    RUN_EFF_TEST_FUNCTION(solver, results_file, random_low_coverage_on_both_sides_test_medium);
    RUN_EFF_TEST_FUNCTION(solver, results_file, random_with_hole_test_medium);
    RUN_EFF_TEST_FUNCTION(solver, results_file, random_zero_coverage_on_both_sides_test_medium);

    RUN_EFF_TEST_FUNCTION(solver, results_file, random_uniform_dist_test_large);
    RUN_EFF_TEST_FUNCTION(solver, results_file, random_low_coverage_on_both_sides_test_large);
    RUN_EFF_TEST_FUNCTION(solver, results_file, random_with_hole_test_large);
    RUN_EFF_TEST_FUNCTION(solver, results_file, random_zero_coverage_on_both_sides_test_large);

    results_file.close();
}

void EfficiencyTester::run_test_and_write_output(qmcp::Solver& solver, std::ofstream& results_file,
                                                 EfficiencyTestFunction test_func,
                                                 const std::string& func_name) {
    EfficiencyTestResult result;
    result = (*test_func)(solver);

    assert(result.size() == 4);

    results_file << func_name << ";";
    results_file << result[0] << ";";
    results_file << result[1] << ";";
    results_file << result[2] << ";";
    results_file << result[3] << std::endl;
}

EfficiencyTestResult EfficiencyTester::random_with_func_dist_test(
    const std::function<double(double)>& dist_func, qmcp::Solver& solver, uint32_t pairs_count,
    uint32_t genome_length, uint32_t m) {
    const uint32_t seed = 12345;
    const uint32_t read_length = 150;
    int64_t solve_time_ms = 0;

    std::mt19937 mt(seed);
    auto input = reads_gen::rand_reads(mt, pairs_count, genome_length, read_length, dist_func);

    bam_api::BamApi bam_api(input);
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto output_indices = solver.solve(m, bam_api);
        auto end = std::chrono::high_resolution_clock::now();
        solve_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }

    return {solve_time_ms, m, genome_length, pairs_count * 2};
}

EfficiencyTestResult EfficiencyTester::random_with_uniform_dist_test(qmcp::Solver& solver,
                                                                     uint32_t pairs_count,
                                                                     uint32_t genome_length,
                                                                     uint32_t m) {
    const uint32_t seed = 12345;
    const uint32_t read_length = 150;
    int64_t solve_time_ms = 0;

    std::mt19937 mt(seed);
    auto input = reads_gen::rand_reads_uniform(mt, pairs_count, genome_length, read_length);

    bam_api::BamApi bam_api(input);
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto output_indices = solver.solve(m, bam_api);
        auto end = std::chrono::high_resolution_clock::now();
        solve_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }

    return {solve_time_ms, m, genome_length, pairs_count * 2};
}

EfficiencyTestResult EfficiencyTester::random_uniform_dist_test_small(qmcp::Solver& solver) {
    const uint32_t pairs_count = 50000;
    const uint32_t genome_length = 3000;
    const uint32_t m = 400;

    return random_with_uniform_dist_test(solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_low_coverage_on_both_sides_test_small(
    qmcp::Solver& solver) {
    const uint32_t pairs_count = 50000;
    const uint32_t genome_length = 3000;
    const uint32_t m = 400;

    auto func = [](double x) { return x - x * x; };
    return random_with_func_dist_test(func, solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_with_hole_test_small(qmcp::Solver& solver) {
    const uint32_t pairs_count = 50000;
    const uint32_t genome_length = 3000;
    const uint32_t m = 400;

    auto func = [](double x) {
        if (x > 0.3684 && x < 0.6316) {
            return 1000.0 * (x * x - x + 0.25) * (x * x - x + 0.25) + 0.2;
        }
        return 0.5;
    };

    return random_with_func_dist_test(func, solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_zero_coverage_on_both_sides_test_small(
    qmcp::Solver& solver) {
    const uint32_t pairs_count = 50000;
    const uint32_t genome_length = 3000;
    const uint32_t m = 400;

    auto func = [](double x) { return -10.0 * (x - 0.5) * (x - 0.5) + 1.0; };
    return random_with_func_dist_test(func, solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_uniform_dist_test_medium(qmcp::Solver& solver) {
    const uint32_t pairs_count = 250000;
    const uint32_t genome_length = 15000;
    const uint32_t m = 2000;

    return random_with_uniform_dist_test(solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_low_coverage_on_both_sides_test_medium(
    qmcp::Solver& solver) {
    const uint32_t pairs_count = 250000;
    const uint32_t genome_length = 15000;
    const uint32_t m = 2000;

    auto func = [](double x) { return x - x * x; };
    return random_with_func_dist_test(func, solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_with_hole_test_medium(qmcp::Solver& solver) {
    const uint32_t pairs_count = 250000;
    const uint32_t genome_length = 15000;
    const uint32_t m = 2000;

    auto func = [](double x) {
        if (x > 0.3684 && x < 0.6316) {
            return 1000.0 * (x * x - x + 0.25) * (x * x - x + 0.25) + 0.2;
        }
        return 0.5;
    };

    return random_with_func_dist_test(func, solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_zero_coverage_on_both_sides_test_medium(
    qmcp::Solver& solver) {
    const uint32_t pairs_count = 250000;
    const uint32_t genome_length = 15000;
    const uint32_t m = 2000;

    auto func = [](double x) { return -10.0 * (x - 0.5) * (x - 0.5) + 1.0; };
    return random_with_func_dist_test(func, solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_uniform_dist_test_large(qmcp::Solver& solver) {
    const uint32_t pairs_count = 500000;
    const uint32_t genome_length = 30000;
    const uint32_t m = 4000;

    return random_with_uniform_dist_test(solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_low_coverage_on_both_sides_test_large(
    qmcp::Solver& solver) {
    const uint32_t pairs_count = 500000;
    const uint32_t genome_length = 30000;
    const uint32_t m = 4000;

    auto func = [](double x) { return x - x * x; };
    return random_with_func_dist_test(func, solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_with_hole_test_large(qmcp::Solver& solver) {
    const uint32_t pairs_count = 500000;
    const uint32_t genome_length = 30000;
    const uint32_t m = 4000;

    auto func = [](double x) {
        if (x > 0.3684 && x < 0.6316) {
            return 1000.0 * (x * x - x + 0.25) * (x * x - x + 0.25) + 0.2;
        }
        return 0.5;
    };

    return random_with_func_dist_test(func, solver, pairs_count, genome_length, m);
}

EfficiencyTestResult EfficiencyTester::random_zero_coverage_on_both_sides_test_large(
    qmcp::Solver& solver) {
    const uint32_t pairs_count = 500000;
    const uint32_t genome_length = 30000;
    const uint32_t m = 4000;

    auto func = [](double x) { return -10.0 * (x - 0.5) * (x - 0.5) + 1.0; };
    return random_with_func_dist_test(func, solver, pairs_count, genome_length, m);
}

}  // namespace test
