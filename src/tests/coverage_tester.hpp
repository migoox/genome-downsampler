#ifndef TESTS_COVERAGE_TESTER
#define TESTS_COVERAGE_TESTER

#include <cstdint>
#include <filesystem>
#include <functional>
#include <utility>
#include <vector>

#include "bam-api/aos_paired_reads.hpp"
#include "qmcp-solver/solver.hpp"
#include "solver_tester.hpp"

namespace fs = std::filesystem;

namespace test {

typedef std::pair<std::vector<uint32_t>, std::vector<uint32_t>> CoverageTestResult;

class CoverageTester : public SolverTester {
   public:
    using CoverageTestFunction = CoverageTestResult (*)(qmcp::Solver&);

    void test(qmcp::Solver& solver, fs::path& outputs_dir_path_) override;

   private:
    // tests
    static CoverageTestResult small_example_test(qmcp::Solver& solver);
    static CoverageTestResult random_uniform_dist_test(qmcp::Solver& solver);
    static CoverageTestResult random_low_coverage_on_both_sides_test(qmcp::Solver& solver);
    static CoverageTestResult random_with_hole_test(qmcp::Solver& solver);
    static CoverageTestResult random_zero_coverage_on_both_sides_test(qmcp::Solver& solver);

    // helpers
    static void run_test_and_write_output(qmcp::Solver& solver, fs::path& outputs_dir_path_,
                                          const std::string& output_filename,
                                          CoverageTestFunction test_func);
    static void write_covers(CoverageTestResult& result, fs::path& output_filepath);
    static CoverageTestResult random_with_func_dist_test(
        const std::function<double(double)>& dist_func, qmcp::Solver& solver);
    static bam_api::AOSPairedReads get_small_aos_example();
    static void cap_cover(std::vector<uint32_t>& cover, uint32_t cap);
    static bool is_out_cover_valid(std::vector<uint32_t>& in_cover,
                                   const std::vector<uint32_t>& out_cover, uint32_t m);
};
}  // namespace test

#endif
