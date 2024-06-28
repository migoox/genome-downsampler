#ifndef TESTS_EFFICIENCY_TESTER
#define TESTS_EFFICIENCY_TESTER

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <vector>

#include "qmcp-solver/solver.hpp"
#include "solver_tester.hpp"

namespace fs = std::filesystem;

namespace test {

// [time ms, m, |S|, |Q|]
typedef std::vector<int64_t> EfficiencyTestResult;

class EfficiencyTester : public SolverTester {
   public:
    using EfficiencyTestFunction = EfficiencyTestResult (*)(qmcp::Solver&);

    void test(qmcp::Solver& solver, fs::path& outputs_dir_path_) override;

   private:
    // tests
    static EfficiencyTestResult random_uniform_dist_test_small(qmcp::Solver& solver);
    static EfficiencyTestResult random_low_coverage_on_both_sides_test_small(qmcp::Solver& solver);
    static EfficiencyTestResult random_with_hole_test_small(qmcp::Solver& solver);
    static EfficiencyTestResult random_zero_coverage_on_both_sides_test_small(qmcp::Solver& solver);
    static EfficiencyTestResult random_uniform_dist_test_medium(qmcp::Solver& solver);
    static EfficiencyTestResult random_low_coverage_on_both_sides_test_medium(qmcp::Solver& solver);
    static EfficiencyTestResult random_with_hole_test_medium(qmcp::Solver& solver);
    static EfficiencyTestResult random_zero_coverage_on_both_sides_test_medium(qmcp::Solver& solver);
    static EfficiencyTestResult random_uniform_dist_test_large(qmcp::Solver& solver);
    static EfficiencyTestResult random_low_coverage_on_both_sides_test_large(qmcp::Solver& solver);
    static EfficiencyTestResult random_with_hole_test_large(qmcp::Solver& solver);
    static EfficiencyTestResult random_zero_coverage_on_both_sides_test_large(qmcp::Solver& solver);

    // helpers
    static void run_test_and_write_output(qmcp::Solver& solver, std::ofstream& results_file,
                                          EfficiencyTestFunction test_func, const std::string& func_name);
    static void write_results(EfficiencyTestResult& result, fs::path& output_filepath);
    static EfficiencyTestResult random_with_func_dist_test(
        const std::function<double(double)>& dist_func, qmcp::Solver& solver, uint32_t pairs_count,
        uint32_t genome_length, uint32_t m);
    static EfficiencyTestResult random_with_uniform_dist_test(qmcp::Solver& solver,
                                                                     uint32_t pairs_count,
                                                                     uint32_t genome_length,
                                                                     uint32_t m);
};
}  // namespace test

#endif
