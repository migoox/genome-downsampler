#ifndef TESTS_COVERAGE_TESTER
#define TESTS_COVERAGE_TESTER

#include <functional>
#include "bam-api/aos_paired_reads.hpp"
#include "qmcp-solver/solver.hpp"
#include "solver_tester.hpp"

namespace test {
class CoverageTester : public SolverTester {
   public:
    void test(qmcp::Solver& solver) override;

   private:
    // tests
    static void small_example_test(qmcp::Solver& solver);
    static void random_uniform_dist_test(qmcp::Solver& solver);
    static void random_low_coverage_on_both_sides_test(qmcp::Solver& solver);
    static void random_with_hole_test(qmcp::Solver& solver);
    static void random_zero_coverage_on_both_sides_test(qmcp::Solver& solver);

    // helpers
    static void random_with_func_dist_test(const std::function<double(double)>& dist_func, qmcp::Solver& solver);
    static bam_api::AOSPairedReads get_small_aos_example();
    static void cap_cover(std::vector<uint32_t>& cover, uint32_t cap);
    static bool is_out_cover_valid(std::vector<uint32_t>& in_cover, const std::vector<uint32_t>& out_cover,
                        uint32_t m);
};
}  // namespace test

#endif
