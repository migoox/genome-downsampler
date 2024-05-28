#include <cstdlib>

#include "app.hpp"
#include "tests/sequential_max_flow_test.hpp"
#include "tests/conjugate_gradient_test.hpp"

int main(int argc, char** argv) {
    // test::small_example_test();
    // test::random_uniform_dist_test();
    // test::random_low_coverage_on_both_sides_test();
    // test::random_with_hole_test();
    // test::random_zero_coverage_on_both_sides_test();
    // test::bam_file_test(
    //   "/home/mytkom/Documents/gpu-programming/data/ESIB_EQA_2023.SARS2.01/reads.bam");
    test::small_example_test_conjugate_gradient();

    // App app;

    // try {
    //     app.Parse(argc, argv);
    // } catch (const CLI::ParseError& e) {
    //     return app.Exit(e);
    // }

    // app.Solve();

    return EXIT_SUCCESS;
}
