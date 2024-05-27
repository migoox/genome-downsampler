#include <htslib/hts.h>
#include <stdio.h>

#include <chrono>
#include <cstdlib>
#include <random>
#include <vector>

#include "app.hpp"
#include "sequential_max_flow_test.hpp"

int main(int argc, char** argv) {
    test::small_example_test();
    test::random_uniform_dist_test();
    test::random_low_coverage_on_both_sides_test();
    test::random_with_hole_test();
    test::random_zero_coverage_on_both_sides_test();
    // test::bam_file_test(
    //     "/home/billyk/Downloads/gpu-programming/data/ESIB_EQA_2023.SARS2.01/reads.bam");
    App app;

    try {
        app.Parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.Exit(e);
    }

    app.Solve();

    return EXIT_SUCCESS;
}
