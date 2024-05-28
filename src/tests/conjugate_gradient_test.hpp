#ifndef CONJUGATE_GRADIENT_TESTS_HPP
#define CONJUGATE_GRADIENT_TESTS_HPP
#include <cstdint>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "qmcp-solver/conjugate_gradient_solver.hpp"
#include "test_helpers.hpp"

namespace test {
void small_example_test_conjugate_gradient() {
    // GIVEN
    const uint32_t m = 4;

    auto input = test_helpers::small_soa_reads_example();
    auto input_cover = bam_api::BamApi::find_cover(input);

    qmcp::ConjugateGradientSolver solver;
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
}  // namespace test
#endif
