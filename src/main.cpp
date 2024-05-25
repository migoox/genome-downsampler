#include <htslib/hts.h>
#include <stdio.h>

#include <chrono>
#include <cstdlib>

#include "app.hpp"
#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/cuda_max_flow_solver.hpp"
#include "qmcp-solver/qmcp-solver.hpp"
#include "qmcp-solver/sequential_cost_scaling_network_solver.hpp"
#include "qmcp-solver/sequential_max_flow_solver.hpp"

bam_api::Read create_read(bool is_first, bam_api::Index start, bam_api::Index end,
                          bam_api::ReadIndex id) {
    return bam_api::Read{id, start, end, 0, is_first};
}

#define ADD_PAIR(where, start1, end1, start2, end2)               \
    (where).push_back(create_read(true, (start1), (end1), id++)); \
    (where).push_back(create_read(false, (start2), (end2), id++));

bam_api::AOSPairedReads aos_reads_1() {
    bam_api::AOSPairedReads result;
    bam_api::ReadIndex id = 0;
    ADD_PAIR(result, 0, 2, 6, 9)
    ADD_PAIR(result, 2, 4, 6, 8)
    ADD_PAIR(result, 1, 3, 7, 10)
    ADD_PAIR(result, 3, 6, 9, 10)
    ADD_PAIR(result, 0, 4, 7, 9)
    ADD_PAIR(result, 4, 6, 9, 10)
    ADD_PAIR(result, 1, 4, 6, 8)
    ADD_PAIR(result, 0, 2, 4, 6)

    result.ref_genome_length = 11;
    return result;
}

int main(int argc, char** argv) {
    qmcp::SequentialMaxFlowSolver seq_solver;
    seq_solver.import_reads(aos_reads_1());
    seq_solver.solve(5);
    // App app;

    // try {
    //     app.Parse(argc, argv);
    // } catch (const CLI::ParseError& e) {
    //     return app.Exit(e);
    // }

    // app.Solve();

    return EXIT_SUCCESS;
}
