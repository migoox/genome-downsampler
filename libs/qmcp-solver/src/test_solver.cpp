#include "qmcp-solver/test_solver.hpp"

#include <memory>

#include "bam-api/bam_api.hpp"
#include "bam-api/read.hpp"
#include "bam-api/soa_paired_reads.hpp"
#include "qmcp-solver/solver.hpp"

std::unique_ptr<qmcp::Solution> qmcp::TestSolver::solve(uint32_t max_coverage,
                                                        bam_api::BamApi& bam_api) {
    const bam_api::SOAPairedReads& paired_reads = bam_api.get_paired_reads_soa();

    std::unique_ptr<Solution> solution = std::make_unique<Solution>();
    solution->reserve(paired_reads.get_reads_count());

    for (bam_api::ReadIndex i = 0; i < paired_reads.get_reads_count(); ++i) {
        solution->push_back(i);
    }

    return solution;
}
