#include "qmcp-solver/test_solver.hpp"

#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"

std::vector<bam_api::BAMReadId> qmcp::TestSolver::solve(uint32_t max_coverage) {
    return bam_api_.get_paired_reads_soa().ids;
}
