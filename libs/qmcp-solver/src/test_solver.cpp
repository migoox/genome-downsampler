#include "qmcp-solver/test_solver.hpp"

#include <memory>

#include "bam-api/bam_api.hpp"
#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/solver.hpp"

std::unique_ptr<qmcp::Solution> qmcp::TestSolver::solve(uint32_t max_coverage,
                                                        bam_api::BamApi& bam_api) {
    return std::make_unique<Solution>(bam_api.get_paired_reads_soa().ids);
}
