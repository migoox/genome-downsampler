#include "qmcp-solver/test_solver.hpp"

#include <iostream>

#include "bam-api/bam_api.hpp"
#include "logging/log.hpp"

void qmcp::TestSolver::import_reads(const std::filesystem::path& input, uint32_t min_seq_length,
                                    uint32_t min_seq_mapq) {
    input_ = input;
    LOG_WITH_LEVEL(logging::kDebug) << "Import, min_len: " << min_seq_length
                                    << ", min_mapq: " << min_seq_mapq << ", input: " << input;

    paired_reads_ = bam_api::BamApi::read_bam_soa(input, min_seq_length, min_seq_mapq);

    LOG_WITH_LEVEL(logging::kInfo) << paired_reads_.ids.size() << " sequences has been imported!";
}

void qmcp::TestSolver::solve(uint32_t max_coverage) {
    LOG_WITH_LEVEL(logging::kDebug) << "Solve (max_coverage set to " << max_coverage << ")";

    solution_ = paired_reads_.ids;

    LOG_WITH_LEVEL(logging::kInfo) << "Solution have " << solution_.size() << " sequences!";
}

void qmcp::TestSolver::export_reads(const std::filesystem::path& output) {
    LOG_WITH_LEVEL(logging::kDebug) << "Output: " << output;

    uint32_t reads_written = 0;
    reads_written = bam_api::BamApi::write_bam(input_, output, solution_);

    LOG_WITH_LEVEL(logging::kInfo)
        << reads_written << " sequences has been written to file " << output;
}
