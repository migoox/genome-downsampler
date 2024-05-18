#include "qmcp-solver/test_solver.hpp"

#include <iostream>

#include "bam-api/bam_api.hpp"
#include "logging/log.hpp"

void qmcp::TestSolver::Import(const std::filesystem::path& input,
                              uint32_t min_seq_length, uint32_t min_seq_mapq) {
    input_ = input;
    LOG(logging::kDebug) << "Import, min_len: " << min_seq_length
              << ", min_mapq: " << min_seq_mapq << ", input: " << input;

    paired_reads_ =
        bam_api::BamApi::read_bam_soa(input, min_seq_length, min_seq_mapq);

    LOG(logging::kInfo) << paired_reads_.ids.size() << " sequences has been imported!";
}

void qmcp::TestSolver::Solve(uint32_t max_coverage) {
    LOG(logging::kDebug) << "Solve (max_coverage set to " << max_coverage << ")";

    solution_ = std::vector<bam_api::ReadIndex>(
        paired_reads_.ids.end() - kReadsFromEndOffset,
        paired_reads_.ids.end() - kReadsFromEndOffset + kReadsCount);

    LOG(logging::kInfo) << "Solution have " << solution_.size() << " sequences!";
}

void qmcp::TestSolver::Export(const std::filesystem::path& output) {
    LOG(logging::kDebug) << "Output: " << output;

    uint32_t reads_written = 0;
    reads_written = bam_api::BamApi::write_bam(input_, output, solution_);

    LOG(logging::kInfo) << reads_written << " sequences has been written to file " << output;
}
