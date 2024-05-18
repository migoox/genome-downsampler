#include "qmcp-solver/test_solver.hpp"

#include <iostream>

#include "bam-api/bam_api.hpp"

void qmcp::TestSolver::Import(const std::filesystem::path& input,
                              uint32_t min_seq_length, uint32_t min_seq_mapq) {
    input_ = input;
    std::cout << "Import, min_len: " << min_seq_length
              << ", min_mapq: " << min_seq_mapq << ", input: " << input
              << std::endl;

    paired_reads_ =
        bam_api::BamApi::read_bam_soa(input, min_seq_length, min_seq_mapq);

    std::cout << paired_reads_.ids.size() << " sequences has been imported!"
              << std::endl;
}

void qmcp::TestSolver::Solve(uint32_t max_coverage) {
    std::cout << "Solve (max_coverage set to " << max_coverage << ")"
              << std::endl;

    solution_ = std::vector<bam_api::ReadIndex>(
        paired_reads_.ids.end() - kReadsFromEndOffset,
        paired_reads_.ids.end() - kReadsFromEndOffset + kReadsCount);

    std::cout << "Solution have " << solution_.size() << " sequences!"
              << std::endl;
}

void qmcp::TestSolver::Export(const std::filesystem::path& output) {
    std::cout << "Output: " << output << std::endl;

    uint32_t reads_written = 0;
    reads_written = bam_api::BamApi::write_bam(input_, output, solution_);

    std::cout << reads_written << " sequences has been written to file "
              << output << std::endl;
}
