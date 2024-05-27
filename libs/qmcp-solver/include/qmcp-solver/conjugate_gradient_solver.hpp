#pragma once

#include <cstdint>

#include "bam-api/bam_paired_reads.hpp"
#include "qmcp-solver/solver.hpp"

namespace qmcp {
class ConjugateGradientSolver : public Solver {
   public:
    void import_reads(const std::filesystem::path& input, uint32_t min_seq_length,
                      uint32_t min_seq_mapq) override;
    void solve(uint32_t max_coverage) override;
    void export_reads(const std::filesystem::path& output) override;

    void make_matrix(int* n_out, int** row_offsets_out, int** columns_out, double** values_out);
    std::vector<double> create_x_vector(uint32_t M);

   private:
    bam_api::SOAPairedReads paired_reads_;
    std::vector<bam_api::ReadIndex> solution_;
    std::filesystem::path input_;
};
}  // namespace qmcp
