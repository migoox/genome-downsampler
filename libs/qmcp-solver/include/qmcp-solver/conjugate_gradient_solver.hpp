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
    std::vector<double> create_b_vector(uint32_t M);
    void find_pairs(bool flag);
    void set_reads(const bam_api::SOAPairedReads& input_sequence);
    const std::vector<bam_api::BAMReadId>& get_output();
    std::vector<bam_api::BAMReadId> output_sequence();

   private:
    bam_api::SOAPairedReads paired_reads_;
    std::vector<bam_api::BAMReadId> solution_;
    std::filesystem::path input_;
    bool find_pairs_ = true;
};
}  // namespace qmcp
