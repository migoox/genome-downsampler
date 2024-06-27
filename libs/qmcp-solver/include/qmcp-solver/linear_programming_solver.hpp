#ifndef QMCP_LINEAR_PROGRAMMING_SOLVER_HPP
#define QMCP_LINEAR_PROGRAMMING_SOLVER_HPP

#include <memory>
#include <vector>

#include "bam-api/read.hpp"
#include "bicgstab/bicgstab.hpp"
#include "logging/log.hpp"
#include "solver.hpp"

namespace qmcp {
class LinearProgrammingSolver : public Solver {
   public:
    std::unique_ptr<Solution> solve(uint32_t max_coverage, bam_api::BamApi& bam_api) override;
    bool uses_quality_of_reads() override { return true; }

   private:
    void make_matrix(int* n_out, int** row_offsets_out, int** columns_out, double** values_out);
    std::vector<double> create_b_vector(uint32_t M);

    std::vector<double> process_bicgstab(int rows, int* row_offsets_out, int* columns_out,
                                         double* values_out, std::vector<double> b,
                                         std::vector<double> x);

    bam_api::SOAPairedReads input_sequence_;
};
}  // namespace qmcp
#endif
