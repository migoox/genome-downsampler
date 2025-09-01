#include "qmcp-solver/cpu_mixed_integer_solver.hpp"

#include <ortools/sat/cp_model.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <future>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "bam-api/read.hpp"
#include "qmcp-solver/qmcp_cpu_cost_scaling_solver.hpp"
#include "qmcp-solver/solver.hpp"

qmcp::CpuMixedIntegerSolver::CpuMixedIntegerSolver() {
    helper_solver = std::make_unique<qmcp::QmcpCpuCostScalingSolver>();
}

uint32_t qmcp::CpuMixedIntegerSolver::get_solution_reads_count(uint32_t max_coverage,
                                                               bam_api::BamApi& bam_api) {
    std::unique_ptr<qmcp::Solution> solution_helper = helper_solver->solve(max_coverage, bam_api);
    return solution_helper->size();
}

// OR - TOOLS NO ITERATION OVER WHOLE SEQUENCE - NOT ENOUGHT MEMORY
std::unique_ptr<qmcp::Solution> qmcp::CpuMixedIntegerSolver::solve(uint32_t max_coverage,
                                                                   bam_api::BamApi& bam_api) {
    using namespace operations_research;
    using namespace sat;

    const bam_api::AOSPairedReads& input_sequence = bam_api.get_paired_reads_aos();
    const uint32_t seq_size = input_sequence.ref_genome_length;
    const uint32_t reads_count = input_sequence.get_reads_count();

    CpModelBuilder model;
    int count = 0;
    std::vector<int> len;
    std::vector<BoolVar> select;

    select.reserve(reads_count);
    for (int i = 0; i < reads_count; ++i) select.push_back(model.NewBoolVar());

    std::vector<std::vector<BoolVar>> covering_reads;
    covering_reads.resize(seq_size);
    for (int x = 0; x < seq_size; x++) covering_reads[x].reserve(reads_count);

    std::vector<int> coverage_count(seq_size, 0);
    for (int nr = 0; nr < reads_count; nr++)
        for (uint32_t pos = input_sequence.reads[nr].start_ind;
             pos <= input_sequence.reads[nr].end_ind; pos++) {
            covering_reads[pos].push_back(select[nr]);
            coverage_count[pos]++;
        }

    for (int pos = 0; pos < seq_size; pos++)
        if (!covering_reads.empty()) {
            const int rhs = std::min(static_cast<int>(max_coverage), coverage_count[pos]);
            model.AddGreaterOrEqual(LinearExpr::Sum(covering_reads[pos]), rhs);
            covering_reads[pos].clear();
            covering_reads[pos].shrink_to_fit();
        }
    covering_reads.clear();
    covering_reads.shrink_to_fit();

    uint32_t reads_in_solution = get_solution_reads_count(max_coverage, bam_api);
    std::vector<uint32_t> read_lengths(reads_count);
    for (int i = 0; i < reads_count; ++i)
        read_lengths[i] = input_sequence.reads[i].end_ind - input_sequence.reads[i].start_ind + 1;

    model.AddEquality(LinearExpr::Sum(select), reads_in_solution);

    LinearExpr total_length;
    for (int i = 0; i < reads_count; ++i) total_length += select[i] * read_lengths[i];

    model.Minimize(total_length);

    const CpSolverResponse response = Solve(model.Build());

    if (response.status() != CpSolverStatus::OPTIMAL &&
        response.status() != CpSolverStatus::FEASIBLE)
        throw std::runtime_error("Didn't solve for optimal");

    std::vector<int64_t> ret_vec(reads_count, 0.0);
    for (int i = 0; i < reads_count; ++i) ret_vec[i] = SolutionIntegerValue(response, select[i]);

    return obtain_sequence(input_sequence, ret_vec);
}

std::unique_ptr<qmcp::Solution> qmcp::CpuMixedIntegerSolver::obtain_sequence(
    const bam_api::AOSPairedReads& sequence, std::vector<int64_t>& solution) {
    auto reduced_reads = std::make_unique<Solution>();
    for (bam_api::ReadIndex read_id = 0; read_id < sequence.get_reads_count(); read_id++)
        if (solution[static_cast<int>(read_id)] > 0) reduced_reads->push_back(read_id);

    return reduced_reads;
}