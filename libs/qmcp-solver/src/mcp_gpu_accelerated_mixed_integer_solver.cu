// #include <__clang_cuda_cmath.h>
// #include <lpsolve/lp_lib.h>

// #include <algorithm>
// #include <cstdint>
// #include <exception>
// #include <iostream>
// #include <memory>
// #include <stdexcept>
// #include <vector>

// #include "bam-api/bam_api.hpp"
// #include "bam-api/read.hpp"
// #include "qmcp-solver/cuda_helpers.cuh"
// #include "qmcp-solver/mcp_cpu_mixed_integer_solver.hpp"
// #include "qmcp-solver/solver.hpp"

// void printReads(const bam_api::SOAPairedReads& input_sequence) {
//     for (int x = 0; x < input_sequence.get_reads_count(); x++) {
//         int y = 0;
//         for (y = 0; y < input_sequence.start_inds[x]; y++) std::cout << "0";
//         for (; y <= input_sequence.end_inds[x]; y++) std::cout << "1";
//         for (; y < input_sequence.ref_genome_length; y++) std::cout << "0";
//         std::cout << "\n";
//     }
// }

// template <typename T>
// void printVector(std::vector<T>& v) {
//     for (auto e : v) std::cout << e << " ";
//     std::cout << "\n";
// }

// __constant__ uint32_t gpu_max_coverage;
// __constant__ uint32_t gpu_seq_size;
// __constant__ uint32_t gpu_reads_count;

// __global__ void prepareMatrixKernel(const bam_api::Index* start_inds,
//                                     const bam_api::Index* end_inds, double** rows, int* b) {
//     int num_of_ones = 0;
//     uint32_t x = blockIdx.x * blockDim.x + threadIdx.x;
//     if (x >= gpu_seq_size) return;
//     num_of_ones = 0;
//     for (int i = 0; i < gpu_reads_count; i++)
//         if (start_inds[i] <= x && end_inds[i] >= x) {
//             rows[x][i + 1] = 1;
//             num_of_ones++;
//         } else
//             rows[i + 1] = 0;
//     b[x] = std::min(static_cast<int>(gpu_max_coverage), num_of_ones);
// }

// // TODO why rows as double
// void prepareMatrix(uint32_t max_coverage, bam_api::BamApi& bam_api, double** rows_out, int*
// b_out) {

// }

// std::unique_ptr<qmcp::Solution> qmcp::McpCpuMixedIntegerSolver::solve(uint32_t max_coverage,
//                                                                       bam_api::BamApi& bam_api) {
//     const bam_api::SOAPairedReads& input_sequence = bam_api.get_paired_reads_soa();
//     // printReads(input_sequence);
//     uint32_t reads_count = input_sequence.get_reads_count();
//     uint32_t seq_size = input_sequence.ref_genome_length;
//     lprec* model = make_lp(0, static_cast<int>(reads_count));
//     if (model == nullptr) throw std::runtime_error("Couldn't construct new model");

//     for (int x = 1; x <= reads_count; x++) set_binary(model, x, TRUE);

//     std::vector<double> row(reads_count + 1, 0);
//     std::vector<double> obj_fn(reads_count + 1, 1);
//     set_obj_fn(model, obj_fn.data());
//     set_minim(model);

//     int num_of_ones = 0;
//     for (int x = 0; x < seq_size; x++) {
//         num_of_ones = 0;
//         for (int i = 0; i < reads_count; i++)
//             if (input_sequence.start_inds[i] <= x && input_sequence.end_inds[i] >= x) {
//                 row[i + 1] = 1;
//                 num_of_ones++;
//             } else
//                 row[i + 1] = 0;
//         int b = std::min(static_cast<int>(max_coverage), num_of_ones);

//         if (b > 0) add_constraint(model, row.data(), GE, b);
//     }
//     // print_lp(model);
//     if (::solve(model) != OPTIMAL) throw std::runtime_error("Didn't solve for optimal");

//     get_variables(model, row.data());
//     return obtain_sequence(input_sequence, row);
// }

// std::unique_ptr<qmcp::Solution> qmcp::McpCpuMixedIntegerSolver::obtain_sequence(
//     const bam_api::SOAPairedReads& sequence, std::vector<double>& solution) {
//     auto reduced_reads = std::make_unique<Solution>();

//     for (bam_api::ReadIndex read_id = 0; read_id < sequence.get_reads_count(); ++read_id)
//         if (solution[static_cast<int>(read_id)] > 0) reduced_reads->push_back(read_id);
//     // printVector(solution);
//     std::cout << "solved\n";
//     return reduced_reads;
// }
