#include <htslib/hts.h>
#include <stdio.h>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "bam-api/bam_api.hpp"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "qmcp-solver/cuda_max_flow_solver.hpp"
#include "qmcp-solver/qmcp-solver.hpp"
#include "qmcp-solver/sequential_cost_scaling_network_solver.hpp"

int main() {
    std::cout << "READING " << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    qmcp::CudaMaxFlowSolver solver = qmcp::CudaMaxFlowSolver(
        "/home/billyk/Downloads/gpu-programming/data/ESIB_EQA_2023.SARS2.01/reads.bam");
    auto stop = std::chrono::high_resolution_clock::now();

    float solve_duration =
        static_cast<float>(
            std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()) /
        1000.F;

    std::cout << "READ TOOK " << solve_duration << " [seconds]" << std::endl;

    std::cout << "SOLVING " << std::endl;
    start = std::chrono::high_resolution_clock::now();
    solver.solve(1000);
    stop = std::chrono::high_resolution_clock::now();

    solve_duration =
        static_cast<float>(
            std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()) /
        1000.F;

    std::cout << "SOLVE TOOK " << solve_duration << " [seconds]" << std::endl;

    std::cout << "EXPORTING " << std::endl;
    start = std::chrono::high_resolution_clock::now();
    solver.export_data(
        "/home/billyk/Downloads/gpu-programming/data/ESIB_EQA_2023.SARS2.01/solve_reads.bam");
    stop = std::chrono::high_resolution_clock::now();

    solve_duration =
        static_cast<float>(
            std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()) /
        1000.F;

    std::cout << "EXPORT TOOK " << solve_duration << " [seconds]" << std::endl;

    return 0;
}
