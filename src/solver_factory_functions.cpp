#include "solver_factory_functions.hpp"

#include <memory>

#include "bam-api/bam_api.hpp"
#include "qmcp-solver/cuda_max_flow_solver.hpp"
#include "qmcp-solver/sequential_cost_scaling_network_solver.hpp"
#include "qmcp-solver/sequential_max_flow_solver.hpp"
#include "qmcp-solver/solver.hpp"
#include "qmcp-solver/test_solver.hpp"

std::unique_ptr<qmcp::Solver> solver_factory_functions::createTestSolver(bam_api::BamApi& bam_api) {
    return std::make_unique<qmcp::TestSolver>(bam_api);
}

std::unique_ptr<qmcp::Solver> solver_factory_functions::createSequentialCostScalingNetworkSolver(
    bam_api::BamApi& bam_api) {
    return std::make_unique<qmcp::SequentialCostScalingNetworkSolver>(bam_api);
}

std::unique_ptr<qmcp::Solver> solver_factory_functions::createCudaMaxFlowSolver(
    bam_api::BamApi& bam_api) {
    return std::make_unique<qmcp::CudaMaxFlowSolver>(bam_api);
}

std::unique_ptr<qmcp::Solver> solver_factory_functions::createSequentialMaxFlowSolver(
    bam_api::BamApi& bam_api) {
    return std::make_unique<qmcp::SequentialMaxFlowSolver>(bam_api);
}
