#ifndef SOLVER_FACTORY_FUNCTIONS_HPP
#define SOLVER_FACTORY_FUNCTIONS_HPP

#include <memory>
#include "bam-api/bam_api.hpp"
#include "qmcp-solver/solver.hpp"

namespace solver_factory_functions {

std::unique_ptr<qmcp::Solver> createTestSolver(bam_api::BamApi& bam_api);
std::unique_ptr<qmcp::Solver> createSequentialCostScalingNetworkSolver(bam_api::BamApi& bam_api);
std::unique_ptr<qmcp::Solver> createCudaMaxFlowSolver(bam_api::BamApi& bam_api);
std::unique_ptr<qmcp::Solver> createSequentialMaxFlowSolver(bam_api::BamApi& bam_api);

}  // namespace solver_factory_functions

#endif
