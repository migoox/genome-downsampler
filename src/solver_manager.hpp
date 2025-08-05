#ifndef SOLVER_MANAGER_HPP
#define SOLVER_MANAGER_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "helpers.hpp"
#include "qmcp-solver/mcp_cpu_cost_scaling_solver.hpp"
#include "qmcp-solver/qmcp_cpu_cost_scaling_solver.hpp"
#include "qmcp-solver/quasi_mcp_cpu_max_flow_solver.hpp"
#include "qmcp-solver/quasi_mcp_cuda_max_flow_solver.hpp"
#include "qmcp-solver/solver.hpp"

class SolverManager {
   public:
    SolverManager() {
        solvers_map_.emplace("quasi-mcp-cpu", std::make_unique<qmcp::QuasiMcpCpuMaxFlowSolver>());
        solvers_map_.emplace("mcp-cpu", std::make_unique<qmcp::McpCpuCostScalingSolver>());
        solvers_map_.emplace("qmcp-cpu", std::make_unique<qmcp::QmcpCpuCostScalingSolver>());
#ifdef CUDA_ENABLED
        solvers_map_.emplace("quasi-mcp-cuda", std::make_unique<qmcp::QuasiMcpCudaMaxFlowSolver>());
#endif

        algorithms_names_ = helpers::get_names_from_map(solvers_map_);
    }

    qmcp::Solver& get(const std::string& solver_name) const {
        return *solvers_map_.at(solver_name);
    }

    bool contains(const std::string& solver_name) {
        return solvers_map_.find(solver_name) != solvers_map_.end();
    }

    const std::vector<std::string>& get_names() const { return algorithms_names_; }

   private:
    std::map<std::string, std::unique_ptr<qmcp::Solver>> solvers_map_;
    std::vector<std::string> algorithms_names_;
};

#endif
