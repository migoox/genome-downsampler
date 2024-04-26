#include "../include/qmcp-solver/sequential_cost_scaling_network_solver.hpp"

#include <ortools/graph/max_flow.h>

#include <iostream>

void qmcp::SequentialCostScalingNetworkSolver::solve() {
    std::cout << "Not implemented!";
    // TODO(implement the function):
    // ortools max flow example basing on the link:
    // https://developers.google.com/optimization/flow/maxflow#c++
    // Worth reading:
    // https://github.com/google/or-tools/blob/stable/ortools/graph/min_cost_flow.h
    operations_research::SimpleMaxFlow max_flow;
    std::vector<int64_t> start_nodes = {0, 0, 0, 1, 1, 2, 2, 3, 3};
    std::vector<int64_t> end_nodes = {1, 2, 3, 2, 4, 3, 4, 2, 4};
    std::vector<int64_t> capacities = {20, 30, 10, 40, 30,  // NOLINT
                                       10, 20, 5,  20};     // NOLINT

    for (int i = 0; i < start_nodes.size(); ++i) {
        max_flow.AddArcWithCapacity(start_nodes[i], end_nodes[i],  // NOLINT
                                    capacities[i]);
    }
    int status = max_flow.Solve(0, 4);
    if (status == operations_research::SimpleMaxFlow::Status::OPTIMAL) {
        LOG(INFO) << "Max flow: " << max_flow.OptimalFlow();
        LOG(INFO) << "";
        LOG(INFO) << "  Arc    Flow / Capacity";
        for (std::size_t i = 0; i < max_flow.NumArcs(); ++i) {
            int i_casted = static_cast<int>(i);
            LOG(INFO) << max_flow.Tail(i_casted) << " -> "
                      << max_flow.Head(i_casted) << "  "
                      << max_flow.Flow(i_casted) << "  / "
                      << max_flow.Capacity(i_casted);
        }
    } else {
        LOG(INFO) << "Solving the max flow problem failed. Solver status: "
                  << status;
    }
}