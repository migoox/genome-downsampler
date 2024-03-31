#ifndef QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#define QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#include "solver.hpp"

namespace qmcp {

class SequenceNetworkSolver : public Solver {
   public:
    SequenceNetworkSolver(/* Algorithm input data */) {}
    void solve() override;
};

}  // namespace qmcp

#endif
