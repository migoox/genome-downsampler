#ifndef QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#define QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#include "solver.hpp"
#include "network_graph.hpp"
#include "../../bam-api/include/bam-api/bam_paired_reads.hpp"

namespace qmcp {

class SequenceNetworkSolver : public Solver {
   public:
    SequenceNetworkSolver(unsigned int M) : M(M) {}
    void solve() override;
    private:
    bam_api::AOSPairedReads sequence;
    unsigned int M;
    
};

}  // namespace qmcp

#endif
