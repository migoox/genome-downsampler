#ifndef QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#define QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#include "solver.hpp"
#include "network_graph.hpp"
#include "../../bam-api/include/bam-api/bam_paired_reads.hpp"

namespace qmcp {

class SequenceNetworkSolver : public Solver {
   public:
    SequenceNetworkSolver(unsigned int M, bam_api::AOSPairedReads& sequence) : M_(M), sequence_(sequence) {}
    void solve() override;
    private:
    bam_api::AOSPairedReads sequence_;
    unsigned int M_;
    
};

}  // namespace qmcp

#endif
