#ifndef QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#define QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#include "solver.hpp"
#include "network_graph.hpp"
#include "../../bam-api/include/bam-api/bam_api.hpp"

namespace qmcp {

class SequenceNetworkSolver : public Solver {
   public:
    SequenceNetworkSolver(const bam_api::BamSequence& sequence,unsigned int M) : sequence(sequence), M(M) {}
    void solve() override;
    private:
    bam_api::BamSequence sequence;
    unsigned int M;
    
};

}  // namespace qmcp

#endif
