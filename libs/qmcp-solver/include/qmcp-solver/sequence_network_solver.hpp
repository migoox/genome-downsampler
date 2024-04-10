#ifndef QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#define QMCP_SEQUENCE_NETWORK_SOLVER_HPP
#include "solver.hpp"
#include "../../bam-api/include/bam-api/bam_api.hpp"

namespace qmcp {

class SequenceNetworkSolver : public Solver {
   public:
    SequenceNetworkSolver(const bam_api::BamSequence& sequence) : sequence(sequence) {}
    void solve() override;
    private:
    bam_api::BamSequence sequence;
    
};

}  // namespace qmcp

#endif
