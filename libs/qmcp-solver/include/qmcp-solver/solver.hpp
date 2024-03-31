#ifndef QMCP_SOLVER_HPP
#define QMCP_SOLVER_HPP

namespace qmcp {

class Solver {
   public:
    virtual ~Solver() = default;
    virtual void solve() = 0;

   private:
};

}  // namespace qmcp

#endif