#ifndef QMCP_SOLVER_HPP
#define QMCP_SOLVER_HPP

#include <cstdint>
namespace qmcp {

class Solver {
   public:
    virtual ~Solver() = default;
    virtual void solve(uint32_t required_cover) = 0;

   private:
};

}  // namespace qmcp

#endif