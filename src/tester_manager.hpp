#ifndef TESTER_MANAGER_HPP
#define TESTER_MANAGER_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "helpers.hpp"
#include "tests/coverage_tester.hpp"
#include "tests/solver_tester.hpp"

class TesterManager {
   public:
    TesterManager() {
        solver_testers_map_.emplace("coverage", std::make_unique<test::CoverageTester>());

        testers_names_ = helpers::get_names_from_map(solver_testers_map_);
    }

    test::SolverTester& get(const std::string& tester_name) {
        return *solver_testers_map_.at(tester_name);
    }

    bool contains(const std::string& tester_name) {
        return solver_testers_map_.find(tester_name) != solver_testers_map_.end();
    }

    const std::vector<std::string>& get_names() const { return testers_names_; }

   private:
    std::map<std::string, std::unique_ptr<test::SolverTester>> solver_testers_map_;
    std::vector<std::string> testers_names_;
};

#endif
