#include <htslib/hts.h>
#include <stdio.h>

#include <chrono>
#include <cstdlib>

#include "app.hpp"
#include "device_launch_parameters.h"
#include "qmcp-solver/cuda_max_flow_solver.hpp"
#include "qmcp-solver/qmcp-solver.hpp"
#include "qmcp-solver/sequential_cost_scaling_network_solver.hpp"

int main(int argc, char** argv) {
    App app;

    try {
        app.Parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.Exit(e);
    }

    app.Solve();

    return EXIT_SUCCESS;
}
