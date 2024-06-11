add_library(
    qmcp-solver STATIC
    src/sequential_cost_scaling_network_solver.cpp
    src/sequential_max_flow_solver.cpp
    src/test_solver.cpp
    src/cuda_max_flow_solver.cu
    src/linear_programming_solver.cpp)

target_include_directories(qmcp-solver PUBLIC ${CMAKE_SOURCE_DIR}/libs/qmcp-solver/include)

target_link_libraries(
    qmcp-solver
    PRIVATE Boost::boost
            ortools
            bam-api
            logging)
