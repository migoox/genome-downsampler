add_library(
    qmcp-solver STATIC
    src/quasi_mcp_cpu_max_flow_solver.cpp
    src/qmcp_cpu_cost_scaling_solver.cpp
    src/mcp_cpu_cost_scaling_solver.cpp
    src/test_solver.cpp
    src/quasi_mcp_cuda_max_flow_solver.cu
    src/linear_programming_solver.cpp
    src/qmcp_linear_programming_solver.cpp)

target_include_directories(qmcp-solver PUBLIC ${CMAKE_SOURCE_DIR}/libs/qmcp-solver/include)

target_link_libraries(qmcp-solver PRIVATE ortools bam-api logging bicgstab)
