add_library(qmcp-solver STATIC src/sequential_cycle_canceling_network_solver.cpp src/sequential_cost_scaling_network_solver.cpp)

target_include_directories(qmcp-solver PUBLIC ${CMAKE_SOURCE_DIR}/libs/qmcp-solver/include)

target_link_libraries(qmcp-solver PRIVATE Boost::boost ortools)
