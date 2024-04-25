add_library(qmcp-solver STATIC src/sequence_network_solver.cpp)

target_include_directories(qmcp-solver PUBLIC ${CMAKE_SOURCE_DIR}/libs/qmcp-solver/include)

target_link_libraries(qmcp-solver PRIVATE Boost::boost ortools)
