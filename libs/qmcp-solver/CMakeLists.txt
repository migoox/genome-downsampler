add_library(qmcp-solver STATIC src/sequence_network_solver.cpp)

target_include_directories(bam-api PUBLIC ${CMAKE_SOURCE_DIR}/libs/qmcp-solver/include)
