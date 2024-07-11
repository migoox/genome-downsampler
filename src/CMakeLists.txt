# Add preprocessing definitions and optional subdirectories
if(WITH_TESTS)
    set(TESTING_ENABLED 1)
    add_subdirectory(tests)
endif()

if(WITH_CUDA)
    set(CUDA_ENABLED 1)
endif()

configure_file("${PROJECT_SOURCE_DIR}/src/config.h.in" "${PROJECT_BINARY_DIR}/config/config.h")
include_directories("${PROJECT_BINARY_DIR}/config")

set(SOURCES ${SOURCES} main.cpp app.cpp)

add_executable(genome-downsampler ${SOURCES})

target_link_libraries(
    genome-downsampler
    PRIVATE CLI11::CLI11
            bam-api
            qmcp-solver
            reads-gen
            logging)
