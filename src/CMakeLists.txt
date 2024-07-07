# Add preprocessing definitions and optional subdirectories
if(WITH_TESTS)
    add_definitions(-DTESTING_ENABLED)
    add_subdirectory(tests)
endif()

if(NO_CUDA)
    add_definitions(-DNO_CUDA)
endif()

set(SOURCES ${SOURCES} main.cpp app.cpp)

add_executable(genome-downsampler ${SOURCES})

target_link_libraries(
    genome-downsampler
    PRIVATE CLI11::CLI11
            bam-api
            qmcp-solver
            reads-gen
            logging)
