option(WITH_TESTS "enables unit tests compilation" OFF)

if(WITH_TESTS)
    add_definitions(-DTESTING_ENABLED)
    add_subdirectory(tests)
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
