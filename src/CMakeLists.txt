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

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Setting Debug flags for genome-downsampler")

    # Disable PIE to make XRay symbolization work
    set_target_properties(genome-downsampler PROPERTIES POSITION_INDEPENDENT_CODE OFF)
    target_link_options(genome-downsampler PRIVATE -no-pie)

    # Enable full XRay instrumentation
    target_compile_definitions(genome-downsampler PRIVATE XRAY_INSTRUMENTATION_ENABLED=1)
    target_compile_options(genome-downsampler PRIVATE
        -g
        -O0
        -fxray-instrument
        -fno-omit-frame-pointer
    )

    # Set CMake global debug flags for consistency (optional)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -fxray-instrument -fno-omit-frame-pointer")
    set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G -g -O0 --extended-lambda")
else()
    message(STATUS "Setting Release flags for genome-downsampler")
    target_compile_options(genome-downsampler PRIVATE -O3)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -O3 --extended-lambda")
endif()

target_link_libraries(
    genome-downsampler
    PRIVATE CLI11::CLI11
            bam-api
            qmcp-solver
            reads-gen
            logging)
