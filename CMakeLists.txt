cmake_minimum_required(VERSION 3.21)

option(WITH_TESTS "Enables unit tests compilation" OFF)
option(WITH_CUDA "Enables CUDA algorithms" OFF)

set(CMAKE_CXX_STANDARD 17)
if(WITH_CUDA)
    set(CMAKE_CUDA_STANDARD 17)
endif()

# Set the target architecture library directory name, if detected.
set(CMAKE_CXX_LIBRARY_ARCHITECTURE x86_64-linux-gnu)
if(WITH_CUDA)
    set(CMAKE_CUDA_LIBRARY_ARCHITECTURE x86_64-linux-gnu)
endif()

# Enable this globally to avoid conflicts between
# -Wpedantic being enabled and -std=c++20 and -std=gnu++20
# for example when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

# The build process is kept in a different directory from
# the main source code, it is kept in ./build/
set(CMAKE_DISABLE_IN_SOURCE_BUILD TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(WITH_CUDA)
    project(genome-downsampler LANGUAGES CXX CUDA)
else()
    project(genome-downsampler LANGUAGES CXX)
endif()

if(UNIX)
    # Add popular cmake configurations paths on UNIX systems,
    # these directories will be searched when the find_package looks
    # for the <Package>Config.cmake files.
    set(CMAKE_PREFIX_PATH
        ${CMAKE_PREFIX_PATH}
        /usr/lib/x86_64-linux-gnu/cmake/
        /usr/lib/cmake/
        /usr/local/lib/cmake/
    )
else()
    message(FATAL_ERROR "Unsupported platform")
endif()

# Set additional Find<Package>.cmake files direction
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

find_package(HTSlib REQUIRED)
find_package(ortools CONFIG REQUIRED)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/deps_fetcher.cmake)
fetch_source_CLI11()

add_subdirectory(src)
add_subdirectory(libs)

include(cmake/generate_clangd_file.cmake)
generate_clangd_file("${PROJ_CXX_FLAGS}")
