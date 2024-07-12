cmake_minimum_required(VERSION 3.21)

option(WITH_TESTS "enables unit tests compilation" OFF)
option(WITH_CUDA "enables CUDA algorithms" ON)

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

# Using clang-tidy this way slows down the compilation
# set(CMAKE_CXX_CLANG_TIDY "clang-tidy")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set the project name and language

if(WITH_CUDA)
    project(genome-downsampler LANGUAGES CXX CUDA)
else()
    project(genome-downsampler LANGUAGES CXX)
endif()

# Uncomment if the repositories included with CPM have some dependencies to build
# set(BUILD_DEPS ON)

# Dependencies with CPM
#   https://github.com/cpm-cmake/CPM.cmake/blob/master/README.md
include(cmake/CPM.cmake)

cpmaddpackage("gh:CLIUtils/CLI11@2.4.1")
cpmaddpackage(
    NAME
    Ccache.cmake
    GITHUB_REPOSITORY
    TheLartians/Ccache.cmake
    VERSION
    1.2)

if(UNIX)
    # Add popular cmake configurations paths on UNIX systems,
    # these directories will be searched when the find_package looks
    # for the <Package>Config.cmake files.
    set(CMAKE_PREFIX_PATH
        ${CMAKE_PREFIX_PATH}
        /usr/lib/x86_64-linux-gnu/cmake/
        /usr/lib/cmake/
        /usr/local/lib/cmake/)
endif()

# Set additional Find<Package>.cmake files direction
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

# Display CMAKE_PREFIX_PATH entries
message("CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")

if(DEFINED ENV{ZLIB_INCLUDE_DIR})
    set(ZLIB_INCLUDE_DIR $ENV{zlib_INCLUDE_DIR})
elseif(UNIX)
    set(ZLIB_INCLUDE_DIR /usr/include/ /usr/local/include/)
endif()

if(DEFINED ENV{ZLIB_LIBRARYDIR})
    set(ZLIB_LIBRARY $ENV{zlib_LIBRARYDIR})
elseif(UNIX)
    set(ZLIB_LIBRARY /usr/lib/ /usr/lib/x86_64-linux-gnu/ /usr/local/lib/)
endif()

# Add ZLIB Library (htslib depends on it)
find_package(ZLIB REQUIRED)

# Add htslib Library
find_package(HTSlib REQUIRED)

# Find the ccache program
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

# Find the ortools package by pure ortoolsConfig.cmake search
find_package(ortools CONFIG REQUIRED)

# Use ccache to speed up the build process
set(USE_CCACHE YES)

# Subdirectories
add_subdirectory(src)
add_subdirectory(libs)
