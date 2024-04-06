cmake_minimum_required(VERSION 3.21)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

# Set the target architecture library directory name, if detected.
set(CMAKE_CXX_LIBRARY_ARCHITECTURE x86_64-linux-gnu)
set(CMAKE_CUDA_LIBRARY_ARCHITECTURE x86_64-linux-gnu)

# Enable this globally to avoid conflicts between
# -Wpedantic being enabled and -std=c++20 and -std=gnu++20
# for example when compiling with PCH enabled
set(CMAKE_CXX_EXTENSIONS OFF)

# The build process is kept in a different directory from
# the main source code
set(CMAKE_DISABLE_IN_SOURCE_BUILD TRUE)

# Using clang-tidy this way slows down the compilation
# set(CMAKE_CXX_CLANG_TIDY "clang-tidy")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

# Use '.a' files (static libraries)
set(Boost_USE_STATIC_LIBS OFF)

# If find_package(Boost) doesn't work, environmental variables
# my be set in order to help the cmake
if(DEFINED ENV{Boost_INCLUDE_DIR})
    set(Boost_INCLUDE_DIR $ENV{Boost_INCLUDE_DIR})
endif()
if(DEFINED ENV{Boost_LIBRARYDIR})
    set(Boost_LIBRARYDIR $ENV{Boost_LIBRARYDIR})
endif()
find_package(Boost 1.74 REQUIRED)

# Find the ccache program
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

# Use ccache to speed up the build process
set(USE_CCACHE YES)
# Set the project name and language
project(cuda-dna-downsampler LANGUAGES CXX CUDA)

# Subdirectories
add_subdirectory(src)
add_subdirectory(libs)
