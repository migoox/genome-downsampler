cmake_minimum_required(VERSION 3.21)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

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

# Set the project name and language
project(cuda-dna-downsampler LANGUAGES CXX CUDA)

# Dependencies with CPM
#   https://github.com/cpm-cmake/CPM.cmake/blob/master/README.md
include(cmake/CPM.cmake)

cpmaddpackage("gh:CLIUtils/CLI11@2.4.1")

# Subdirectories
add_subdirectory(src)
add_subdirectory(libs)
