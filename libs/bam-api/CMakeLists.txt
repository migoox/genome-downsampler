add_library(bam-api STATIC src/bam_api.cpp)

target_include_directories(bam-api PUBLIC ${CMAKE_SOURCE_DIR}/libs/bam-api/include)
