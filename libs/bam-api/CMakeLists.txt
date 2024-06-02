add_library(bam-api STATIC
  src/bam_api.cpp
  src/read.cpp
  src/paired_reads.cpp src/aos_paired_reads.cpp src/soa_paired_reads.cpp
  src/amplicon.cpp src/amplicon_set.cpp
  src/bam_api_config_builder.cpp
)

target_include_directories(bam-api PUBLIC ${CMAKE_SOURCE_DIR}/libs/bam-api/include)

# bam-api depends on htslib
target_link_libraries(bam-api PRIVATE ${HTSlib_LIBRARIES} logging)
