add_executable(dna-downsampler main.cu)

target_link_libraries(dna-downsampler PRIVATE CLI11::CLI11)

target_include_directories(dna-downsampler PRIVATE ${CMAKE_SOURCE_DIR}/libs/bam-api/include
                                                   ${CMAKE_SOURCE_DIR}/libs/qmcp-solver/include)
