add_executable(cuda-dna-downsampler main.cu)

target_link_libraries(cuda-dna-downsampler PRIVATE CLI11::CLI11 bam-api qmcp-solver)
