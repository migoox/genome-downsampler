add_executable(dna-downsampler main.cu)

target_link_libraries(dna-downsampler PRIVATE CLI11::CLI11)
