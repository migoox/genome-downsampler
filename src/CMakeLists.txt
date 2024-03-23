add_executable(dna-downsampler main.cpp)

target_link_libraries(
    dna-downsampler 
    PRIVATE CLI11::CLI11
)