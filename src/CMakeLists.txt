add_subdirectory(tests)
set(SOURCES
    ${SOURCES}
    main.cpp
    app.cpp
)

add_executable(cuda-dna-downsampler ${SOURCES})

target_link_libraries(
    cuda-dna-downsampler
    PRIVATE CLI11::CLI11
            bam-api
            qmcp-solver
            reads-gen
            logging)
