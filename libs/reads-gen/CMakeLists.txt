add_library(reads-gen STATIC src/reads_gen.cpp)

target_include_directories(reads-gen PUBLIC ${CMAKE_SOURCE_DIR}/libs/reads-gen/include)
target_link_libraries(reads-gen PRIVATE bam-api)
target_compile_options(${PROJECT_NAME} PRIVATE ${PROJ_CXX_FLAGS})
