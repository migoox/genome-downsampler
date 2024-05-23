add_library(cg STATIC  src/cg.cpp)

target_include_directories(cg PUBLIC ${CMAKE_SOURCE_DIR}/libs/cg/include)

target_link_libraries(
    cg
    PRIVATE cudart cublas cusparse logging)
