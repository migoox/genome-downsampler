add_library(bicgstab STATIC  src/bicgstab.cpp)

target_include_directories(bicgstab PUBLIC ${CMAKE_SOURCE_DIR}/libs/bicgstab/include)

target_link_libraries(
    bicgstab
    PRIVATE cudart cublas cusparse logging)