add_library(logging STATIC src/log.cpp)

target_include_directories(logging PUBLIC ${CMAKE_SOURCE_DIR}/libs/logging/include)
