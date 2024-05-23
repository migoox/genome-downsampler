add_library(sample-data-generator STATIC src/sample_data_generator.cpp)

target_include_directories(sample-data-generator PUBLIC ${CMAKE_SOURCE_DIR}/libs/sample-data-generator/include)
