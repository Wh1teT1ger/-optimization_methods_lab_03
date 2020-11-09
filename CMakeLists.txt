cmake_minimum_required(VERSION 3.16)
project(simplex_method)
project(simplex_method_lib)

set(CMAKE_CXX_STANDARD 17)

add_library(branch_and_bound_method_lib STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/sources/branch_and_bound_method.cpp
        )

add_executable(branch_and_bound_method sources/main.cpp)

target_include_directories(branch_and_bound_method
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
        )
target_include_directories(branch_and_bound_method_lib
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
        )

target_link_libraries(branch_and_bound_method branch_and_bound_method_lib)