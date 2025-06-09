cmake_minimum_required(VERSION 3.10)

project(physics LANGUAGES C)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(physics STATIC)

target_sources(physics PRIVATE
    physics.c
)

target_include_directories(physics PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/../string
    ${CMAKE_CURRENT_LIST_DIR}/../logging
    ${CMAKE_CURRENT_LIST_DIR}/../vector
)

target_link_libraries(physics
    vector
)
