cmake_minimum_required(VERSION 3.9)

project(physics LANGUAGES C)

add_library(physics STATIC)

target_sources(physics PRIVATE
    src/physics.c
)

target_include_directories(physics PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/include
    ${CMAKE_CURRENT_LIST_DIR}/../vector/include
    ${CMAKE_CURRENT_LIST_DIR}/../logging/include
    ${CMAKE_CURRENT_LIST_DIR}/../string/include
)

target_link_libraries(physics
    vector
)
