cmake_minimum_required(VERSION 3.10)
include_guard(GLOBAL)

project(physics LANGUAGES C)

find_package(GSL REQUIRED)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(physics STATIC)
add_subdirectory(../string ${CMAKE_BINARY_DIR}/string_from_physics)
add_subdirectory(../logging ${CMAKE_BINARY_DIR}/logging_from_physics)
add_subdirectory(../vector ${CMAKE_BINARY_DIR}/vector_from_physics)

target_sources(physics PRIVATE
	physics.c
)

target_include_directories(physics PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(physics
	PUBLIC
		string
		logging
		vector
	PRIVATE
		${GSL_LIBRARIES}
		m
)

target_include_directories(physics PRIVATE
	${GSL_INCLUDE_DIRS}
)
