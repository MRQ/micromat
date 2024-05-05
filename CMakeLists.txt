project(micromat)

cmake_minimum_required(VERSION 3.3)

add_library(micromat INTERFACE)
target_include_directories(micromat INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

if(BUILD_TESTING)
	add_subdirectory(tests)
endif()
