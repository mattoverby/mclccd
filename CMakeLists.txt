# Copyright Matt Overby 2021.
# Distributed under the MIT License.
cmake_minimum_required(VERSION 3.0)

project(mclccd_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(DownloadProject)

option(MCL_DEBUG "Compile in test mode" OFF)
option(MCL_BUILD_EXAMPLES "Compile examples" ON)

# I know you should not set build type directly and instead let the user
# set it, but this is way more convenient.
if (MCL_DEBUG)
	add_definitions(-DMCL_DEBUG)
	set(CMAKE_BUILD_TYPE Debug)
	if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wundef -Wextra -Wno-long-long -Wlogical-op -Wshadow -Wmaybe-uninitialized")
	endif()
else()
	set(CMAKE_BUILD_TYPE Release)
endif()

# thread building blocks (CPU threading)
find_package(TBB REQUIRED)
include_directories(SYSTEM ${TBB_INCLUDE_DIRS})

# Eigen
find_package(EIGEN)
if(EIGEN_FOUND)
	include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})
else()
	download_project(PROJ eigen3
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_TAG master
		UPDATE_DISCONNECTED 1
		QUIET)
	include_directories(SYSTEM ${eigen3_SOURCE_DIR})
endif()

# Create the library
file(GLOB SRC
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/*.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/ccd_internal/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/ccd_internal/*.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/ccd_internal/*.cpp)
add_library(mclccd ${SRC})

if(MCL_BUILD_EXAMPLES)

	# Get libigl (and glad/glfw/eigen) from github using DownloadProject
	option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF) # off = header only
	option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
	option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
	option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
	option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
	option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
	option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
	option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
	option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
	option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         OFF)
	option(LIBIGL_WITH_XML               "Use XML"                      OFF)
	download_project(PROJ libigl
		GIT_REPOSITORY https://github.com/libigl/libigl.git
		GIT_TAG main
		UPDATE_DISCONNECTED 1
		QUIET)
	add_subdirectory(${libigl_SOURCE_DIR} ${libigl_BINARY_DIR})

	include_directories(SYSTEM ${glad_SOURCE_DIR}/include)
	include_directories(SYSTEM ${libigl_SOURCE_DIR}/include)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

	# Warning: some libigl classes don't play nice with newer Eigen
	# but it's best to use newest version when you can.
#	include_directories(SYSTEM ${libigl_SOURCE_DIR}/external/eigen)

	add_definitions(-DMCL_CCD_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

	# Examples
	add_executable(dillo ${CMAKE_CURRENT_SOURCE_DIR}/examples/dillo.cpp)
	target_link_libraries(dillo mclccd glad glfw ${TBB_LIBRARIES})

endif()
