# Copyright Matt Overby 2021.
# Distributed under the MIT License.
cmake_minimum_required(VERSION 3.0)
project(mclccd LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FetchContent)

option(MCLCCD_DEBUG "Compile in test mode" OFF)
option(MCLCCD_BUILD_EXAMPLES "Compile examples" ON)

if (MCLCCD_DEBUG)
	set(CMAKE_BUILD_TYPE Debug)
	if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wundef -Wextra -Wno-long-long -Wlogical-op -Wshadow -Wmaybe-uninitialized")
	endif()
else()
	set(CMAKE_BUILD_TYPE Release)
endif()

# thread building blocks
# TODO: build this with fetchcontent if not found
find_package(TBB REQUIRED)

# Eigen, download if not found
find_package(EIGEN)
if(NOT EIGEN_FOUND)
	FetchContent_Declare(
		Eigen3
		QUIET
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_TAG 3.4.0
	)
	FetchContent_MakeAvailable(Eigen3)
endif()

# Create the library
file(GLOB MCLBVH_INCLUDE
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/*.hpp)
file(GLOB MCLBVH_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/ccd_internal/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/ccd_internal/*.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/ccd_internal/*.cpp)
add_library(${PROJECT_NAME} ${MCLBVH_INCLUDE} ${MCLBVH_SRC})
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen TBB::tbb)
target_include_directories(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
add_library(mcl::ccd ALIAS ${PROJECT_NAME})

# build examples
if(MCLCCD_BUILD_EXAMPLES)

    # Libigl for rendering/interface
    option(LIBIGL_WITH_OPENGL "Use OpenGL" ON)
    option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW" ON)
    include(libigl)

    set(MCLCCD_DEP_LIBRARIES
        igl::core
        igl::opengl_glfw
        mcl::ccd)

	add_definitions(-DMCLCCD_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

	add_executable(dillo ${CMAKE_CURRENT_SOURCE_DIR}/examples/dillo.cpp)
	add_executable(traverser ${CMAKE_CURRENT_SOURCE_DIR}/examples/traverser.cpp)
	target_link_libraries(dillo ${MCLCCD_DEP_LIBRARIES})
	target_link_libraries(traverser ${MCLCCD_DEP_LIBRARIES})

endif()

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(FILES ${MCLBVH_INCLUDE} DESTINATION include/MCL)