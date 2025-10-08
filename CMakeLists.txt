# Copyright Matt Overby 2021.
# Distributed under the MIT License.
cmake_minimum_required(VERSION 3.0)
project(mclccd LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FetchContent)

option(MCLCCD_DEBUG "Compile in test mode" OFF)
option(MCLCCD_TESTING_ONLY "only compile unit tests" OFF)

if (MCLCCD_DEBUG)
	set(CMAKE_BUILD_TYPE Debug)
	if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wundef -Wextra -Wno-long-long -Wlogical-op -Wshadow -Wmaybe-uninitialized")
	endif()
else()
	if (NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE "Release")
	endif()
endif()

# thread building blocks
# TODO: build this with fetchcontent if not found
find_package(TBB REQUIRED)

# Eigen, download if not found
find_package(Eigen3)
if(NOT EIGEN_FOUND)
	option(BUILD_TESTING "Enable creation of Eigen tests." OFF)
	FetchContent_Declare(
		Eigen3
		QUIET
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_TAG 3.4.0
		GIT_SHALLOW TRUE
	)
	FetchContent_MakeAvailable(Eigen3)
endif()

# Create the library
file(GLOB MCLBVH_INCLUDE
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/*.hpp)
file(GLOB MCLBVH_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/include/MCL/*.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third-party/*.h
	${CMAKE_CURRENT_SOURCE_DIR}/third-party/*.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/third-party/*.cpp)
add_library(${PROJECT_NAME} ${MCLBVH_INCLUDE} ${MCLBVH_SRC})
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen TBB::tbb)
target_include_directories(${PROJECT_NAME} INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third-party>)
add_library(mcl::ccd ALIAS ${PROJECT_NAME})

# build examples
if(NOT MCLCCD_TESTING_ONLY)

    # Libigl for rendering/interface
    include(libigl)
	igl_include(glfw)

	# Set define for loading objs
	add_definitions(-DMCLCCD_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

	add_executable(dillo ${CMAKE_CURRENT_SOURCE_DIR}/examples/dillo.cpp)
	add_executable(traverser ${CMAKE_CURRENT_SOURCE_DIR}/examples/traverser.cpp)
	target_link_libraries(dillo igl::glfw mcl::ccd)
	target_link_libraries(traverser igl::glfw mcl::ccd)

endif()

# unit testing
enable_testing()
add_executable(test_2D ${CMAKE_CURRENT_SOURCE_DIR}/test/test_2D.cpp)
target_link_libraries(test_2D PUBLIC mcl::ccd)
add_test(NAME Test2D COMMAND test_2D)

# install
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(FILES ${MCLBVH_INCLUDE} DESTINATION include/MCL)