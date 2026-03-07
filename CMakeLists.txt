# Copyright Matt Overby 2021.
# Distributed under the MIT License.

cmake_minimum_required(VERSION 3.0)
project(mclccd LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(MCLCCD_DEBUG "Compile in test mode" OFF)
option(MCLCCD_EXAMPLES "Compile visual examples" OFF)
option(MCLCCD_UNIT_TESTS "Compile unit tests" OFF)

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
find_package(TBB REQUIRED)

# Eigen linalg
find_package(Eigen3 REQUIRED)

# Create the library
set(MCLBVH_INCLUDE
	${CMAKE_CURRENT_SOURCE_DIR}/src/MCL/BVHLeaf.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MCL/BVHTree.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MCL/KdBVH.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MCL/NarrowPhase.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MCL/Projection.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/third-party/CTCD.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/third-party/Distance.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/third-party/rpoly.h
	${CMAKE_CURRENT_SOURCE_DIR}/third-party/tt_isect.hpp)
set(MCLBVH_SRC
	${CMAKE_CURRENT_SOURCE_DIR}/src/MCL/BVHLeaf.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MCL/BVHTree.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/MCL/NarrowPhase.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/third-party/CTCD.cpp)
add_library(${PROJECT_NAME} ${MCLBVH_INCLUDE} ${MCLBVH_SRC})
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen TBB::tbb)
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
add_library(mcl::ccd ALIAS ${PROJECT_NAME})

# build viewer tests
if(MCLCCD_EXAMPLES)
	include(libigl)
	igl_include(glfw)
	add_definitions(-DMCLCCD_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
	add_executable(dillo ${CMAKE_CURRENT_SOURCE_DIR}/test/dillo.cpp)
	add_executable(traverser ${CMAKE_CURRENT_SOURCE_DIR}/test/traverser.cpp)
	target_link_libraries(dillo igl::glfw mcl::ccd)
	target_link_libraries(traverser igl::glfw mcl::ccd)
endif()

# unit testing
if (MCLCCD_UNIT_TESTS)
	enable_testing()
	add_executable(test_2D ${CMAKE_CURRENT_SOURCE_DIR}/test/test_2D.cpp)
	target_link_libraries(test_2D PUBLIC mcl::ccd)
	add_test(NAME Test2D COMMAND test_2D)
endif()

# install
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(FILES ${MCLBVH_INCLUDE} DESTINATION include/MCL)