# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(dankgraph)
# We build using c++11
set(CMAKE_CXX_STANDARD 11)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O2")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O2")

# We use OpenMP for parallelism
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
  GIT_TAG "d52aa9a71513d132e30c09491b5899af449ebb94"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# DYNAMIC (full build using its cmake config)
ExternalProject_Add(dynamic
  GIT_REPOSITORY "https://github.com/xxsds/DYNAMIC.git"
  GIT_TAG "e46c0245bafb9d5b157c72be7dc71142c8ae9239"
  # we don't actually install dynamic... it's header only
  #CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(dynamic INSTALL_DIR)
set(dynamic_INCLUDE "${INSTALL_DIR}/src/dynamic/include")

# gfakluge (full build using its cmake config)
ExternalProject_Add(gfakluge
  GIT_REPOSITORY "https://github.com/edawson/gfakluge.git"
  GIT_TAG "8c30259b86fcfd592477e5f21e3999ddac850fdd"
  BUILD_IN_SOURCE TRUE
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(gfakluge INSTALL_DIR)
set(gfakluge_INCLUDE "${INSTALL_DIR}/src/gfakluge")
set(gfakluge_LIB "${INSTALL_DIR}/src/gfakluge")


set(CMAKE_BUILD_TYPE Debug) 

# set up our target executable and specify its dependencies and includes
add_executable(dankgraph
  ${CMAKE_SOURCE_DIR}/src/graph.cpp
  ${CMAKE_SOURCE_DIR}/src/dynamic_structs.cpp
  ${CMAKE_SOURCE_DIR}/src/path.cpp
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  )
add_dependencies(dankgraph sdsl-lite)
add_dependencies(dankgraph dynamic)
target_include_directories(dankgraph PUBLIC
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${dynamic_INCLUDE}"
  )
target_link_libraries(dankgraph
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  )

if (APPLE)
elseif (TRUE)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()
