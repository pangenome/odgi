# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(dankgraph)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

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

# TODO: We're using INSTALL_DIR very wrong. We *should* be actually installing
# the external projects into their prefixes and working with the installed
# files. Instead we're building but not installing them and trying to work with
# the non-installed build trees.
# 
# Hence the blanked out INSTALL_COMMANDs to suppress the install step.
#
# We need to NOT blank out UPDATE_COMMAND or we can never change the Git revision we point to.
# The cost of this is that we have to re-configure on every build if we do update.

# sdsl-lite (full build using its cmake config)
ExternalProject_Add(sdsl-lite
  GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
  GIT_TAG "d52aa9a71513d132e30c09491b5899af449ebb94"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# DYNAMIC (full build using its cmake config)
ExternalProject_Add(dynamic
  GIT_REPOSITORY "https://github.com/vgteam/DYNAMIC.git"
  GIT_TAG "615d8be5276bcd4c5a3d8e31679b4f8e81b2eefc"
  # we don't actually install dynamic... it's header only
  #CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(dynamic INSTALL_DIR)
set(dynamic_INCLUDE "${INSTALL_DIR}/src/dynamic/include")

# gfakluge (now header only)
ExternalProject_Add(gfakluge
  GIT_REPOSITORY "https://github.com/vgteam/gfakluge.git"
  GIT_TAG "d99d927ab286c3dfaeb4db15c15d306d1fab0806"
  BUILD_IN_SOURCE TRUE
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(gfakluge INSTALL_DIR)
set(gfakluge_INCLUDE "${INSTALL_DIR}/src/gfakluge/src")
set(gfakluge_tinyFA_INCLUDE "${INSTALL_DIR}/src/gfakluge/src/tinyFA")
set(gfakluge_LIB "${INSTALL_DIR}/src/gfakluge")

# libhandlegraph (full build using its cmake config)
ExternalProject_Add(handlegraph
  GIT_REPOSITORY "https://github.com/vgteam/libhandlegraph.git"
  GIT_TAG "4fb0a4b750995c35e45b8c510397339af4ddc2fc"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  INSTALL_COMMAND "")
ExternalProject_Get_property(handlegraph INSTALL_DIR)
set(handlegraph_INCLUDE "${INSTALL_DIR}/src/handlegraph/src/include")
set(handlegraph_LIB "${INSTALL_DIR}/src/handlegraph-build")

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  GIT_REPOSITORY "https://github.com/Taywee/args.git"
  GIT_TAG "3de44ec671db452cc0c4ef86399b108939768abb"
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# parallel bsort (header-only, so here we just check if we can build it)
ExternalProject_Add(bsort
  GIT_REPOSITORY "https://github.com/ekg/bsort.git"
  GIT_TAG "b14c3302fa6379cc5d49fdd023d375032e801c41"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(bsort SOURCE_DIR)
set(bsort_INCLUDE "${SOURCE_DIR}/src")
set(bsort_LIB "${SOURCE_DIR}/lib")

# BBHash perfect hasher
ExternalProject_Add(bbhash
  GIT_REPOSITORY "https://github.com/vgteam/BBHash.git"
  GIT_TAG "ccc22805eae48aca2a01a7ff5ef8e002857020d3"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(bbhash SOURCE_DIR)
set(bbhash_INCLUDE "${SOURCE_DIR}")

# backwardscpp
ExternalProject_Add(backwardscpp
  GIT_REPOSITORY "https://github.com/bombela/backward-cpp.git"
  GIT_TAG "b353585d4bd19f9e51816bc11d2e53b0e5b4a760"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  INSTALL_COMMAND "")
ExternalProject_Get_property(backwardscpp SOURCE_DIR)
set(backwardscpp_INCLUDE "${SOURCE_DIR}")

# sparsepp
ExternalProject_Add(sparsepp
  GIT_REPOSITORY "https://github.com/edawson/sparsepp.git"
  GIT_TAG "1c5a285e87b2244383efe0398f9992acd61234e9"
  BUILD_IN_SOURCE TRUE
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(sparsepp INSTALL_DIR)
set(sparsepp_INCLUDE "${INSTALL_DIR}/src/sparsepp/sparsepp/")
set(sparsepp_LIB "${INSTALL_DIR}/src/sparsepp/sparsepp/")

# ska
ExternalProject_Add(ska
  GIT_REPOSITORY "https://github.com/skarupke/flat_hash_map.git"
  GIT_TAG "2c4687431f978f02a3780e24b8b701d22aa32d9c"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ska SOURCE_DIR)
set(ska_INCLUDE "${SOURCE_DIR}")

set(CMAKE_BUILD_TYPE Release)

# set up our target executable and specify its dependencies and includes
add_executable(dsgvg
  ${CMAKE_SOURCE_DIR}/src/graph.cpp
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  ${CMAKE_SOURCE_DIR}/src/crash.cpp
  ${CMAKE_SOURCE_DIR}/src/position.cpp
  ${CMAKE_SOURCE_DIR}/src/threads.cpp
  ${CMAKE_SOURCE_DIR}/src/node.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/driver.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/handle.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/subcommand.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/build_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/test_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/stats_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/sort_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/view_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/kmers_main.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/topological_sort.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/kmer.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/is_single_stranded.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/remove_high_degree.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/prune.cpp
  )
add_dependencies(dsgvg sdsl-lite)
add_dependencies(dsgvg dynamic)
add_dependencies(dsgvg gfakluge)
add_dependencies(dsgvg handlegraph)
add_dependencies(dsgvg tayweeargs)
add_dependencies(dsgvg bbhash)
add_dependencies(dsgvg sparsepp)
add_dependencies(dsgvg ska)
add_dependencies(dsgvg bsort)
target_include_directories(dsgvg PUBLIC
  "${CMAKE_SOURCE_DIR}/src"
  "${sdsl-lite_INCLUDE}"
  "${sdsl-lite-divsufsort_INCLUDE}"
  "${dynamic_INCLUDE}"
  "${gfakluge_INCLUDE}"
  "${gfakluge_tinyFA_INCLUDE}"
  "${handlegraph_INCLUDE}"
  "${tayweeargs_INCLUDE}"
  "${backwardscpp_INCLUDE}"
  "${sparsepp_INCLUDE}"
  "${ska_INCLUDE}"
  "${bbhash_INCLUDE}"
  "${bsort_INCLUDE}")
target_link_libraries(dsgvg
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "${handlegraph_LIB}/libhandlegraph.a"
  "-ldl"
#  "-ltcmalloc_minimal"
  )

#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltcmalloc")

if (APPLE)
elseif (TRUE)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()
