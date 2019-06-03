# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(odgi)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

# Use openmp for parallelism, but it's configured differently on OSX
find_package(OpenMP)
if (OPENMP_FOUND)
  if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # assumes clang build
    # we can't reliably detect when we're using clang, so for the time being we assume
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Xpreprocessor -fopenmp")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
  endif()
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
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
set(sdsl-lite_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/include")
set(sdsl-lite-divsufsort_INCLUDE "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
set(sdsl-lite_LIB "${INSTALL_DIR}/src/sdsl-lite-build/lib")
set(sdsl-lite-divsufsort_LIB "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib")

# DYNAMIC (full build using its cmake config)
ExternalProject_Add(dynamic
  GIT_REPOSITORY "https://github.com/vgteam/DYNAMIC.git"
  GIT_TAG "625da93483a50cc9888270325017df4e7be9a36b"
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
  UPDATE_COMMAND ""
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
  GIT_TAG "2f7300c96399a76cad23861bedbfcfb00bed3e13"
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(handlegraph INSTALL_DIR)
set(handlegraph_INCLUDE "${INSTALL_DIR}/src/handlegraph/src/include")
set(handlegraph_LIB "${INSTALL_DIR}/src/handlegraph-build")

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  GIT_REPOSITORY "https://github.com/Taywee/args.git"
  GIT_TAG "3de44ec671db452cc0c4ef86399b108939768abb"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# parallel bsort (header-only, so here we just check if we can build it)
ExternalProject_Add(bsort
  GIT_REPOSITORY "https://github.com/ekg/bsort.git"
  GIT_TAG "c3ab0d3308424030e0a000645a26d2c10a59a124"
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
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(backwardscpp SOURCE_DIR)
set(backwardscpp_INCLUDE "${SOURCE_DIR}")

# sparsepp
ExternalProject_Add(sparsepp
  GIT_REPOSITORY "https://github.com/edawson/sparsepp.git"
  GIT_TAG "1c5a285e87b2244383efe0398f9992acd61234e9"
  BUILD_IN_SOURCE TRUE
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  UPDATE_COMMAND ""
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
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ska SOURCE_DIR)
set(ska_INCLUDE "${SOURCE_DIR}")

# intervaltree
ExternalProject_Add(intervaltree
  GIT_REPOSITORY "https://github.com/ekg/intervaltree.git"
  GIT_TAG "e8082c74a6f5c18de99d8b4cc4a55e2e62a1150d"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(intervaltree SOURCE_DIR)
set(intervaltree_INCLUDE "${SOURCE_DIR}")

# lodepng
ExternalProject_Add(lodepng
  GIT_REPOSITORY "https://github.com/ekg/lodepng.git"
  GIT_TAG "6b7f7b2b1c0eaab7042b0137fa19129618d1bb9f"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(lodepng SOURCE_DIR)
set(lodepng_INCLUDE "${SOURCE_DIR}")
set(lodepng_LIB "${SOURCE_DIR}/lib")

# structures
ExternalProject_Add(structures
  GIT_REPOSITORY "https://github.com/jeizenga/structures.git"
  GIT_TAG "71722385665aa724b466889bbec1afbfd375201f"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(structures SOURCE_DIR)
set(structures_INCLUDE "${SOURCE_DIR}/src/include")

ExternalProject_Add(sonlib
  GIT_REPOSITORY "https://github.com/vgteam/pinchSonLib.git"
  GIT_TAG "d25257e7ba615081078c41ac82dd70fffdc8baad"
  BUILD_IN_SOURCE TRUE
  INSTALL_COMMAND ""
  UPDATE_COMMAND ""
  BUILD_COMMAND "make"
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(sonlib SOURCE_DIR)
set(sonLib_inc_INCLUDE "${SOURCE_DIR}/sonLib/C/inc")
set(sonLib_lib_INCLUDE "${SOURCE_DIR}/sonLib/lib")
set(sonLib_LIB "${SOURCE_DIR}/sonLib/lib")

ExternalProject_Add(picosha256
  GIT_REPOSITORY "https://github.com/okdshin/PicoSHA2.git"
  GIT_TAG "b699e6c900be6e00152db5a3d123c1db42ea13d0"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(picosha256 SOURCE_DIR)
set(picosha256_INCLUDE "${SOURCE_DIR}")

set(CMAKE_BUILD_TYPE Release)

# set up our target executable and specify its dependencies and includes
add_executable(odgi
  ${CMAKE_SOURCE_DIR}/src/graph.cpp
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  ${CMAKE_SOURCE_DIR}/src/crash.cpp
  ${CMAKE_SOURCE_DIR}/src/position.cpp
  ${CMAKE_SOURCE_DIR}/src/threads.cpp
  ${CMAKE_SOURCE_DIR}/src/split.cpp
  ${CMAKE_SOURCE_DIR}/src/node.cpp
  ${CMAKE_SOURCE_DIR}/src/subgraph.cpp
  ${CMAKE_SOURCE_DIR}/src/cactus.cpp
  #${CMAKE_SOURCE_DIR}/src/snarls.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/driver.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/handle.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/fuzz.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/simplify.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/subcommand.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/build_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/test_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/stats_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/sort_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/view_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/kmers_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/viz_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/paths_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/cactus_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/prune_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/simplify_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/subset_main.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/topological_sort.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/kmer.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/hash.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/is_single_stranded.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/remove_high_degree.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/prune.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/coverage.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/cycle_breaking_sort.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/eades_algorithm.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/dagify.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/strongly_connected_components.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/weakly_connected_components.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/dfs.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/find_shortest_paths.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/id_ordered_paths.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/simple_components.cpp
  )
add_dependencies(odgi sdsl-lite)
add_dependencies(odgi dynamic)
add_dependencies(odgi gfakluge)
add_dependencies(odgi handlegraph)
add_dependencies(odgi tayweeargs)
add_dependencies(odgi bbhash)
add_dependencies(odgi sparsepp)
add_dependencies(odgi ska)
add_dependencies(odgi bsort)
add_dependencies(odgi lodepng)
add_dependencies(odgi structures)
add_dependencies(odgi sonlib)
add_dependencies(odgi picosha256)
target_include_directories(odgi PUBLIC
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
  "${intervaltree_INCLUDE}"
  "${lodepng_INCLUDE}"
  "${bbhash_INCLUDE}"
  "${bsort_INCLUDE}"
  "${structures_INCLUDE}"
  "${sonLib_lib_INCLUDE}"
  "${sonLib_inc_INCLUDE}"
  "${picosha256_INCLUDE}")
target_link_libraries(odgi
  "${sdsl-lite_LIB}/libsdsl.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort.a"
  "${sdsl-lite-divsufsort_LIB}/libdivsufsort64.a"
  "${handlegraph_LIB}/libhandlegraph.a"
  "${lodepng_LIB}/liblodepng.a"
  "${sonLib_LIB}/stPinchesAndCacti.a"
  "${sonLib_LIB}/3EdgeConnected.a"
  "${sonLib_LIB}/sonLib.a"
  "-ldl"
#  "-ltcmalloc_minimal"
  )

#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltcmalloc")

if (APPLE)
elseif (TRUE)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()
