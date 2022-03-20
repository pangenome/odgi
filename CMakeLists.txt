# Specify the minimum version for CMake
#
# Accepted cmake flags are:
#
#   (default builds all with -fPIC and -g, not optimal!)
#
# For best performance with gcc>=11 compile the binary with
#
#     cmake -DPIC=OFF -DEXTRA_FLAGS="-Ofast -march=native -pipe -msse4.2 -funroll-all-loops" ..
# To build the smallest possible executable use (from 8.6Mb to 5.5Mb at a cost of 25% speed)
#
#     cmake -DPIC=OFF -DEXTRA_FLAGS="-Os -ffunction-sections -fdata-sections -march=native -pipe -msse4.2 -Wl,--gc-sections" ..
#
# These architecture specific flags can not be the default.
#
# For the performance of the odgi binary it pays to compile -DPIC=OFF and run performance guided optimization (PGO).
#
# For more information see ./INSTALL.md

cmake_minimum_required(VERSION 3.1)

# Project's name
project(odgi)
# Enforce c++17
set(CMAKE_CXX_STANDARD 17)

# Command line switches. Compile with cmake -DINLINE_HANDLEGRAPH_SOURCES=ON
option(PIC "Compile all odgi sources with -fPIC - required for shared libs" ON)
option(INLINE_HANDLEGRAPH_SOURCES "Compile handlegraph sources inline" OFF)

include(ExternalProject)
include(FeatureSummary)

find_package(PkgConfig REQUIRED)
find_package(pybind11 CONFIG)
find_package(OpenMP)

feature_summary(
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  WHAT REQUIRED_PACKAGES_NOT_FOUND)

pkg_check_modules(SDSLLITE sdsl-lite)
pkg_check_modules(LIBDIVSUFSORT libdivsufsort)

if (NOT CMAKE_BUILD_TYPE) # default build to Release
  set(CMAKE_BUILD_TYPE Release)
endif()

if (PIC)
  message(STATUS "Compiling odgi sources with PIC=ON")
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif (PIC)

# Set optimization through command line; see INSTALL.md
# set(EXTRA_FLAGS "-Ofast -march=native -pipe -msse4.2 -funroll-all-loops") #  -fprofile-generate=../pgo")

if (${CMAKE_BUILD_TYPE} MATCHES "Release")
  if (EXTRA_FLAGS)
    set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG") # reset CXX_FLAGS to replace -O3 with -Ofast
  endif (EXTRA_FLAGS)
  # Use all standard-compliant optimizations - always add these:
  set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
  set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
endif ()

# set(CMAKE_BUILD_TYPE Debug) -- don't uncomment this, instead run
#   cmake -DCMAKE_BUILD_TYPE=Debug ..

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  # Debug use the defaults - so commenting out:
  # set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -fsanitize=address")
  # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -fsanitize=address")
endif ()


message(STATUS "ODGI CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # assumes clang build
  # we can't reliably detect when we're using clang, so for the time being we assume
  # TODO: can't we though?
  # adapted from https://stackoverflow.com/questions/46414660/macos-cmake-and-openmp
  # find_package(OpenMP) does not work reliably on macOS, so we do its work ourselves
  set (OpenMP_C "${CMAKE_C_COMPILER}")
  set (OpenMP_C_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
  set (OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set (OpenMP_CXX "${CMAKE_CXX_COMPILER}")
  set (OpenMP_CXX_FLAGS " -Xpreprocessor -fopenmp -I/opt/local/include/libomp -I/usr/local/include -L/opt/local/lib/libomp -L/usr/local/lib")
  set (OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
  set (OpenMP_libomp_LIBRARY "omp")
  set (OpenMP_libgomp_LIBRARY "gomp")
  set (OpenMP_libiomp5_LIBRARY "iomp5")
  # and now add the OpenMP parameters to the compile flags
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -lomp")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  find_package(OpenMP REQUIRED)
  # add the flags it detects to the compile flags
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
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

if (NOT SDSLLITE_FOUND)
  # sdsl-lite (full build using its cmake config)
  ExternalProject_Add(sdsl-lite
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sdsl-lite"
    CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
    UPDATE_COMMAND ""
    INSTALL_COMMAND "")
  ExternalProject_Get_property(sdsl-lite INSTALL_DIR)
  set(SDSLLITE_INCLUDE_DIRS "${INSTALL_DIR}/src/sdsl-lite-build/include")
  set(LIBDIVSUFSORT_INCLUDE_DIRS "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/include")
  set(SDSLLITE_LINK_LIBRARIES "${INSTALL_DIR}/src/sdsl-lite-build/lib/libsdsl.a")
  set(LIBDIVSUFSORT_LINK_LIBRARIES "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib/libdivsufsort.a" "${INSTALL_DIR}/src/sdsl-lite-build/external/libdivsufsort/lib/libdivsufsort64.a")
endif (NOT SDSLLITE_FOUND)

# DYNAMIC (full build using its cmake config)
ExternalProject_Add(dynamic
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/DYNAMIC"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(dynamic SOURCE_DIR)
set(dynamic_INCLUDE "${SOURCE_DIR}/include")

# hopscotch_map (required by DYNAMIC)
ExternalProject_Add(hopscotch_map
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/hopscotch-map"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(hopscotch_map SOURCE_DIR)
set(hopscotch_map_INCLUDE "${SOURCE_DIR}/include")

# gfakluge (now header only)
ExternalProject_Add(gfakluge
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/gfakluge"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(gfakluge SOURCE_DIR)
set(gfakluge_INCLUDE "${SOURCE_DIR}/src")
set(gfakluge_tinyFA_INCLUDE "${SOURCE_DIR}/src/tinyFA")
ExternalProject_Get_property(gfakluge INSTALL_DIR)
set(gfakluge_LIB "${INSTALL_DIR}/src/gfakluge")

if (NOT INLINE_HANDLEGRAPH_SOURCES)
  # libhandlegraph (full build using its cmake config)
  ExternalProject_Add(handlegraph
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/libhandlegraph"
    CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
  ExternalProject_Get_property(handlegraph INSTALL_DIR)
  set(handlegraph_INCLUDE "${INSTALL_DIR}/include")
  set(handlegraph_LIB "${INSTALL_DIR}/lib")

else (NOT INLINE_HANDLEGRAPH_SOURCES)
  MESSAGE(STATUS "ODGI inlining handlegraph sources")
  set(handlegraph_INCLUDE "${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/include")

  set(handlegraph_sources
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/append_graph.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/apply_orientations.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/are_equivalent.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/copy_graph.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/count_walks.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/dagify.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/deletable_handle_graph.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/dijkstra.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/eades_algorithm.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/extend.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/find_shortest_paths.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/find_tips.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/handle_graph.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/is_acyclic.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/is_single_stranded.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/mutable_handle_graph.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/path_handle_graph.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/path_position_handle_graph.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/ranked_handle_graph.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/reverse_complement.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/serializable.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/snarl_decomposition.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/split_strands.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/strongly_connected_components.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/topological_sort.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/trivially_serializable.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/types.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/unchop.cpp
  ${CMAKE_SOURCE_DIR}/deps/libhandlegraph/src/weakly_connected_components.cpp
)
endif (NOT INLINE_HANDLEGRAPH_SOURCES)

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/args"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

# BBHash perfect hasher
ExternalProject_Add(bbhash
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/BBHash"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(bbhash SOURCE_DIR)
set(bbhash_INCLUDE "${SOURCE_DIR}")

# sparsepp
ExternalProject_Add(sparsepp
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sparsepp"
  BUILD_IN_SOURCE TRUE
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(sparsepp SOURCE_DIR)
set(sparsepp_INCLUDE "${SOURCE_DIR}/sparsepp/")
ExternalProject_Get_property(sparsepp INSTALL_DIR)
set(sparsepp_LIB "${INSTALL_DIR}/src/sparsepp/sparsepp/")

# ska
ExternalProject_Add(ska
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/flat_hash_map"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ska SOURCE_DIR)
set(ska_INCLUDE "${SOURCE_DIR}")

# intervaltree
ExternalProject_Add(intervaltree
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/intervaltree"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(intervaltree SOURCE_DIR)
set(intervaltree_INCLUDE "${SOURCE_DIR}")

# cgranges
ExternalProject_Add(cgranges
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/cgranges"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(cgranges SOURCE_DIR)
set(cgranges_INCLUDE "${SOURCE_DIR}/cpp")

# mmmulti (memory mapped multimap, multiset, and interval tree)
ExternalProject_Add(mmmulti
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmmulti"
  BUILD_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(mmmulti SOURCE_DIR)
set(mmmulti_INCLUDE "${SOURCE_DIR}/src")

# In-place Parallel Super Scalar Samplesort (IPS⁴o), header only
ExternalProject_Add(ips4o
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/ips4o"
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(ips4o SOURCE_DIR)
set(ips4o_INCLUDE "${SOURCE_DIR}")

# atomic queue library
ExternalProject_Add(atomicqueue
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/atomic_queue/include/atomic_queue"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(atomicqueue SOURCE_DIR)
set(atomicqueue_INCLUDE "${SOURCE_DIR}/")

# lodepng
ExternalProject_Add(lodepng
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/lodepng"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(lodepng SOURCE_DIR)
set(lodepng_INCLUDE "${SOURCE_DIR}")
set(lodepng_LIB "${SOURCE_DIR}/lib")

# structures
ExternalProject_Add(structures
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/structures"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(structures SOURCE_DIR)
set(structures_INCLUDE "${SOURCE_DIR}/src/include")

ExternalProject_Add(picosha256
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/PicoSHA2"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(picosha256 SOURCE_DIR)
set(picosha256_INCLUDE "${SOURCE_DIR}")

# SGD based graph layout
ExternalProject_Add(sgd2
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/sgd2"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(sgd2 SOURCE_DIR)
set(sgd2_INCLUDE "${SOURCE_DIR}/src")

#ExternalProject_Add(mondriaan
#  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mondriaan"
#  UPDATE_COMMAND ""
#  INSTALL_COMMAND "")
#ExternalProject_Get_property(mondriaan SOURCE_DIR)
#set(mondriaan_INCLUDE "${SOURCE_DIR}/src")
#set(mondriaan_LIB "${SOURCE_DIR}/lib")

ExternalProject_Add(libbf
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/libbf"
  CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
ExternalProject_Get_property(libbf INSTALL_DIR)
#message(STATUS "libbf target " ${INSTALL_DIR})
set(libbf_INCLUDE "${INSTALL_DIR}/include")
set(libbf_LIB "${INSTALL_DIR}/lib")

# httplib for HTTP server
ExternalProject_Add(httplib
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/cpp-httplib"
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        BUILD_COMMAND ""
        CONFIGURE_COMMAND "")
ExternalProject_Get_property(httplib SOURCE_DIR)
set(httplib_INCLUDE "${SOURCE_DIR}")

# cpp_random_distributions for Zipfian distribution
ExternalProject_Add(random_distributions
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/cpp_random_distributions"
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        BUILD_COMMAND ""
        CONFIGURE_COMMAND "")
ExternalProject_Get_property(random_distributions SOURCE_DIR)
set(random_distributions_INCLUDE "${SOURCE_DIR}")

# dirtyzipf pow-approximate Zipf distribution
ExternalProject_Add(dirtyzipf
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/dirtyzipf"
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        BUILD_COMMAND ""
        CONFIGURE_COMMAND "")
ExternalProject_Get_property(dirtyzipf SOURCE_DIR)
set(dirtyzipf_INCLUDE "${SOURCE_DIR}")

# header-only pseudorandom number generator library
ExternalProject_Add(xoshiro
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/Xoshiro-cpp"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(xoshiro SOURCE_DIR)
set(xoshiro_INCLUDE "${SOURCE_DIR}")

ExternalProject_Add(atomicbitvector
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/atomicbitvector/include"
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        BUILD_COMMAND ""
        CONFIGURE_COMMAND "")
ExternalProject_Get_property(atomicbitvector SOURCE_DIR)
set(atomicbitvector_INCLUDE "${SOURCE_DIR}")

#add_subdirectory(deps/mmmulti/deps/mio)
ExternalProject_Add(mio
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/deps/mmmulti/deps/mio"
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
        BUILD_COMMAND ""
        CONFIGURE_COMMAND "")
ExternalProject_Get_property(mio SOURCE_DIR)
set(mio_INCLUDE "${SOURCE_DIR}/include")

# set up our target executable and specify its dependencies and includes
add_library(odgi_objs OBJECT
  ${CMAKE_SOURCE_DIR}/src/odgi.cpp
  ${CMAKE_SOURCE_DIR}/src/odgi-api.cpp
  ${CMAKE_SOURCE_DIR}/src/reclaimer.cpp
  ${CMAKE_SOURCE_DIR}/src/utils.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/subgraph/region.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/subgraph/extract.cpp
  ${CMAKE_SOURCE_DIR}/src/position.cpp
  ${CMAKE_SOURCE_DIR}/src/gfa_to_handle.cpp
  ${CMAKE_SOURCE_DIR}/src/split.cpp
  ${CMAKE_SOURCE_DIR}/src/node.cpp
  ${CMAKE_SOURCE_DIR}/src/subgraph.cpp
  ${CMAKE_SOURCE_DIR}/src/version.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/depth_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/overlap_main.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/driver.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/handle.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/fuzz.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/simplify.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/sort.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/pathindex.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/edge.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/extract.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/stepindex.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/subcommand.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/build_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/test_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/stats_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/cover_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/explode_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/squeeze_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/sort_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/view_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/kmers_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/unitig_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/viz_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/paths_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/prune_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/unchop_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/normalize_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/extract_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/position_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/degree_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/bin_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/matrix_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/chop_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/crush_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/groom_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/layout0_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/layout_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/draw_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/flatten_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/break_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/pathindex_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/panpos_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/server_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/version_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/untangle_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/tips_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/stepindex_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/heaps_main.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/pav_main.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/topological_sort.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/kmer.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/hash.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/is_single_stranded.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/remove_high_degree.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/prune.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/depth.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/degree.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/cycle_breaking_sort.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/random_order.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/eades_algorithm.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/dagify.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/dagify_sort.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/split_strands.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/strongly_connected_components.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/weakly_connected_components.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/dfs.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/bfs.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/find_shortest_paths.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/id_ordered_paths.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/simple_components.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/bin_path_info.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/bin_path_depth.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/sgd_layout.cpp
#  ${CMAKE_SOURCE_DIR}/src/algorithms/mondriaan_sort.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/matrix_writer.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/temp_file.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/linear_index.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/linear_sgd.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/break_cycles.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/xp.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/cut_tips.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/merge.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/normalize.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/simplify_siblings.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/chop.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/unchop.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/perfect_neighbors.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/cover.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/path_sgd.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/path_sgd_layout.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/draw.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/layout.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/atomic_image.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/remove_isolated.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/expand_context.cpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/validate_main.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/untangle.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/stepindex.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/groom.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/crush_n.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/heaps.cpp
  ${CMAKE_SOURCE_DIR}/src/unittest/edge.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/tips.cpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/path_jaccard.cpp
  ${handlegraph_sources}
)

set(odgi_DEPS
    dynamic
    hopscotch_map
    gfakluge
    tayweeargs
    bbhash
    sparsepp
    ska
    intervaltree
    cgranges
    lodepng
    structures
    picosha256
    sgd2
    libbf
    httplib
    random_distributions
    dirtyzipf
    mmmulti
    ips4o
    xoshiro
    atomicqueue
    mio)

if (NOT SDSLLITE_FOUND)
    list(APPEND odgi_DEPS sdsl-lite)
endif (NOT SDSLLITE_FOUND)

if (NOT INLINE_HANDLEGRAPH_SOURCES)
    list(APPEND odgi_DEPS handlegraph)
endif (NOT INLINE_HANDLEGRAPH_SOURCES)

add_dependencies(odgi_objs ${odgi_DEPS})

set(odgi_INCLUDES
  "${CMAKE_SOURCE_DIR}/src"
  "${CMAKE_SOURCE_DIR}/src/algorithms"
  "${SDSLLITE_INCLUDE_DIRS}"
  "${LIBDIVSUFSORT_INCLUDE_DIRS}"
  "${dynamic_INCLUDE}"
  "${hopscotch_map_INCLUDE}"
  "${gfakluge_INCLUDE}"
  "${gfakluge_tinyFA_INCLUDE}"
  "${handlegraph_INCLUDE}"
  "${tayweeargs_INCLUDE}"
  "${sparsepp_INCLUDE}"
  "${ska_INCLUDE}"
  "${intervaltree_INCLUDE}"
  "${cgranges_INCLUDE}"
  "${mmmulti_INCLUDE}"
  "${ips4o_INCLUDE}"
  "${atomicqueue_INCLUDE}"
  "${lodepng_INCLUDE}"
  "${bbhash_INCLUDE}"
  "${structures_INCLUDE}"
  "${picosha256_INCLUDE}"
  "${sgd2_INCLUDE}"
#  "${mondriaan_INCLUDE}"
  "${libbf_INCLUDE}"
  "${httplib_INCLUDE}"
  "${random_distributions_INCLUDE}"
  "${dirtyzipf_INCLUDE}"
  "${xoshiro_INCLUDE}"
  "${atomicbitvector_INCLUDE}"
  "${mio_INCLUDE}")

set(odgi_LIBS
  ${SDSLLITE_LINK_LIBRARIES}
  ${LIBDIVSUFSORT_LINK_LIBRARIES}
  "${lodepng_LIB}/liblodepng.a"
#  "${mondriaan_LIB}/libmondriaan.a"
  "${libbf_LIB}/libbf.a"
  "-ldl"
  "-latomic"
  jemalloc)
  #"-lefence") # for malloc error checking
  #"-ltcmalloc") # for heap profiling

if (NOT INLINE_HANDLEGRAPH_SOURCES)
    list(APPEND odgi_LIBS "${handlegraph_LIB}/libhandlegraph.a")
endif (NOT INLINE_HANDLEGRAPH_SOURCES)

set(odgi_HEADERS
  ${CMAKE_SOURCE_DIR}/include/odgi_git_version.hpp
  ${CMAKE_SOURCE_DIR}/src/hash_map.hpp
  ${CMAKE_SOURCE_DIR}/src/odgi.hpp
  ${CMAKE_SOURCE_DIR}/src/odgi-api.h
  ${CMAKE_SOURCE_DIR}/src/node.hpp
  ${CMAKE_SOURCE_DIR}/src/bmap.hpp
  ${CMAKE_SOURCE_DIR}/src/subgraph.hpp
  ${CMAKE_SOURCE_DIR}/src/split.hpp
  ${CMAKE_SOURCE_DIR}/src/varint.hpp
  ${CMAKE_SOURCE_DIR}/src/dna.hpp
  ${CMAKE_SOURCE_DIR}/src/phf.hpp
  ${CMAKE_SOURCE_DIR}/src/bgraph.hpp
  ${CMAKE_SOURCE_DIR}/src/gfa_to_handle.hpp
  ${CMAKE_SOURCE_DIR}/src/subcommand/subcommand.hpp
  ${CMAKE_SOURCE_DIR}/src/io_helper.hpp
  ${CMAKE_SOURCE_DIR}/src/version.hpp
  ${CMAKE_SOURCE_DIR}/src/btypes.hpp
  ${CMAKE_SOURCE_DIR}/src/dynamic_types.hpp
  ${CMAKE_SOURCE_DIR}/src/dynamic_structs.hpp
  ${CMAKE_SOURCE_DIR}/src/position.hpp
  ${CMAKE_SOURCE_DIR}/src/dset64.hpp
  ${CMAKE_SOURCE_DIR}/src/lockfree_hashtable.hpp
  ${CMAKE_SOURCE_DIR}/src/reclaimer.hpp
  ${CMAKE_SOURCE_DIR}/src/colorbrewer.hpp
  ${CMAKE_SOURCE_DIR}/src/unittest/driver.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/linear_index.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/random_order.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/cycle_breaking_sort.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/prune.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/reverse_complement.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/bin_path_info.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/bin_path_depth.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/dfs.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/chop.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/unchop.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/apply_bulk_modifications.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/dagify.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/weakly_connected_components.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/extract_containing_graph.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/eades_algorithm.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/perfect_neighbors.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/temp_file.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/distance_to_tail.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/simple_components.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/cut_tips.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/break_cycles.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/bfs.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/is_acyclic.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/id_sort.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/linear_sgd.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/matrix_writer.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/find_shortest_paths.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/extend.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/split_strands.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/extract_extending_graph.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/extract_connecting_graph.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/cover.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/path_sgd.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/path_sgd_layout.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/kmer.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/expand_context.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/id_ordered_paths.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/normalize.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/merge.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/count_walks.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/remove_high_degree.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/a_star.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/dagify_sort.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/xp.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/remove_isolated.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/sgd_term.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/simplify_siblings.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/sgd_layout.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/topological_sort.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/depth.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/degree.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/sorted_id_ranges.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/strongly_connected_components.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/hash.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/hilbert.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/groom.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/distance_to_head.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/is_single_stranded.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/shortest_cycle.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/untangle.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/progress.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/tips.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/tips_bed_writer_thread.hpp
  ${CMAKE_SOURCE_DIR}/src/algorithms/path_jaccard.hpp)

target_include_directories(odgi_objs PUBLIC ${odgi_INCLUDES})

add_library(libodgi_static STATIC $<TARGET_OBJECTS:odgi_objs>)
set_target_properties(libodgi_static PROPERTIES OUTPUT_NAME "odgi")
set_target_properties(libodgi_static PROPERTIES PUBLIC_HEADER "${odgi_HEADERS}")

if (NOT PIC)
  MESSAGE(STATUS "libodgi.so requires -DPIC=ON")
else (NOT PIC)
  add_library(libodgi_shared SHARED $<TARGET_OBJECTS:odgi_objs> src/algorithms/fonts/field8.h src/algorithms/fonts/font5x8.h)
  set_target_properties(libodgi_shared PROPERTIES OUTPUT_NAME "odgi")
  set_target_properties(libodgi_shared PROPERTIES PUBLIC_HEADER "${odgi_HEADERS}")
endif (NOT PIC)

add_executable(odgi
  $<TARGET_OBJECTS:odgi_objs>
  ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_link_libraries(odgi ${odgi_LIBS})
set_target_properties(odgi PROPERTIES OUTPUT_NAME "odgi")

if (NOT PIC)
  MESSAGE(STATUS "Can not build python bindings with PIC=OFF")
else (NOT PIC)
  if (NOT pybind11_FOUND)
    add_subdirectory(deps/pybind11)
  endif (NOT pybind11_FOUND)

  # Build Python modules, first the FFI
  pybind11_add_module(odgi_ffi "${CMAKE_SOURCE_DIR}/src/pythonffi.cpp")
  add_dependencies(odgi_ffi ${odgi_DEPS} libodgi_static)
  target_include_directories(odgi_ffi PUBLIC ${odgi_INCLUDES})
  target_link_libraries(odgi_ffi PRIVATE "${CMAKE_SOURCE_DIR}/lib/libodgi.a" "${odgi_LIBS}")
  set_target_properties(odgi_ffi PROPERTIES OUTPUT_NAME "odgi_ffi")
  install(TARGETS odgi_ffi LIBRARY DESTINATION lib)

  # Build original Python module
  pybind11_add_module(odgi_pybind11 "${CMAKE_SOURCE_DIR}/src/pythonmodule.cpp")
  add_dependencies(odgi_pybind11 ${odgi_DEPS} libodgi_static)
  target_include_directories(odgi_pybind11 PUBLIC ${odgi_INCLUDES})
  target_link_libraries(odgi_pybind11 PRIVATE "${CMAKE_SOURCE_DIR}/lib/libodgi.a" "${odgi_LIBS}")
  set_target_properties(odgi_pybind11 PROPERTIES OUTPUT_NAME "odgi")
  install(TARGETS odgi_pybind11 LIBRARY DESTINATION lib)
endif (NOT PIC)

file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)
execute_process(COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/generate_git_version.sh ${CMAKE_SOURCE_DIR}/include)

install(TARGETS odgi DESTINATION bin)
install(TARGETS libodgi_static ARCHIVE DESTINATION lib LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include/odgi)
if (PIC)
  install(TARGETS libodgi_shared ARCHIVE DESTINATION lib LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include/odgi)
endif (PIC)

if (APPLE)
elseif (TRUE)
  if (BUILD_STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
  endif()
endif()

# ---- Add tests
enable_testing()

SET(CTEST_ENVIRONMENT
"LD_LIBRARY_PATH=/gnu/store/6hcvbvicszj935bfnh99q6gll99cl22p-profile/lib"
)

function(add_pydoctest TEST_FILE)
  add_test(
      NAME ${TEST_FILE}
      COMMAND python3 -m doctest -o NORMALIZE_WHITESPACE -o REPORT_UDIFF python/${TEST_FILE}.md
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test)
  set_tests_properties(${TEST_FILE} PROPERTIES ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/lib;LD_LIBRARY_PATH=/gnu/store/6hcvbvicszj935bfnh99q6gll99cl22p-profile/lib")

endfunction()

add_pydoctest(odgi_ffi)
add_test(NAME odgi-test COMMAND odgi test)
