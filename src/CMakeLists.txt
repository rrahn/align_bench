# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /src/CMakeLists.txt
#
# CMakeLists.txt file for align_bench.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (align_bench_app CXX)
message (STATUS "Configuring src/align_bench_app")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# if (MSVC)
#     message (STATUS "Yara is currently not build on Microsoft Visual because it requires OpenMP 3.0.")
#     return ()
# endif ()

# Update the cmake flags required for SeqAn.

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
# Search SeqAn and select dependencies.
# find_package (ZLIB)
# find_package (BZip2)

message (STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

set (TBB_INCLUDE_DIR "/opt/local/include/tbb/")
set (TBB_LIBRARY "/opt/local/lib/")

find_package (TBB)
find_package (OpenMP)
find_package (SeqAn REQUIRED)

# Warn if TBB wasn't find.
if (NOT TBB_FOUND)
    message (STATUS "   TBB not found")
else(NOT TBB_FOUND)
    message (STATUS "   Using TBB Version ${TBB_VERSION}")
    list (APPEND SEQAN_DEFINITIONS "${TBB_DEFINITIONS};-DSEQAN_TBB")
    list (APPEND SEQAN_LIBRARIES "${TBB_LIBRARIES}")
    list (APPEND SEQAN_INCLUDE_DIRS "${TBB_INCLUDE_DIRS}")
endif (NOT TBB_FOUND)

# Warn if OpenMP was not found.
if (NOT OPENMP_FOUND)
  message (STATUS "  No OpenMP support for parallel alignment computation.")
endif (NOT OPENMP_FOUND)

# Warn if Zlib was not found.
# if (NOT ZLIB_FOUND)
#   message (STATUS "  Zlib not found: building yara without gzipped input and BAM output.")
# endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
# if (NOT BZIP2_FOUND)
#   message (STATUS "  BZip2 not found: building yara without bzipped input.")
# endif (NOT BZIP2_FOUND)

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

set (SEQAN_APP_VERSION "0.0.1")

# option (YARA_LARGE_CONTIGS "Set to OFF to disable support for more than 32k contigs or contigs longer than 4Gbp." ON)
# if (YARA_LARGE_CONTIGS AND NOT SEQAN_TRAVIS_BUILD)
#   add_definitions (-DYARA_LARGE_CONTIGS=1)
# endif (YARA_LARGE_CONTIGS AND NOT SEQAN_TRAVIS_BUILD)

# ----------------------------------------------------------------------------
# Disable Compiler Warnings
# ----------------------------------------------------------------------------

# Disables warning about conversion from 64 bit to 32 bit integrals.
# Note that there is a runtime exception thrown if the type is not sufficient.
# Unfortunately the VC does not detect it.
if (MSVC)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4244")  # When compiling on 32 bit windows.
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4267")  # When compiling on 64 bit windows.
endif (MSVC)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable(align_bench align_bench.cpp
                           align_bench_options.hpp
                           sequence_generator.hpp
                           benchmark_executor.hpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (align_bench ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation - requires SeqAn build system.
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
# if (NOT SEQAN_PREFIX_SHARE_DOC)
#   seqan_setup_install_vars (align_bench)
# endif (NOT SEQAN_PREFIX_SHARE_DOC)

# # Install yara in ${PREFIX}/bin directory
# install (TARGETS align_bench
#          DESTINATION bin)

# # Install non-binary files for the package to "." for app builds and
# # ${PREFIX}/share/doc/yara for SeqAn release builds.
# install (FILES LICENSE README.rst
#          DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (yara)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable yara_mapper in CTD structure.
# set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES}
#                           align_bench CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

# if (SEQAN_BUILD_SYSTEM STREQUAL "APP:align_bench")
#   set (CPACK_PACKAGE_NAME "align_bench")
#   set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Align_Bench - Benchmark tool for alignments.")
#   set (CPACK_DEBIAN_PACKAGE_MAINTAINER "René Rahn <rene.rahn@fu-berlin.de>")
#   set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
#   set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.rst")

#   seqan_configure_cpack_app (align_bench "align_bench")
# endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:align_bench")
