# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (align_bench CXX)

# ===========================================================================
# Setup Modules and Find Packages.
# ===========================================================================

# Note that the find modules in util/make take precedence
# over the built-in ones from CMake.  This is used in
# FindBZip2, for example which we slightly modify.
#set (CMAKE_MODULE_PATH
    #$ENV{SEQAN_ROOT}/util/cmake/FindTBB
#    $ENV{SEQAN_ROOT}/util/cmake)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/lib/seqan/util/cmake")
message (STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

set (SEQAN_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/lib/seqan/include")
message (STATUS "SEQAN_INCLUDE_PATH: ${SEQAN_INCLUDE_PATH}")

set (CMAKE_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/lib")
message (STATUS "CMAKE_INCLUDE_PATH = ${CMAKE_INCLUDE_PATH}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ----------------------------------------------------------------------------
# Prepare global variables
# ----------------------------------------------------------------------------

option (USE_UME_SIMD "Use UME::Simd if available." OFF)

# ----------------------------------------------------------------------------
# Include Src Directory
# ----------------------------------------------------------------------------

message (STATUS "Configuring src")
add_subdirectory(src)

# ===========================================================================
# Register Libraries (for GUI project generators)
# ===========================================================================

#seqan_setup_library ()

# ===========================================================================
# Include CPack
# ===========================================================================

#include (package)

# ===========================================================================
# Check Optional Libraries, Print Usage
# ===========================================================================

#include (FeatureSummary)

# ===========================================================================
# Setup Common Tool Description for Generic Workflow Nodes
# ===========================================================================

#include (SeqAnCtdSetup)
