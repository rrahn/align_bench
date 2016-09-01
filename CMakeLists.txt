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


message (STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

#set (SEQAN_INCLUDE_PATH $ENV{SEQAN_ROOT}/include)
message (STATUS "SEQAN_INCLUDE_PATH: ${SEQAN_INCLUDE_PATH}")

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
