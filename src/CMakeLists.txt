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

# Update the cmake flags required for SeqAn.

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
# Search SeqAn and select dependencies.
find_package (ZLIB)
find_package (BZip2)

message (STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

#set (TBB_INCLUDE_DIR "/opt/local/include/tbb/")
#set (TBB_LIBRARY "/opt/local/lib/")

#find_package (TBB)
find_package (OpenMP)
find_package (Umesimd)
find_package (SeqAn REQUIRED)

# Warn if TBB wasn't find.
#if (NOT TBB_FOUND)
#    message (STATUS "   TBB not found")
#else(NOT TBB_FOUND)
#    message (STATUS "   Using TBB Version ${TBB_VERSION}")
#    list (APPEND SEQAN_DEFINITIONS "${TBB_DEFINITIONS};-DSEQAN_TBB")
#    list (APPEND SEQAN_LIBRARIES "${TBB_LIBRARIES}")
#    list (APPEND SEQAN_INCLUDE_DIRS "${TBB_INCLUDE_DIRS}")
#endif (NOT TBB_FOUND)

if (NOT UMESIMD_FOUND)
    message (STATUS "   UMESimd not found")
    set (USE_UME_SIMD OFF CACHE INTERNAL "Use ume simd if available." FORCE)
else (NOT UMESIMD_FOUND)
    message (STATUS "   Using UME Simd Version ${UMESIMD_VERSION_STRING}")
    list (APPEND SEQAN_INCLUDE_DIRS "${UMESIMD_INCLUDE_DIR}")
endif (NOT UMESIMD_FOUND)

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

# ----------------------------------------------------------------------------
# machine specific optimizations
# ----------------------------------------------------------------------------

# Disable machine specific optimizations on 32 bit targets and when using MSVC compiler.
if ((CMAKE_SIZEOF_VOID_P EQUAL 4) OR COMPILER_MSVC)
    set (SEQAN_ARCH_SSE4 FALSE)
    set (SEQAN_ARCH_AVX2 FALSE)
    set (SEQAN_ARCH_AVX512_KNL FALSE)
    set (SEQAN_ARCH_AVX512_SKX FALSE)
    set (SEQAN_ARCH_AVX512_CNL FALSE)
endif()

# Enable SSE4 if AVX[\d]+ is set. (Other parts in our build system expect it
# to be set and it is basically the synonym for 'SIMD is enabled')
if (SEQAN_ARCH_AVX2 OR SEQAN_ARCH_AVX512_KNL OR SEQAN_ARCH_AVX512_SKX OR SEQAN_ARCH_AVX512_CNL)
    set (SEQAN_ARCH_SSE4 TRUE)
endif ()

#if ((NOT USE_UME_SIMD) AND (SEQAN_ARCH_AVX512_KNL OR SEQAN_ARCH_AVX512_SKX OR SEQAN_ARCH_AVX512_CNL))
#    message (STATUS "Requested target platform requires UMESimd. Fall back to SSE4.")
#    set (SEQAN_ARCH_AVX512_KNL FALSE)
#    set (SEQAN_ARCH_AVX512_SKX FALSE)
#    set (SEQAN_ARCH_AVX512_CNL FALSE)
#endif()

if (COMPILER_MSVC)
    message (STATUS "Compiler is msvc")
    set(SEQAN_SIMD_SSE4_OPTIONS /arch:AVX)
    set(SEQAN_SIMD_AVX2_OPTIONS /arch:AVX2)
    set(SEQAN_SIMD_AVX512_KNL_OPTIONS "")
    set(SEQAN_SIMD_AVX512_SKX_OPTIONS "")
    set(SEQAN_SIMD_AVX512_CNL_OPTIONS "")
elseif (COMPILER_WINTEL)
    message (STATUS "Compiler is wintel")
    set(SEQAN_SIMD_SSE4_OPTIONS /QxSSE4.2)
    set(SEQAN_SIMD_AVX2_OPTIONS /QxCORE-AVX2)
    set(SEQAN_SIMD_AVX512_KNL_OPTIONS /QxMIC-AVX512)
    set(SEQAN_SIMD_AVX512_SKX_OPTIONS /QxCORE-AVX512)
    set(SEQAN_SIMD_AVX512_CNL_OPTIONS /QxCORE-AVX512)
elseif (COMPILER_LINTEL)
    message (STATUS "Compiler is lintel")
    set(SEQAN_SIMD_SSE4_OPTIONS -xSSE4.2)
    set(SEQAN_SIMD_AVX2_OPTIONS -xCORE-AVX2)
    set(SEQAN_SIMD_AVX512_KNL_OPTIONS -xMIC-AVX512)
    set(SEQAN_SIMD_AVX512_SKX_OPTIONS -xCORE-AVX512)
    set(SEQAN_SIMD_AVX512_CNL_OPTIONS -xCORE-AVX512)
    set(SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wcontext-limit=1000")
else()
    message (STATUS "Compiler is gnu/clang")
    set(SEQAN_SIMD_SSE4_OPTIONS -msse4)
    set(SEQAN_SIMD_AVX2_OPTIONS -mavx2)
    set(SEQAN_SIMD_AVX512_KNL_OPTIONS "-mavx512f -mavx512cd -mavx512er -mavx512pf")
    set(SEQAN_SIMD_AVX512_SKX_OPTIONS "-mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl")
    set(SEQAN_SIMD_AVX512_CNL_OPTIONS "${SEQAN_SIMD_AVX512_SKX_OPTIONS} -mavx512ifma -mavx512vbmi")
endif()

if (SEQAN_ARCH_NATIVE)
    message (STATUS "Building binaries optimized for this specific CPU. They might not work elsewhere.")
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=native")
    if (COMPILER_LINTEL)
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -ipo -no-prec-div -fp-model fast=2 -xHOST")
    endif ()
elseif (SEQAN_ARCH_SSE4)

    if (NOT ${CMAKE_SYSTEM_NAME} STREQUAL "OpenBSD" AND NOT COMPILER_LINTEL)
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -mpopcnt")
    endif ()
    if (COMPILER_LINTEL)
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -ipo -no-prec-div -fp-model fast=2")
    endif ()

    if (SEQAN_ARCH_AVX512_CNL)
        message (STATUS "Building optimized binaries up to AVX512 CNL and POPCNT.")
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} ${SEQAN_SIMD_AVX512_CNL_OPTIONS}")
        set (SEQAN_DEFINITIONS "${SEQAN_DEFINITIONS};-DALIGN_BENCH_SIMD_AVX512_CNL")
    elseif (SEQAN_ARCH_AVX512_SKX)
        message (STATUS "Building optimized binaries up to AVX512 SKX and POPCNT.")
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} ${SEQAN_SIMD_AVX512_SKX_OPTIONS}")
        set (SEQAN_DEFINITIONS "${SEQAN_DEFINITIONS};-DALIGN_BENCH_SIMD_AVX512_SKX")
    elseif (SEQAN_ARCH_AVX512_KNL)
        message (STATUS "Building optimized binaries up to AVX512 KNL and POPCNT.")
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} ${SEQAN_SIMD_AVX512_KNL_OPTIONS}")
        set (BLOCK_OFFSET_OPTIMIZATION OFF CACHE INTERNAL "Select block offset optimization." FORCE)
        set (SEQAN_DEFINITIONS "${SEQAN_DEFINITIONS};-DALIGN_BENCH_SIMD_AVX512_KNL")
    elseif (SEQAN_ARCH_AVX2)
        message (STATUS "Building optimized binaries up to AVX2 and POPCNT.")
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} ${SEQAN_SIMD_AVX2_OPTIONS}")
        set (SEQAN_DEFINITIONS "${SEQAN_DEFINITIONS};-DALIGN_BENCH_SIMD_AVX2")
    else ()
        message (STATUS "Building optimized binaries up to SSE4 and POPCNT.")
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} ${SEQAN_SIMD_SSE4_OPTIONS}")
        set (SEQAN_DEFINITIONS "${SEQAN_DEFINITIONS};-DALIGN_BENCH_SIMD_SSE41")
    endif ()

    # Enable correct simd library.
    if (USE_UME_SIMD)
        message (STATUS "UME::simd library is enabled.")
        set (SEQAN_DEFINITIONS "${SEQAN_DEFINITIONS};-DSEQAN_UMESIMD_ENABLED")
    else ()
        message (STATUS "SeqAn::simd library is enabled")
        set (SEQAN_DEFINITIONS "${SEQAN_DEFINITIONS};-DSEQAN_SEQANSIMD_ENABLED")
    endif ()
endif ()

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
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SEQAN_CXX_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable(align_bench_seq align_bench_seq.cpp
                           align_bench_options.hpp
                           align_bench_parser.hpp
                           align_bench_configure.hpp
                           timer.hpp
                           sequence_generator.hpp
                           benchmark_executor.hpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (align_bench_seq ${SEQAN_LIBRARIES})

add_executable(align_bench_par align_bench_par.cpp
                           align_bench_options.hpp
                           align_bench_parser.hpp
                           align_bench_configure.hpp
                           timer.hpp
                           sequence_generator.hpp
                           benchmark_executor.hpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (align_bench_par ${SEQAN_LIBRARIES})

add_executable(align_bench_wave align_bench_wave.cpp
                           align_bench_options.hpp
                           align_bench_parser.hpp
                           align_bench_configure.hpp
                           timer.hpp
                           sequence_generator.hpp
                           benchmark_executor.hpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (align_bench_wave ${SEQAN_LIBRARIES})

add_executable(bam2fa bam_fasta_converter.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (bam2fa ${SEQAN_LIBRARIES})

add_executable(simd_test simd_test.cpp
                         align_bench_options.hpp
                         align_bench_parser.hpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (simd_test ${SEQAN_LIBRARIES})
