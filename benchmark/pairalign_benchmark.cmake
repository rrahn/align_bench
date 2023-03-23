cmake_minimum_required (VERSION 3.14)

include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)

# Initialise seqan2 library target
find_package (ZLIB)
find_package (BZip2)
find_package (OpenMP)
set (SEQAN_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/lib/seqan/include")
find_package (SeqAn REQUIRED HINTS "${CMAKE_SOURCE_DIR}/lib/seqan/util/cmake" NO_DEFAULT_PATH)

add_library (seqan2_seqan2 INTERFACE IMPORTED)
target_compile_definitions (seqan2_seqan2 INTERFACE ${SEQAN_DEFINITIONS})
target_compile_options (seqan2_seqan2 INTERFACE "${SEQAN_CXX_FLAGS_LIST}")
target_link_libraries (seqan2_seqan2 INTERFACE ${SEQAN_LIBRARIES})
target_include_directories (seqan2_seqan2 INTERFACE ${SEQAN_INCLUDE_DIRS})
add_library (seqan2::seqan2 ALIAS seqan2_seqan2)

# Initialise data root directory
if (NOT DATA_ROOT_DIR)
    set(DATA_ROOT_DIR "${CMAKE_BINARY_DIR}")
endif()

# Set directories for test output files, input data and binaries.
set (DATA_DIR "${DATA_ROOT_DIR}/data/")
message(STATUS "Using DATA_ROOT_DIR: ${DATA_ROOT_DIR}")

file (MAKE_DIRECTORY ${DATA_DIR})
file (MAKE_DIRECTORY ${DATA_ROOT_DIR}/output)

add_definitions (-DOUTPUTDIR=\"${DATA_ROOT_DIR}/output/\")
add_definitions (-DDATADIR=\"${DATA_DIR}\")
add_definitions (-DBINDIR=\"${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/\")

# Define cmake configuration flags to configure and build external projects with the same flags as specified for
# this project.
# set (APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "")
# list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
# list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
# list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
# list (APPEND APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# Set the seqan3 specific external project cmake args here as well, since we use the seqan3 external procject to
# fetch goolge test and others.
# set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "${APP_TEMPLATE_EXTERNAL_PROJECT_CMAKE_ARGS}")

# Add more cmake tooling from this project.
find_path (PAIRALIGN_CMAKE_MODULE_DIR NAMES pairalign_datasources.cmake HINTS "${CMAKE_SOURCE_DIR}/cmake/")
list(APPEND CMAKE_MODULE_PATH "${PAIRALIGN_CMAKE_MODULE_DIR}")

# find_path (SEQAN3_TEST_CMAKE_MODULE_DIR NAMES pairalign_test_component.cmake
#                                         HINTS "${CMAKE_SOURCE_DIR}/lib/seqan3/test/cmake/")
# list(APPEND CMAKE_MODULE_PATH "${SEQAN3_TEST_CMAKE_MODULE_DIR}")

# Build tests just before their execution, because they have not been built with "all" target.
# The trick is here to provide a cmake file as a directory property that executes the build command.
# file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake"
#             "execute_process(COMMAND ${CMAKE_COMMAND} --build . --target benchmark_test)")
# set_directory_properties (PROPERTIES TEST_INCLUDE_FILE "${CMAKE_CURRENT_BINARY_DIR}/build_test_targets.cmake")

# Define the test targets. All depending targets are built just before the test execution.
# add_custom_target (api_test)
# add_custom_target (cli_test)
# add_custom_target (benchmark_test)

# Test executables and libraries should not mix with the application files.
unset (CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
unset (CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset (CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# Define some helper interface libraries for the test
find_path (PAIRALIGN_BENCH_INCLUDE_DIR NAMES pairalign/benchmark/units.hpp
                                       HINTS "${CMAKE_SOURCE_DIR}/benchmark/include/")
message(STATUS "Found bench include dir: ${PAIRALIGN_BENCH_INCLUDE_DIR}")

#--------------------------------------------------------------------
# Cmake interface targets
#--------------------------------------------------------------------

add_library (pairalign_base INTERFACE)
target_include_directories (pairalign_base INTERFACE "${PAIRALIGN_BENCH_INCLUDE_DIR}")
target_compile_options (pairalign_base INTERFACE "-pedantic"  "-Wall" "-Wextra" "-Wno-error")
target_compile_features (pairalign_base INTERFACE cxx_std_20)
target_link_libraries (pairalign_base INTERFACE "seqan2::seqan2" "pthread")
add_library (pairalign::base ALIAS pairalign_base)

# add_library (pairalign_test_unit INTERFACE)
# target_include_directories (pairalign_test_unit INTERFACE "GTest::gtest" "GTest::gtest_main" "pairalign::test")
# target_link_libraries (pairalign_test_unit INTERFACE "GTest::gtest" "GTest::gtest_main" "pairalign::test")
# add_library (pairalign::test::unit ALIAS pairalign_test_unit)

add_library (pairalign_benchmark INTERFACE)
target_include_directories (pairalign_benchmark INTERFACE "benchmark_main" "benchmark" "pairalign::base")
target_link_libraries (pairalign_benchmark INTERFACE "benchmark_main" "benchmark" "pairalign::base")
add_library (pairalign::benchmark ALIAS pairalign_benchmark)

# add_library (pairalign_test_asan INTERFACE)
# target_compile_options(pairalign_test_asan INTERFACE "${JST_SANITIZER_FLAGS}")
# target_link_options(pairalign_test_asan INTERFACE "${JST_SANITIZER_FLAGS}")
# target_link_libraries (pairalign_test_asan INTERFACE "pairalign::test::unit")
# add_library (pairalign::test::asan ALIAS pairalign_test_asan)

# add_library (pairalign_test_tsan INTERFACE)
# target_compile_options(pairalign_test_tsan INTERFACE "${JST_SANITIZER_FLAGS}")
# target_link_options(pairalign_test_tsan INTERFACE "${JST_SANITIZER_FLAGS}")
# target_link_libraries (pairalign_test_tsan INTERFACE "pairalign::test::unit")
# add_library (pairalign::test::tsan ALIAS pairalign_test_tsan)

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules.
# ----------------------------------------------------------------------------

include (pairalign_datasources)
include (${CMAKE_SOURCE_DIR}/data/datasources.cmake)

include (pairalign_require_benchmark)
include (pairalign_require_ccache)
include (pairalign_add_subdirectories)
include (pairalign_test_component)
