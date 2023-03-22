cmake_minimum_required (VERSION 3.14)

enable_testing ()

set (PAIRALIGN_BENCHMARK_REPETITIONS "1" CACHE STRING "Set --benchmark_repetitions= for each benchmark.")

macro (pairalign_benchmark benchmark_cpp link_libraries)
    file (RELATIVE_PATH benchmark "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${benchmark_cpp}")
    pairalign_test_component (target "${benchmark}" TARGET_NAME)
    pairalign_test_component (test_name "${benchmark}" TEST_NAME)

    add_executable (${target} ${benchmark_cpp})
    target_link_libraries (${target} pairalign::benchmark ${link_libraries})
    add_test (NAME "${test_name}" COMMAND ${target} "--benchmark_repetitions=${PAIRALIGN_BENCHMARK_REPETITIONS}"
                                                    "--benchmark_out=${target}_out.json")

    unset (benchmark)
    unset (target)
    unset (test_name)
endmacro ()

# Download and build Googlebenchmark module.
pairalign_require_ccache ()
pairalign_require_benchmark ()

add_subdirectories ()
