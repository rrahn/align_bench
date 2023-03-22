cmake_minimum_required (VERSION 3.14)

## Create seqan library target

set (SEQAN_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/lib/seqan/include")
message (STATUS "SEQAN_INCLUDE_PATH: ${SEQAN_INCLUDE_PATH}")
find_package (SeqAn REQUIRED HINTS "${CMAKE_SOURCE_DIR}/lib/seqan/util/cmake" NO_DEFAULT_PATH)

add_library (seqan2_seqan2 INTERFACE IMPORTED)
target_compile_definitions (seqan2_seqan2 INTERFACE ${SEQAN_DEFINITIONS})
target_compile_options (seqan2_seqan2 INTERFACE "${SEQAN_CXX_FLAGS_LIST}")
target_link_libraries (seqan2_seqan2 INTERFACE ${SEQAN_LIBRARIES})
target_include_directories (seqan2_seqan2 INTERFACE ${SEQAN_INCLUDE_DIRS})
add_library (seqan2::seqan2 ALIAS seqan2_seqan2)

set (PAIRALIGN_BENCHMARK_REPETITIONS 10)

## Register benchmark sources

pairalign_benchmark(seqan_global_affine_unitary_scalar_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_global_affine_unitary_scalar_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_semiglobal_affine_unitary_scalar_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_semiglobal_affine_unitary_scalar_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_overlap_affine_unitary_scalar_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_overlap_affine_unitary_scalar_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_local_affine_unitary_scalar_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_local_affine_unitary_scalar_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_global_affine_unitary_simd_int32_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_global_affine_unitary_simd_int32_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_semiglobal_affine_unitary_simd_int32_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_semiglobal_affine_unitary_simd_int32_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_overlap_affine_unitary_simd_int32_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_overlap_affine_unitary_simd_int32_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_local_affine_unitary_simd_int32_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_local_affine_unitary_simd_int32_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_global_affine_unitary_simd_int16_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_global_affine_unitary_simd_int16_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_semiglobal_affine_unitary_simd_int16_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_semiglobal_affine_unitary_simd_int16_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_overlap_affine_unitary_simd_int16_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_overlap_affine_unitary_simd_int16_benchmark FILES sim_reads_n1K_rl150.fq)

pairalign_benchmark(seqan_local_affine_unitary_simd_int16_benchmark.cpp "seqan2::seqan2")
target_use_datasources (seqan_local_affine_unitary_simd_int16_benchmark FILES sim_reads_n1K_rl150.fq)
