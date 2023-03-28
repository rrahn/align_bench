#include <benchmark/benchmark.h>

#include <seqan/align_parallel.h>

#include "seqan_affine_profile_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(seqan_affine_profile_bench_fixture,
                     pairwise_aligner_local_affine_profile_simd_int32_as_ho,
                     &AS500,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {

    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;

    seqan::String<score_type> res{};
    auto seq1 = sequence1();
    auto seq2 = sequence2();
    for (auto _ : state) {
        res = seqan::localAlignmentScore(exec_policy_t{}, seq1, seq2, matrix_scheme());
    }
}

BENCHMARK_TEMPLATE_F(seqan_affine_profile_bench_fixture,
                     pairwise_aligner_local_affine_profile_simd_int32_as_ht,
                     &ASUniProt,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {
    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;

    seqan::String<score_type> res{};
    auto seq1 = sequence1();
    auto seq2 = sequence2();
    for (auto _ : state) {
        res = seqan::localAlignmentScore(exec_policy_t{}, seq1, seq2, matrix_scheme());
    }
}

BENCHMARK_MAIN();
