#include <benchmark/benchmark.h>

#include <ranges>

#include <seqan/align_parallel.h>

#include "seqan_affine_profile_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(seqan_affine_profile_bench_fixture,
                     seqan_semiglobal_affine_profile_simd_int16_as_ho,
                     &AS500,
                     seqan::AminoAcid,
                     int16_t)(benchmark::State& state) {
    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;
    using end_gaps_t = seqan::AlignConfig<false, true, true, false>;

    seqan::String<score_type> res{};
    auto seq1 = sequence1();
    auto seq2 = sequence2();
    for (auto _ : state) {
        res = seqan::globalAlignmentScore(exec_policy_t{}, seq1, seq2, matrix_scheme(), end_gaps_t{});
    }
}

BENCHMARK_TEMPLATE_F(seqan_affine_profile_bench_fixture,
                     seqan_semiglobal_affine_profile_simd_int16_as_ht,
                     &ASUniProt,
                     seqan::AminoAcid,
                     int16_t)(benchmark::State& state) {
    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;
    using end_gaps_t = seqan::AlignConfig<false, true, true, false>;

    seqan::String<score_type> res{};
    auto seq1 = sequence1();
    auto seq2 = sequence2();
    for (auto _ : state) {
        res = seqan::globalAlignmentScore(exec_policy_t{}, seq1, seq2, matrix_scheme(), end_gaps_t{});
    }
}

BENCHMARK_MAIN();
