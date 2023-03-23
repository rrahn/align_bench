#include <benchmark/benchmark.h>

#include <seqan/align_parallel.h>

#include "seqan_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(seqan_affine_bench_fixture,
                     seqan_global_affine_matrix_simd_int16_as_ho,
                     &DS150,
                     seqan::Dna5,
                     int16_t)(benchmark::State& state) {

    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;
    using end_gaps_t = seqan::AlignConfig<false, false, false, false>;

    seqan::String<score_type> res{};
    for (auto _ : state) {
        res = seqan::globalAlignmentScore(exec_policy_t{}, sequence1(), sequence2(), matrix_scheme(), end_gaps_t{});
    }

    benchmark::DoNotOptimize(res);
}

BENCHMARK_TEMPLATE_F(seqan_affine_bench_fixture,
                     seqan_global_affine_matrix_simd_int16_as_ht,
                     &DS400_800,
                     seqan::Dna5,
                     int16_t)(benchmark::State& state) {

    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;
    using end_gaps_t = seqan::AlignConfig<false, false, false, false>;

    seqan::String<score_type> res{};
    for (auto _ : state) {
        res = seqan::globalAlignmentScore(exec_policy_t{}, sequence1(), sequence2(), matrix_scheme(), end_gaps_t{});
    }

    benchmark::DoNotOptimize(res);
}

BENCHMARK_MAIN();
