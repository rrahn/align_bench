#include <benchmark/benchmark.h>

#include <seqan/align_parallel.h>

#include "seqan_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(seqan_affine_bench_fixture,
                     seqan_semiglobal_affine_unitary_simd_int32_ds_ho,
                     &DS150,
                     seqan::Dna5,
                     int32_t)(benchmark::State& state) {

    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;
    using end_gaps_t = seqan::AlignConfig<false, true, true, false>;

    seqan::String<score_type> res{};
    for (auto _ : state) {
        res = seqan::globalAlignmentScore(exec_policy_t{}, sequence1(), sequence2(), scoring_scheme(), end_gaps_t{});
    }

    benchmark::DoNotOptimize(res);
}

BENCHMARK_TEMPLATE_F(seqan_affine_bench_fixture,
                     seqan_semiglobal_affine_unitary_simd_int32_ds_ht,
                     &DS400_800,
                     seqan::Dna5,
                     int32_t)(benchmark::State& state) {

    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;
    using end_gaps_t = seqan::AlignConfig<false, true, true, false>;

    seqan::String<score_type> res{};
    for (auto _ : state) {
        res = seqan::globalAlignmentScore(exec_policy_t{}, sequence1(), sequence2(), scoring_scheme(), end_gaps_t{});
    }

    benchmark::DoNotOptimize(res);
}

// Run the benchmark
BENCHMARK_MAIN();
