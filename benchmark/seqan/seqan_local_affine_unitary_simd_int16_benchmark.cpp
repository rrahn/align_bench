#include <benchmark/benchmark.h>

#include <seqan/align_parallel.h>

#include "seqan_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(seqan_affine_bench_fixture,
                     seqan_local_affine_unitary_simd_int16_ds_ho,
                     &DS150,
                     seqan::Dna5,
                     int16_t)(benchmark::State& state) {

    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;

    seqan::String<score_type> res{};
    for (auto _ : state) {
        res = seqan::localAlignmentScore(exec_policy_t{}, sequence1(), sequence2(), scoring_scheme());
    }

    benchmark::DoNotOptimize(res);
}

BENCHMARK_TEMPLATE_F(seqan_affine_bench_fixture,
                     seqan_local_affine_unitary_simd_int16_ds_ht,
                     &DS400_800,
                     seqan::Dna5,
                     int16_t)(benchmark::State& state) {

    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>;

    seqan::String<score_type> res{};
    for (auto _ : state) {
        res = seqan::localAlignmentScore(exec_policy_t{}, sequence1(), sequence2(), scoring_scheme());
    }

    benchmark::DoNotOptimize(res);
}

// Run the benchmark
BENCHMARK_MAIN();
