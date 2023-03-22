#include <benchmark/benchmark.h>

#include <seqan/align_parallel.h>

#include "affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(affine_bench_fixture,
                     seqan_global_affine_unitary_scalar,
                     seqan::Dna5,
                     int32_t)(benchmark::State& state) {

    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Serial>;
    using end_gaps_t = seqan::AlignConfig<false, false, false, false>;

    seqan::String<score_type> res{};
    for (auto _ : state) {
        res = seqan::globalAlignmentScore(exec_policy_t{}, sequence1(), sequence2(), scoring_scheme(), end_gaps_t{});
    }

    benchmark::DoNotOptimize(res);
}

BENCHMARK_MAIN();
