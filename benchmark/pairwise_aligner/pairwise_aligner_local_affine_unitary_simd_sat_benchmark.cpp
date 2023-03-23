#include <benchmark/benchmark.h>

#include <ranges>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_local.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd_saturated.hpp>

#include "pa_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_local_affine_unitary_simd_sat,
                     seqan::Dna5,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    constexpr size_t simd_lanes = simd_score<int8_t>::size;
    auto align_config =
        cfg::score_model_unitary_simd_saturated(
            cfg::method_local(gap_cost()),
            static_cast<score_type>(match_cost()), static_cast<score_type>(mismatch_cost())
        );

    run(state, cfg::configure_aligner(align_config), std::integral_constant<size_t, simd_lanes>{});
}

// Run the benchmark
BENCHMARK_MAIN();
