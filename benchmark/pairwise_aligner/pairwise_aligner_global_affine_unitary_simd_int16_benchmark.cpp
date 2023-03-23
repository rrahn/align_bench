#include <benchmark/benchmark.h>

#include <ranges>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>

#include "pa_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_global_affine_unitary_simd_int16,
                     seqan::Dna5,
                     int16_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    constexpr size_t simd_lanes = simd_score<score_type>::size;

    auto align_config =
        cfg::score_model_unitary_simd(cfg::method_global(gap_cost(), cfg::leading_end_gap{}, cfg::trailing_end_gap{}),
                                      static_cast<score_type>(match_cost()), static_cast<score_type>(mismatch_cost()));

    run(state, cfg::configure_aligner(align_config), std::integral_constant<size_t, simd_lanes>{});
}

BENCHMARK_MAIN();
