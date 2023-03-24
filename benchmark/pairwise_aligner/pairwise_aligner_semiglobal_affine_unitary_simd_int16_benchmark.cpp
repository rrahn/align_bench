#include <benchmark/benchmark.h>

#include <ranges>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>

#include "pa_affine_bulk_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_bulk_bench_fixture,
                     pairwise_aligner_semiglobal_affine_unitary_simd_int16_ds_ho,
                     &DS150,
                     seqan::Dna5,
                     int16_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    constexpr size_t simd_lanes = simd_score<score_type>::size_v;
    auto align_config =
        cfg::score_model_unitary_simd(
            cfg::method_global(
                gap_cost(),
                cfg::leading_end_gap{
                    .first_column = cfg::end_gap::free
                },
                cfg::trailing_end_gap{
                    .last_column = cfg::end_gap::free
                }
            ),
            static_cast<score_type>(match_cost()), static_cast<score_type>(mismatch_cost())
        );

    run(state, cfg::configure_aligner(align_config), std::integral_constant<size_t, simd_lanes>{});
}

BENCHMARK_TEMPLATE_F(pa_affine_bulk_bench_fixture,
                     pairwise_aligner_semiglobal_affine_unitary_simd_int16_ds_ht,
                     &DS400_800,
                     seqan::Dna5,
                     int16_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    constexpr size_t simd_lanes = simd_score<score_type>::size_v;
    auto align_config =
        cfg::score_model_unitary_simd(
            cfg::method_global(
                gap_cost(),
                cfg::leading_end_gap{
                    .first_column = cfg::end_gap::free
                },
                cfg::trailing_end_gap{
                    .last_column = cfg::end_gap::free
                }
            ),
            static_cast<score_type>(match_cost()), static_cast<score_type>(mismatch_cost())
        );

    run(state, cfg::configure_aligner(align_config), std::integral_constant<size_t, simd_lanes>{});
}

// Run the benchmark
BENCHMARK_MAIN();
