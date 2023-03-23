#include <benchmark/benchmark.h>

#include <ranges>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_matrix_simd_1xN.hpp>

#include "pa_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_overlap_affine_profile_simd_int32_as_ho,
                     &AS500,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    constexpr size_t simd_lanes = simd_score<score_type>::size_v;

    auto align_config =
        cfg::score_model_matrix_simd_1xN(
                cfg::method_global(gap_cost(),
                                   cfg::leading_end_gap{.first_column = cfg::end_gap::free, .first_row = cfg::end_gap::free},
                                   cfg::trailing_end_gap{.last_column = cfg::end_gap::free, .last_row = cfg::end_gap::free}),
                cost_matrix());

    run_profile(state, cfg::configure_aligner(align_config), std::integral_constant<size_t, simd_lanes>{});
}

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_overlap_affine_profile_simd_int32_as_ht,
                     &ASUniProt,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    constexpr size_t simd_lanes = simd_score<score_type>::size_v;

    auto align_config =
        cfg::score_model_matrix_simd_1xN(
                cfg::method_global(gap_cost(),
                                   cfg::leading_end_gap{.first_column = cfg::end_gap::free, .first_row = cfg::end_gap::free},
                                   cfg::trailing_end_gap{.last_column = cfg::end_gap::free, .last_row = cfg::end_gap::free}),
                cost_matrix());

    run_profile(state, cfg::configure_aligner(align_config), std::integral_constant<size_t, simd_lanes>{});
}

BENCHMARK_MAIN();
