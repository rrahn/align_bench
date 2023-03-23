#include <benchmark/benchmark.h>

#include <ranges>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_local.hpp>
#include <pairwise_aligner/configuration/score_model_matrix_simd_saturated_NxN.hpp>

#include "pa_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_local_affine_matrix_simd_sat_ds_ho,
                     &AS500,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    constexpr size_t simd_lanes = simd_score<int8_t>::size_v;

    auto align_config =
        cfg::score_model_matrix_simd_saturated_NxN(
            cfg::method_local(gap_cost()),
            cost_matrix());

    run(state, cfg::configure_aligner(align_config), std::integral_constant<size_t, simd_lanes>{});
}

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_local_affine_unitary_simd_sat_ds_ht,
                     &ASUniProt,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    constexpr size_t simd_lanes = simd_score<int8_t>::size_v;

    auto align_config =
        cfg::score_model_matrix_simd_saturated_NxN(
            cfg::method_local(gap_cost()),
            cost_matrix());

    run(state, cfg::configure_aligner(align_config), std::integral_constant<size_t, simd_lanes>{});
}
BENCHMARK_MAIN();
