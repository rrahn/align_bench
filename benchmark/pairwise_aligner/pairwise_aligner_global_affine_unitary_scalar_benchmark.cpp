#include <benchmark/benchmark.h>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>

#include "pa_affine_single_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_single_bench_fixture,
                     pairwise_aligner_global_affine_unitary_scalar_ds_ho,
                     &DS150,
                     seqan::Dna5,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    auto align_config =
        cfg::score_model_unitary(
            cfg::method_global(
                gap_cost(),
                cfg::leading_end_gap{},
                cfg::trailing_end_gap{}
            ),
            match_cost(), mismatch_cost());

    run(state, cfg::configure_aligner(align_config));
}

BENCHMARK_TEMPLATE_F(pa_affine_single_bench_fixture,
                     pairwise_aligner_global_affine_unitary_scalar_ds_ht,
                     &DS400_800,
                     seqan::Dna5,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    auto align_config =
        cfg::score_model_unitary(
            cfg::method_global(
                gap_cost(),
                cfg::leading_end_gap{},
                cfg::trailing_end_gap{}
            ),
            match_cost(), mismatch_cost());

    run(state, cfg::configure_aligner(align_config));
}

BENCHMARK_MAIN();
