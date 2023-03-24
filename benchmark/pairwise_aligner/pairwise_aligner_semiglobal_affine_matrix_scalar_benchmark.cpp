#include <benchmark/benchmark.h>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_matrix.hpp>

#include "pa_affine_single_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_single_bench_fixture,
                     pairwise_aligner_semiglobal_affine_matrix_scalar_as_ho,
                     &AS500,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    auto align_config =
        cfg::score_model_matrix(
            cfg::method_global(gap_cost(),
                               cfg::leading_end_gap{.first_column = cfg::end_gap::free},
                               cfg::trailing_end_gap{.last_column = cfg::end_gap::free}),
            cost_matrix());

    run(state, cfg::configure_aligner(align_config));
}

BENCHMARK_TEMPLATE_F(pa_affine_single_bench_fixture,
                     pairwise_aligner_semiglobal_affine_matrix_scalar_as_ht,
                     &ASUniProt,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    auto align_config =
        cfg::score_model_matrix(
            cfg::method_global(gap_cost(),
                               cfg::leading_end_gap{.first_column = cfg::end_gap::free},
                               cfg::trailing_end_gap{.last_column = cfg::end_gap::free}),
            cost_matrix());

    run(state, cfg::configure_aligner(align_config));
}

BENCHMARK_MAIN();
