#include <benchmark/benchmark.h>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_matrix.hpp>

#include "pa_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_overlap_affine_matrix_scalar_as_ho,
                     &AS500,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    auto align_config =
        cfg::score_model_matrix(
            cfg::method_global(gap_cost(),
                               cfg::leading_end_gap{.first_column = cfg::end_gap::free, .first_row = cfg::end_gap::free},
                               cfg::trailing_end_gap{.last_column = cfg::end_gap::free, .last_row = cfg::end_gap::free}),
            cost_matrix());

    auto aligner = cfg::configure_aligner(align_config);

    auto seq1 = sequence1();
    auto seq2 = sequence2();

    std::vector<score_type> results{};

    results.resize(seq1.size());
    for (auto _ : state) {
        for (int32_t i = 0; i < std::ranges::ssize(seq1); ++i) {
           results[i] = aligner.compute(seq1[i], seq2[i]).score();
        }
    }
    benchmark::DoNotOptimize(results);
}

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_overlap_affine_matrix_scalar_as_ht,
                     &ASUniProt,
                     seqan::AminoAcid,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    auto align_config =
        cfg::score_model_matrix(
            cfg::method_global(gap_cost(),
                               cfg::leading_end_gap{.first_column = cfg::end_gap::free, .first_row = cfg::end_gap::free},
                               cfg::trailing_end_gap{.last_column = cfg::end_gap::free, .last_row = cfg::end_gap::free}),
            cost_matrix());

    auto aligner = cfg::configure_aligner(align_config);

    auto seq1 = sequence1();
    auto seq2 = sequence2();

    std::vector<score_type> results{};

    results.resize(seq1.size());
    for (auto _ : state) {
        for (int32_t i = 0; i < std::ranges::ssize(seq1); ++i) {
           results[i] = aligner.compute(seq1[i], seq2[i]).score();
        }
    }
    benchmark::DoNotOptimize(results);
}

BENCHMARK_MAIN();
