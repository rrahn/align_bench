#include <benchmark/benchmark.h>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary.hpp>

#include "pa_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_semiglobal_affine_unitary_scalar,
                     seqan::Dna5,
                     int32_t)(benchmark::State& state) {
    using namespace seqan::pairwise_aligner;
    auto align_config =
        cfg::score_model_unitary(
            cfg::method_global(
                gap_cost(),
                cfg::leading_end_gap{.first_column = cfg::end_gap::free},
                cfg::trailing_end_gap{.last_column = cfg::end_gap::free}),
            match_cost(), mismatch_cost());

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

// Run the benchmark
BENCHMARK_MAIN();
