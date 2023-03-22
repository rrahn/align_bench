#include <benchmark/benchmark.h>

#include <ranges>

#include <pairwise_aligner/configuration/configure_aligner.hpp>
#include <pairwise_aligner/configuration/method_global.hpp>
#include <pairwise_aligner/configuration/score_model_unitary_simd.hpp>

#include "pa_affine_bench_fixture.hpp"

BENCHMARK_TEMPLATE_F(pa_affine_bench_fixture,
                     pairwise_aligner_semiglobal_affine_unitary_simd_int32,
                     seqan::Dna5,
                     int32_t)(benchmark::State& state) {

    constexpr size_t simd_lanes = seqan::pairwise_aligner::simd_score<score_type>::size;
    auto align_config =
        seqan::pairwise_aligner::cfg::score_model_unitary_simd(
            seqan::pairwise_aligner::cfg::method_global(
                gap_cost(),
                seqan::pairwise_aligner::cfg::leading_end_gap{
                    .first_column = seqan::pairwise_aligner::cfg::end_gap::free
                },
                seqan::pairwise_aligner::cfg::trailing_end_gap{
                    .last_column = seqan::pairwise_aligner::cfg::end_gap::free
                }
            ),
            static_cast<score_type>(match_cost()), static_cast<score_type>(mismatch_cost())
        );

    auto aligner = seqan::pairwise_aligner::cfg::configure_aligner(align_config);

    auto seq1 = sequence1();
    auto seq2 = sequence2();

    std::vector<score_type> scores{};
    scores.resize(seq1.size());

    for (auto _ : state) {
        for (std::ptrdiff_t i = 0; i < std::ranges::ssize(scores); i += simd_lanes) {
            std::ranges::subrange chunk1{
                std::ranges::next(std::ranges::begin(seq1), i),
                std::ranges::next(std::ranges::begin(seq1), i + simd_lanes, std::ranges::end(seq1))
            };
            std::ranges::subrange chunk2{
                std::ranges::next(std::ranges::begin(seq2), i),
                std::ranges::next(std::ranges::begin(seq2), i + simd_lanes, std::ranges::end(seq2))
            };
            auto results = aligner.compute(chunk1, chunk2);
            for (int32_t j = 0; j < std::ranges::ssize(results); ++j) {
                scores[i] = results[j].score();
            }
        }
    }
    benchmark::DoNotOptimize(scores);
}

BENCHMARK_MAIN();
