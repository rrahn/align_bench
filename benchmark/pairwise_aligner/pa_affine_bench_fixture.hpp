#include "pa_bench_fixture.hpp"

#include <pairwise_aligner/configuration/gap_model_affine.hpp>

#include <pairalign/benchmark/data_sources.hpp>
template <auto * data, typename alphabet_t, typename score_t = int32_t>
class pa_affine_bench_fixture : public pa_bench_fixture<data, alphabet_t> {

public:
    using score_type = score_t;

    score_type match_cost() const noexcept {
        return 5;
    }

    score_type mismatch_cost() const noexcept {
        return -4;
    }

    constexpr auto gap_cost() const noexcept {
        return seqan::pairwise_aligner::cfg::gap_model_affine(-10, -1);
    }

    template <typename aligner_t, size_t simd_lanes_v>
    void run(benchmark::State & state, aligner_t aligner, std::integral_constant<size_t, simd_lanes_v>) const {
        auto seq1 = this->sequence1();
        auto seq2 = this->sequence2();

        std::vector<score_type> scores{};
        scores.resize(seq1.size());

        for (auto _ : state) {
            for (std::ptrdiff_t i = 0; i < std::ranges::ssize(scores); i += simd_lanes_v) {
                std::ranges::subrange chunk1{
                    std::ranges::next(std::ranges::begin(seq1), i),
                    std::ranges::next(std::ranges::begin(seq1), i + simd_lanes_v, std::ranges::end(seq1))
                };
                std::ranges::subrange chunk2{
                    std::ranges::next(std::ranges::begin(seq2), i),
                    std::ranges::next(std::ranges::begin(seq2), i + simd_lanes_v, std::ranges::end(seq2))
                };
                auto results = aligner.compute(chunk1, chunk2);
                for (int32_t j = 0; j < std::ranges::ssize(results); ++j) {
                    scores[i] = results[j].score();
                }
            }
        }
        benchmark::DoNotOptimize(scores);
    }

};
