#pragma once

#include "pa_bench_fixture.hpp"

template <auto * data, typename alphabet_t, typename score_t = int32_t>
class pa_affine_profile_bench_fixture : public pa_bench_fixture<data, alphabet_t, score_t> {
protected:
    using typename pa_bench_fixture<data, alphabet_t, score_t>::sequence_t;
public:
    using typename pa_bench_fixture<data, alphabet_t, score_t>::score_type;

    pa_affine_profile_bench_fixture() = default;
    virtual ~pa_affine_profile_bench_fixture() = default;

    template <typename aligner_t, size_t simd_lanes_v>
    void run(benchmark::State & state, aligner_t aligner, std::integral_constant<size_t, simd_lanes_v>) const {
        auto seq1 = this->sequence1()[0];
        auto seq2 = this->sequence2();

        std::vector<score_type> scores{};
        scores.resize(seq2.size());

        for (auto _ : state) {
            for (std::ptrdiff_t i = 0; i < std::ranges::ssize(scores); i += simd_lanes_v) {
                std::ranges::subrange chunk2{
                    std::ranges::next(std::ranges::begin(seq2), i),
                    std::ranges::next(std::ranges::begin(seq2), i + simd_lanes_v, std::ranges::end(seq2))
                };
                auto results = aligner.compute(seq1, chunk2);
                for (int32_t j = 0; j < std::ranges::ssize(results); ++j) {
                    scores[i] = results[j].score();
                }
            }
        }
        benchmark::DoNotOptimize(scores);
    }
protected:

    void tearDownImpl(benchmark::State & state) const override {
        using alignment_instance_t = std::pair<sequence_t, sequence_t>;
        std::vector<alignment_instance_t> sequences{};
        sequences.resize(std::ranges::size(this->sequence1()));
        for (int32_t i = 0; i < std::ranges::ssize(sequences); ++i) {
            sequences[i] = std::pair{this->sequence1()[0], this->sequence2()[i]};
        }

        state.counters["#Cells"] = pairalign::units::compute_total_cells(sequences);
        state.counters["CUPS"] = pairalign::units::cups(sequences);
    }

};
