#pragma once

#include "pa_bench_fixture.hpp"

template <auto * data, typename alphabet_t, typename score_t = int32_t>
class pa_affine_single_bench_fixture : public pa_bench_fixture<data, alphabet_t, score_t> {
protected:
    using typename pa_bench_fixture<data, alphabet_t, score_t>::sequence_t;
public:
    using typename pa_bench_fixture<data, alphabet_t, score_t>::score_type;

    pa_affine_single_bench_fixture() = default;
    virtual ~pa_affine_single_bench_fixture() = default;

    template <typename aligner_t>
    void run(benchmark::State & state, aligner_t aligner) const {
        auto seq1 = this->sequence1();
        auto seq2 = this->sequence2();

        std::vector<score_type> results{};
        results.resize(seq1.size());
        for (auto _ : state) {
            for (int32_t i = 0; i < std::ranges::ssize(seq1); ++i) {
                results[i] = aligner.compute(seq1[i], seq2[i]).score();
            }
        }
        benchmark::DoNotOptimize(results);
    }
};
