#include <string_view>
#include <iostream>

#include <seqan/align.h>

#include <pairalign/benchmark/data_sources.hpp>

#include "seqan_bench_fixture.hpp"

template <auto * data, typename alphabet_t, typename score_t = int32_t>
class seqan_affine_profile_bench_fixture : public seqan_bench_fixture<alphabet_t, data> {

    using base_t = seqan_bench_fixture<alphabet_t, data>;
    using typename base_t::sequence_collection_type;
    using matrix_scheme_t = seqan::Score<score_t, seqan::ScoreMatrix<seqan::AminoAcid, seqan::ScoreSpecBlosum62>>;

public:

    seqan_affine_profile_bench_fixture() = default;
    virtual ~seqan_affine_profile_bench_fixture() = default;

    using score_type = score_t;
    matrix_scheme_t matrix_scheme() const noexcept {
        return matrix_scheme_t{-1, -11};
    }

    sequence_collection_type sequence1() const noexcept {
        sequence_collection_type new_collection = base_t::sequence1();
        for (std::ptrdiff_t i = 1; i < static_cast<std::ptrdiff_t>(seqan::length(new_collection)); ++i) {
            new_collection[i] = new_collection[0];
        }
        return new_collection;
    }
protected:

    virtual void tearDownImpl(benchmark::State & state) const override {
        using std_sequence = std::vector<alphabet_t>;
        using alignment_instance_t = std::pair<std_sequence, std_sequence>;
        std::vector<alignment_instance_t> sequences{};
        sequences.resize(seqan::length(sequence1()));

        sequence_collection_type seq1 = sequence1();
        sequence_collection_type seq2 = this->sequence2();

        for (int32_t i = 0; i < std::ranges::ssize(sequences); ++i) {
            std_sequence tmp1{};
            std_sequence tmp2{};
            tmp1.resize(seqan::length(seq1[i]));
            tmp2.resize(seqan::length(seq2[i]));
            seqan::arrayCopy(seqan::begin(seq1[i]), seqan::end(seq1[i]), tmp1.data());
            seqan::arrayCopy(seqan::begin(seq2[i]), seqan::end(seq2[i]), tmp2.data());
            sequences[i] = std::pair{std::move(tmp1), std::move(tmp2)};
        }

        state.counters["CUPS"] = pairalign::units::cups(sequences);
    }
};
