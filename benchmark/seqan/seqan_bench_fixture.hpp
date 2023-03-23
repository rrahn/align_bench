#include <benchmark/benchmark.h>

#include <algorithm>
#include <utility>

#include <seqan/sequence.h>

#include <pairalign/benchmark/fastq_reader.hpp>
#include <pairalign/benchmark/units.hpp>

template<typename alphabet_t, auto * data>
class seqan_bench_fixture : public benchmark::Fixture {

    using sequence_t = seqan::String<alphabet_t>;
    using sequence_collection_type = seqan::StringSet<sequence_t>;

    sequence_collection_type _sequence_collection1;
    sequence_collection_type _sequence_collection2;

    public:
    void SetUp(const ::benchmark::State&) {
        seqan::clear(_sequence_collection1);
        seqan::clear(_sequence_collection2);

        _sequence_collection1 = pairalign::input::load_fastq<alphabet_t>(*data);

        std::vector<sequence_t> tmp_collection2{};
        tmp_collection2.resize(seqan::length(_sequence_collection1));
        for (int32_t i = 0; i < std::ranges::ssize(tmp_collection2); ++i) {
            tmp_collection2[i] = _sequence_collection1[i];
        }
        std::ranges::rotate(tmp_collection2, std::ranges::next(std::ranges::begin(tmp_collection2)));

        seqan::reserve(_sequence_collection2, tmp_collection2.size());
        std::ranges::for_each(tmp_collection2, [&] (auto seq) {
            seqan::appendValue(_sequence_collection2, std::move(seq));
        });
    }

    void TearDown(benchmark::State & state) {
        using std_sequence = std::vector<alphabet_t>;
        using alignment_instance_t = std::pair<std_sequence, std_sequence>;
        std::vector<alignment_instance_t> sequences{};
        sequences.resize(seqan::length(sequence1()));

        for (int32_t i = 0; i < std::ranges::ssize(sequences); ++i) {
            std_sequence tmp1{};
            std_sequence tmp2{};
            tmp1.resize(seqan::length(sequence1()[i]));
            tmp2.resize(seqan::length(sequence2()[i]));
            seqan::arrayCopy(seqan::begin(sequence1()[i]), seqan::end(sequence1()[i]), tmp1.data());
            seqan::arrayCopy(seqan::begin(sequence2()[i]), seqan::end(sequence2()[i]), tmp2.data());
            sequences[i] = std::pair{std::move(tmp1), std::move(tmp2)};
        }

        state.counters["CUPS"] = pairalign::units::cups(sequences);
    }

    sequence_collection_type const & sequence1() const noexcept {
        return _sequence_collection1;
    }

    sequence_collection_type const & sequence2() const noexcept {
        return _sequence_collection2;
    }
};
