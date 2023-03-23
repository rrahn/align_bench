#include <benchmark/benchmark.h>

#include <algorithm>
#include <utility>
#include <string>
#include <vector>

#include <pairalign/benchmark/fastq_reader.hpp>
#include <pairalign/benchmark/units.hpp>

template<auto * data, typename alphabet_t>
class pa_bench_fixture : public benchmark::Fixture {

    using sequence_t = std::string;
    using sequence_collection_type = std::vector<sequence_t>;

    sequence_collection_type _sequence_collection1;
    sequence_collection_type _sequence_collection2;

    public:
    void SetUp(const ::benchmark::State&) {
        _sequence_collection1.clear();
        _sequence_collection2.clear();
        auto tmpColl = pairalign::input::load_fastq<alphabet_t>(*data);
        _sequence_collection1.resize(seqan::length(tmpColl));
        for (int32_t i = 0; i < std::ranges::ssize(_sequence_collection1); ++i) {
            seqan::String<char> tmp{tmpColl[i]};
            _sequence_collection1[i] = std::string{seqan::toCString(tmp), seqan::length(tmp)};
        }

        _sequence_collection2 = _sequence_collection1;
        std::ranges::rotate(_sequence_collection2, std::ranges::next(std::ranges::begin(_sequence_collection2)));
    }

    void TearDown(benchmark::State & state) {
        tearDownImpl(state);
    }

    sequence_collection_type const & sequence1() const noexcept {
        return _sequence_collection1;
    }

    sequence_collection_type const & sequence2() const noexcept {
        return _sequence_collection2;
    }

    protected:

        virtual void tearDownImpl(benchmark::State & state) {
            using alignment_instance_t = std::pair<sequence_t, sequence_t>;
            std::vector<alignment_instance_t> sequences{};
            sequences.resize(std::ranges::size(sequence1()));
            for (int32_t i = 0; i < std::ranges::ssize(sequences); ++i) {
                sequences[i] = std::pair{sequence1()[i], sequence2()[i]};
            }

            state.counters["CUPS"] = pairalign::units::cups(sequences);
        }
};
