#include <benchmark/benchmark.h>

#include <algorithm>
#include <utility>

#include <seqan/sequence.h>

#include <pairalign/benchmark/units.hpp>

#include "fastq_reader.hpp"

template<typename sequence_t>
class bench_fixture : public benchmark::Fixture {

    using sequence_collection_type = seqan::StringSet<sequence_t>;

    sequence_collection_type _sequence_collection1;
    sequence_collection_type _sequence_collection2;

    public:
    void SetUp(const ::benchmark::State&) {

        pairalign::input::load_fastq(DATADIR"sim_reads_n1K_rl150.fq", _sequence_collection1);
        _sequence_collection2 = _sequence_collection1;
        std::ranges::rotate(_sequence_collection2, std::ranges::next(std::ranges::begin(_sequence_collection1)));
    }

    void TearDown(benchmark::State & state) {
        using alignment_instance_t = std::pair<sequence_t, sequence_t>;
        std::vector<alignment_instance_t> sequences{};
        sequences.resize(std::ranges::size(sequence1()));
        for (int32_t i = 0; i < std::ranges::ssize(sequences); ++i) {
            sequences[i] = std::pair{sequence1()[i], sequence2()[i]};
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
