#pragma once

#include <filesystem>
#include <iterator>
#include <string>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

namespace pairalign::input {

template <typename alphabet_t>
inline seqan::StringSet<seqan::String<alphabet_t>> load_fastq(std::filesystem::path const & path) {
    using sequence_t = seqan::String<alphabet_t>;
    using collection_t = seqan::StringSet<sequence_t>;

    std::string tmp_id{};
    std::string tmp_qual{};
    seqan::SeqFileIn seqFileIn(path.c_str());
    collection_t sequences{};

    size_t const max_count = 10000;
    size_t counter = 0;

    while (!seqan::atEnd(seqFileIn) && counter < max_count) {
        sequence_t seq{};
        seqan::readRecord(tmp_id, seq, tmp_qual, seqFileIn);
        seqan::clear(tmp_id);
        seqan::clear(tmp_qual);
        if (seqan::length(seq) < 50)
            continue;

        ++counter;
        seqan::appendValue(sequences, std::move(seq));
    }
    return sequences;
}
} // namespace pairalign::seqan
