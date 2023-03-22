
#include <filesystem>
#include <string>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

namespace pairalign::input {


inline seqan::StringSet<seqan::String<seqan::Dna5>> load_fastq(std::filesystem::path const & path) {
    seqan::StringSet<seqan::String<seqan::Dna5>> sequences{};

    std::string tmp_id{};
    std::string tmp_qual{};
    seqan::SeqFileIn seqFileIn(path.c_str());
    while (!seqan::atEnd(seqFileIn)) {
        seqan::resize(sequences, std::ranges::size(sequences) + 1);
        seqan::readRecord(tmp_id, seqan::back(sequences), tmp_qual, seqFileIn);
        seqan::clear(tmp_id);
        seqan::clear(tmp_qual);
    }
    return sequences;
}
} // namespace pairalign::seqan
