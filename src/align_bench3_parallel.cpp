#include <chrono>
#include <iostream>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/pairwise_combine.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/range/views/zip.hpp>

int main(int const argc, char const ** argv)
{
    seqan3::sequence_file_input input_fasta{argv[1]};

    std::vector<seqan3::dna5_vector> sequences{};

    for (auto const & record : input_fasta)
        sequences.push_back(std::move(seqan3::get<seqan3::field::seq>(record)));

    seqan3::nucleotide_scoring_scheme scoring_scheme{seqan3::match_score{6}, seqan3::mismatch_score{-4}};
    seqan3::gap_scheme gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}};

    auto configuration = seqan3::align_cfg::mode{seqan3::global_alignment} |
                         seqan3::align_cfg::scoring{scoring_scheme} |
                         seqan3::align_cfg::gap{gap_scheme} |
                         seqan3::align_cfg::result{seqan3::with_score, seqan3::using_score_type<int16_t>} |
                         seqan3::align_cfg::vectorise;
                        //   |
                        //  seqan3::align_cfg::parallel{4};

    // auto vectorised_configuration = configuration | seqan3::align_cfg::vectorise;

    std::vector<int32_t> scores{};

    auto start_time = std::chrono::high_resolution_clock::now();

    auto sequence_pairs = seqan3::views::pairwise_combine(sequences);
    // Configure the alignment algorithm.
    // auto && [algorithm, adapted_config] = seqan3::detail::alignment_configurator::configure<decltype(sequence_pairs)>(configuration);

    // size_t index{0};
    // for (auto && [seq1, seq2] : sequence_pairs)
    //     scores.push_back(algorithm(index, seq1, seq2).score());

    // for (auto && indexed_sequence_pair : seqan3::views::zip(sequence_pairs, std::views::iota(0)) | seqan3::views::chunk(1000))
    //     for (auto const & res : algorithm(indexed_sequence_pair))
    //         scores.push_back(res.score());

        // (algorithm(indexed_sequence_pair)[0].score());
    for (auto const & alignment_result : seqan3::align_pairwise(sequence_pairs, configuration))
        scores.push_back(alignment_result.score());

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count() << "s\n";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms\n";

    for (auto const & res : scores)
        std::cerr << res << ",\n";

    std::cout << "\nDone!\n";
}
