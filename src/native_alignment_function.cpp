#include <vector>
#include <string>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/pairwise_combine.hpp>
#include <seqan3/range/views/pairwise_combine.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>

int main(int const argc, char const ** argv)
{

    seqan3::sequence_file_input input{argv[1]};

    std::vector<seqan3::dna5_vector> sequences{};
    for (auto && rec : input)
        sequences.push_back(std::move(seqan3::get<seqan3::field::seq>(rec)));

    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<int32_t> results{};
    int32_t gap_extension{-1};
    int32_t gap_open{-11};
    int32_t match{6};
    int32_t mismatch{-4};
    size_t sequence_index = 0;

    std::vector<int32_t> optimal_column{};
    std::vector<int32_t> horizontal_column{};

    // Iterate over the sequences
    for (auto && [seq1, seq2] : seqan3::views::pairwise_combine(sequences | seqan3::views::take(1000)))
    {
        // Initialise matrix
        optimal_column.clear();
        horizontal_column.clear();
        optimal_column.resize(seq2.size() + 1, 0);
        horizontal_column.resize(seq2.size() + 1, 0);
        int32_t index = 0;
        int32_t diagonal{};
        int32_t vertical{gap_open};

        // Initialise the first column.
        for (auto && [opt, hor] : seqan3::views::zip(optimal_column, horizontal_column) | seqan3::views::drop(1))
        {
            opt = vertical;
            hor = opt + gap_open;
            vertical += gap_extension;
        }

        // Compute the matrix
        for (auto it_col = seq1.begin(); it_col != seq1.end(); ++it_col)
        {
            // Initialise first cell of optimal_column.
            auto opt_it = optimal_column.begin();
            auto hor_it = horizontal_column.begin();

            diagonal = *opt_it;  // cache the diagonal for next cell
            *opt_it =  gap_open + gap_extension * (it_col - seq1.begin()); // initialise the horizontal score
            *hor_it = *opt_it; // initialise the horizontal score
            vertical = *opt_it + gap_open; // initialise the vertical value
            // std::cout << "vert: " << *opt_it << "\n";
            // Go to next cell.
            ++opt_it;
            ++hor_it;
            // std::cout << "diagonal: ";
            for (auto it_row = seq2.begin(); it_row != seq2.end(); ++it_row, ++opt_it, ++hor_it)
            {
                // Precompute the diagonal score.
                int32_t tmp = diagonal + ((*it_col == *it_row) ? match : mismatch);

                // std::cout << diagonal << ": " <<   tmp << " ";

                tmp = (tmp < vertical) ? vertical : tmp;
                tmp = (tmp < *hor_it) ? *hor_it : tmp;

                // Store the current max score.
                diagonal = *opt_it; // cache the next diagonal before writing it
                *opt_it = tmp; // store the temporary result

                tmp += gap_open;  // add gap open costs
                vertical += gap_extension;
                *hor_it += gap_extension;

                // store the vertical and horizontal value in the next path
                vertical = (vertical < tmp) ? tmp : vertical;
                *hor_it = (*hor_it < tmp) ? tmp : *hor_it;
            }
        }
        results.push_back(optimal_column.back());
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count() << "s\n";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms\n";

    for (auto const & res : results)
        std::cerr << res << ",\n";

    std::cout << "\nDone!\n";
}
