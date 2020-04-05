#include <vector>
#include <string>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/view_to_simd.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/views/pairwise_combine.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>

int main(int const argc, char const ** argv)
{
    using score_t = seqan3::simd_type_t<int16_t>;

    seqan3::sequence_file_input input{argv[1]};

    std::vector<seqan3::dna5_vector> sequences{};
    for (auto && rec : input)
        sequences.push_back(std::move(seqan3::get<seqan3::field::seq>(rec)));

    using time_t = typename std::chrono::high_resolution_clock::time_point;
    using duration_type = typename std::chrono::high_resolution_clock::duration;
    // duration_t duration{};

    duration_type duration_prep{};
    duration_type duration_algo{};
    duration_type duration_res{};

    time_t start_total = std::chrono::high_resolution_clock::now();

    std::vector<int32_t> results{};
    score_t gap_extension{seqan3::simd::fill<score_t>(-1)};
    score_t gap_open{seqan3::simd::fill<score_t>(-11)};
    score_t match{seqan3::simd::fill<score_t>(6)};
    score_t mismatch{seqan3::simd::fill<score_t>(-4)};
    size_t sequence_index = 0;

    std::vector<score_t, seqan3::aligned_allocator<score_t>> optimal_column{};
    std::vector<score_t, seqan3::aligned_allocator<score_t>> horizontal_column{};

    auto convert_to_simd = [] (auto && sequences)
    {
         std::vector<score_t, seqan3::aligned_allocator<score_t, alignof(score_t)>> simd_sequence{};
         simd_sequence.reserve(std::ranges::size(sequences[0]));

        for (auto && simd_vector_chunk : sequences | seqan3::views::to_simd<score_t>(0x8000))
            for (auto && simd_vector : simd_vector_chunk)
                simd_sequence.push_back(std::move(simd_vector));

        return simd_sequence;
    };

    // Iterate over the sequences
    // TODO: make chunks
    auto all_pairs = seqan3::views::pairwise_combine(sequences);

    using seq1_reference_t = std::tuple_element_t<0, std::ranges::range_reference_t<decltype(all_pairs)>>;
    using seq2_reference_t = std::tuple_element_t<1, std::ranges::range_reference_t<decltype(all_pairs)>>;

    std::vector<decltype(seqan3::views::type_reduce(std::declval<seq1_reference_t>()))> batch1{};
    batch1.reserve(seqan3::simd_traits<score_t>::length);
    std::vector<decltype(seqan3::views::type_reduce(std::declval<seq2_reference_t>()))> batch2{};
    batch2.reserve(seqan3::simd_traits<score_t>::length);

    auto it = all_pairs.begin();
    size_t chunk = 0;
    while (it != all_pairs.end())
    {
        time_t start = std::chrono::high_resolution_clock::now();

        batch1.clear();
        batch2.clear();
        for (size_t i = 0; i < seqan3::simd_traits<score_t>::length && it != all_pairs.end(); ++i, ++it)
        {
            batch1.push_back(seqan3::views::type_reduce(std::get<0>(*it)));
            batch2.push_back(seqan3::views::type_reduce(std::get<1>(*it)));
        }

        auto seq1_simd = convert_to_simd(batch1);
        auto seq2_simd = convert_to_simd(batch2);

        duration_prep += std::chrono::high_resolution_clock::now() - start;
        start = std::chrono::high_resolution_clock::now();

        // Initialise matrix
        optimal_column.clear();
        horizontal_column.clear();
        optimal_column.resize(seq2_simd.size() + 1, seqan3::simd::fill<score_t>(0));
        horizontal_column.resize(seq2_simd.size() + 1, seqan3::simd::fill<score_t>(0));

        score_t diagonal{seqan3::simd::fill<score_t>(0)};
        score_t vertical{gap_open};
        horizontal_column[0] = gap_open;

        // Initialise the first column.
        for (auto && [opt, hor] : seqan3::views::zip(optimal_column, horizontal_column) | seqan3::views::drop(1))
        {
            opt = vertical;
            hor = opt + gap_open;
            vertical += gap_extension;
        }

        // Compute the matrix
        for (auto it_col = seq1_simd.begin(); it_col != seq1_simd.end(); ++it_col)
        {
            // Initialise first cell of optimal_column.
            auto opt_it = optimal_column.begin();
            auto hor_it = horizontal_column.begin();

            diagonal = *opt_it;  // cache the diagonal for next cell
            *opt_it = *hor_it; // initialise the horizontal score
            *hor_it += gap_extension; // initialise the horizontal score
            vertical = *opt_it + gap_open; // initialise the vertical value
            // std::cout << "vert: " << *opt_it << "\n";
            // Go to next cell.
            ++opt_it;
            ++hor_it;
            // std::cout << "diagonal: ";
            for (auto it_row = seq2_simd.begin(); it_row != seq2_simd.end(); ++it_row, ++opt_it, ++hor_it)
            {
                // Precompute the diagonal score.
                score_t horizontal = *hor_it;
                score_t tmp = diagonal + (((*it_col ^ *it_row) <= seqan3::simd::fill<score_t>(0)) ? match : mismatch);

                // std::cout << diagonal << ": " <<   tmp << " ";

                tmp = (tmp < vertical) ? vertical : tmp;
                tmp = (tmp < horizontal) ? horizontal : tmp;

                // Store the current max score.
                diagonal = *opt_it; // cache the next diagonal before writing it
                *opt_it = tmp; // store the temporary result

                tmp += gap_open;  // add gap open costs
                vertical += gap_extension;
                horizontal += gap_extension;

                // store the vertical and horizontal value in the next path
                vertical = (vertical < tmp) ? tmp : vertical;
                *hor_it = (horizontal < tmp) ? tmp : horizontal;
            }
        }

        duration_algo += std::chrono::high_resolution_clock::now() - start;
        start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < std::ranges::size(batch1); ++i)
            results.push_back(optimal_column.back()[i]);

        duration_res += std::chrono::high_resolution_clock::now() - start;
    }

    auto duration_total = std::chrono::high_resolution_clock::now() - start_total;
    // auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Preparation:    " << std::chrono::duration_cast<std::chrono::milliseconds>(duration_prep).count() << " ms\n";
    std::cout << "Algorithm:      " << std::chrono::duration_cast<std::chrono::milliseconds>(duration_algo).count() << " ms\n";
    std::cout << "Extract result: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration_res).count()  << " ms\n";
    std::cout << "Total:          " << std::chrono::duration_cast<std::chrono::milliseconds>(duration_total).count()  << " ms\n";

    for (auto const & res : results)
        std::cerr << res << ",\n";

    std::cout << "\nDone!\n";
}
