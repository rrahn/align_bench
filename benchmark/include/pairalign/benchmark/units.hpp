// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides predefined custom units for google benchmark.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <benchmark/benchmark.h>

#include <ranges>

namespace pairalign::units
{

/*!\brief This returns a counter which represents how many bytes were processed per second.
 *
 * \param  bytes The total number of bytes processed of a complete benchmark run.
 * \return       Returns a benchmark Counter which represents bytes/s.
 */
inline benchmark::Counter bytes_per_second(size_t bytes)
{
    return benchmark::Counter(bytes, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

//!\brief Calculates the number of cell updates for given sequences for a specific alignment config.
template <typename sequences_range_t>
inline size_t compute_total_cells(sequences_range_t const & sequences_range)
{
    auto count_cells = [&](auto && seq1, auto && seq2)
    {
        size_t const columns = std::ranges::size(seq1) + 1;
        size_t const rows = std::ranges::size(seq2) + 1;

        return columns * rows;
    };

    size_t matrix_cells = 0u;
    for (auto && [seq1, seq2] : sequences_range)
        matrix_cells += count_cells(seq1, seq2);

    return matrix_cells;
}

/*!\brief This returns a counter which represents how many cell updates were done for a matrix.
 *
 * \param  cells The total number of cells processed of a complete benchmark run.
 * \return       Returns a benchmark Counter which represents CUPS (cell updates per second).
 */
inline benchmark::Counter cups(size_t cells)
{
    return benchmark::Counter(cells, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1000);
}

template <typename ai_collection_t>
    requires std::ranges::range<ai_collection_t>
inline benchmark::Counter cups(ai_collection_t const & alignment_instances)
{
    return cups(compute_total_cells(alignment_instances));
}

} // namespace pairalign::units
