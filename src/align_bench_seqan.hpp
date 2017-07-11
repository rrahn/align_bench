// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef ALIGN_BENCH_SEQAN_HPP
#define ALIGN_BENCH_SEQAN_HPP

#include <seqan/basic.h>
#include <seqan/align_parallel_2.h>

#include "benchmark_executor.hpp"

using namespace seqan;


template <typename TAlignObject>
inline void writeAlignment(AlignBenchOptions const & options,
                           TAlignObject const & align)
{
    if (options.alignOut == "stdout")
    {
        std::cout << align << std::endl;
        return;
    }

    std::ofstream alignOut;
    alignOut.open(options.alignOut.c_str());
    if (!alignOut.good())
    {
        std::cerr << "Could not open file << " << options.alignOut.c_str() << ">>!" << std::endl;
        return;
    }
    alignOut << align;
    alignOut.close();
}

template <typename TStream, typename TResultVec>
inline void writeScores_(TStream & stream,
                        TResultVec const & vec)
{
    std::for_each(std::begin(vec), std::end(vec), [&](auto & sc) { stream << sc << ",\n"; });
    stream << "\n";
}

inline void writeScores(AlignBenchOptions const & options)
{
    if (options.alignOut == "stdout")
    {
        writeScores_(std::cout, options.stats.scores);
    } else
    {
        std::ofstream alignOut;
        alignOut.open(options.alignOut.c_str());
        if (!alignOut.good())
        {
            std::cerr << "Could not open file << " << options.alignOut.c_str() << ">>!" << std::endl;
            return;
        }
        writeScores_(alignOut, options.stats.scores);
    }
}

template <typename TExecPolicy,
          typename TSet1,
          typename TSet2,
          typename TScore>
inline void
BenchmarkExecutor::runAlignment(AlignBenchOptions & options,
                                TExecPolicy const & execPolicy,
                                TSet1 const & set1,
                                TSet2 const & set2,
                                TScore const & scoreMat)
{
    decltype(globalAlignmentScore(execPolicy, set1, set2, scoreMat)) res;
    switch (options.method)
    {
        case AlignMethod::GLOBAL:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
#if defined(ALIGN_BENCH_BANDED)
            options.stats.isBanded = (options.isBanded) ? "yes" : "no";
            if (options.isBanded)
                res = globalAlignmentScore(execPolicy, set1, set2, scoreMat, options.lower, options.upper);
            else
#endif // ALIGN_BENCH_BANDED
                res = globalAlignmentScore(execPolicy, set1, set2, scoreMat);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::LOCAL:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            res = localAlignmentScore(execPolicy, set1, set2, scoreMat);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::SEMIGLOBAL:
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            res = globalAlignmentScore(execPolicy, set1, set2, scoreMat, AlignConfig<true, false, false, true>{});
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
    }
    writeScores(options);
}

#endif // ALIGN_BENCH_SEQAN_HPP
