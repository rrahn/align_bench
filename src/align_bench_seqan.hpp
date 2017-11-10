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
#include <seqan/align_parallel.h>

#include "benchmark_executor.hpp"

using namespace seqan;


template <typename TGapsH, typename TGapsV>
inline void writeAlignment(AlignBenchOptions const & options,
                           TGapsH const & gapsH,
                           TGapsV const & gapsV)
{
    auto printGaps = [&](auto & stream)
    {
        for (unsigned i = 0; i < length(gapsH); ++i)
        {
            stream << "Alignment no. " << i << std::endl;
            stream << "Score: " << options.stats.scores[i] << std::endl;
            stream << gapsH[i] << std::endl;
            stream << gapsV[i] << std::endl;
        }
    };
    if (options.alignOut == "stdout")
    {
        printGaps(std::cout);
        return;
    }

    std::ofstream alignOut;
    alignOut.open(options.alignOut.c_str());
    if (!alignOut.good())
    {
        std::cerr << "Could not open file << " << options.alignOut.c_str() << ">>!" << std::endl;
        return;
    }
    printGaps(alignOut);
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

#if defined(ALIGN_BENCH_TRACE)
template <typename TExecPolicy,
          typename TSet1,
          typename TSet2,
          typename TScore>
inline void
BenchmarkExecutor::runAlignmentTrace(AlignBenchOptions & options,
                                     TExecPolicy const & execPolicy,
                                     TSet1 & set1,
                                     TSet2 & set2,
                                     TScore const & scoreMat)
{
    options.stats.isBanded = "no";
    using TSeqH = typename Value<TSet1>::Type;
    using TSeqV = typename Value<TSet2>::Type;
    StringSet<Gaps<TSeqH>> gapsSet1;
    StringSet<Gaps<TSeqV>> gapsSet2;

    auto fillGaps = [](auto & gaps, auto & sequences)
    {
        resize(gaps, length(sequences), Exact{});
        for (unsigned i = 0; i < length(sequences); ++i)
        {
            assignSource(gaps[i], sequences[i]);
        }
    };
    fillGaps(gapsSet1, set1);
    fillGaps(gapsSet2, set2);

    switch (options.method)
    {
        case AlignMethod::GLOBAL:
        {
            resize(options.stats.scores, length(gapsSet1), Exact());
            start(mTimer);
            auto res = globalAlignment(execPolicy, gapsSet1, gapsSet2, scoreMat);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::LOCAL:
        {
            resize(options.stats.scores, length(gapsSet1), Exact());
            start(mTimer);
            auto res = localAlignment(execPolicy, gapsSet1, gapsSet2, scoreMat);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::SEMIGLOBAL:
        {
            resize(options.stats.scores, length(gapsSet1), Exact());
            start(mTimer);
            auto res = globalAlignment(execPolicy, gapsSet1, gapsSet2, scoreMat, AlignConfig<true, false, false, true>{});
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::OVERLAP:
        {
            resize(options.stats.scores, length(gapsSet1), Exact());
            start(mTimer);
            auto res = globalAlignment(execPolicy, gapsSet1, gapsSet2, scoreMat, AlignConfig<true, true, true, true>{});
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
    }
    writeAlignment(options, gapsSet1, gapsSet2);
}

template <typename TExecPolicy,
          typename TSet1,
          typename TSet2,
          typename TScore>
inline void
BenchmarkExecutor::runAlignmentBandedTrace(AlignBenchOptions & options,
                                           TExecPolicy const & execPolicy,
                                           TSet1 & set1,
                                           TSet2 & set2,
                                           TScore const & scoreMat)
{
    options.stats.isBanded = "yes";

    using TSeqH = typename Value<TSet1>::Type;
    using TSeqV = typename Value<TSet2>::Type;
    StringSet<Gaps<TSeqH>> gapsSet1;
    StringSet<Gaps<TSeqV>> gapsSet2;

    auto fillGaps = [](auto & gaps, auto & sequences)
    {
        resize(gaps, length(sequences), Exact{});
        for (unsigned i = 0; i < length(sequences); ++i)
        {
            assignSource(gaps[i], sequences[i]);
        }
    };
    fillGaps(gapsSet1, set1);
    fillGaps(gapsSet2, set2);

    switch (options.method)
    {
        case AlignMethod::GLOBAL:
        {
            resize(options.stats.scores, length(gapsSet1), Exact());
            start(mTimer);
            auto res = globalAlignment(execPolicy, gapsSet1, gapsSet2, scoreMat, options.lower, options.upper);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::LOCAL:
        {
            resize(options.stats.scores, length(gapsSet1), Exact());
            start(mTimer);
            auto res = localAlignment(execPolicy, gapsSet1, gapsSet2, scoreMat, options.lower, options.upper);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::SEMIGLOBAL:
        {
            resize(options.stats.scores, length(gapsSet1), Exact());
            start(mTimer);
            auto res = globalAlignment(execPolicy, gapsSet1, gapsSet2, scoreMat, AlignConfig<true, false, false, true>{}, options.lower, options.upper);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::OVERLAP:
        {
            resize(options.stats.scores, length(gapsSet1), Exact());
            start(mTimer);
            auto res = globalAlignment(execPolicy, gapsSet1, gapsSet2, scoreMat, AlignConfig<true, true, true, true>{}, options.lower, options.upper);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
    }
    writeAlignment(options, gapsSet1, gapsSet2);
}

#else // ALIGN_BENCH_TRACE
template <typename TExecPolicy,
          typename TSet1,
          typename TSet2,
          typename TScore>
inline void
BenchmarkExecutor::runAlignment(AlignBenchOptions & options,
                                TExecPolicy const & execPolicy,
                                TSet1 & set1,
                                TSet2 & set2,
                                TScore const & scoreMat)
{
    options.stats.isBanded = "no";
    switch (options.method)
    {
        case AlignMethod::GLOBAL:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            auto res = globalAlignmentScore(execPolicy, set1, set2, scoreMat);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::LOCAL:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            auto res = localAlignmentScore(execPolicy, set1, set2, scoreMat);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::SEMIGLOBAL:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            auto res = globalAlignmentScore(execPolicy, set1, set2, scoreMat, AlignConfig<true, false, false, true>{});
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::OVERLAP:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            auto res = globalAlignmentScore(execPolicy, set1, set2, scoreMat, AlignConfig<true, true, true, true>{});
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
    }
    writeScores(options);
}

#if defined(ALIGN_BENCH_BANDED)
template <typename TExecPolicy,
          typename TSet1,
          typename TSet2,
          typename TScore>
inline void
BenchmarkExecutor::runAlignmentBanded(AlignBenchOptions & options,
                                      TExecPolicy const & execPolicy,
                                      TSet1 const & set1,
                                      TSet2 const & set2,
                                      TScore const & scoreMat)
{
    options.stats.isBanded = "yes";
    switch (options.method)
    {
        case AlignMethod::GLOBAL:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            auto res = globalAlignmentScore(execPolicy, set1, set2, scoreMat, options.lower, options.upper);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::LOCAL:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            auto res = localAlignmentScore(execPolicy, set1, set2, scoreMat, options.lower, options.upper);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::SEMIGLOBAL:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            auto res = globalAlignmentScore(execPolicy, set1, set2, scoreMat, AlignConfig<true, false, false, true>{}, options.lower, options.upper);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
        case AlignMethod::OVERLAP:
        {
            resize(options.stats.scores, length(set1), Exact());
            start(mTimer);
            auto res = globalAlignmentScore(execPolicy, set1, set2, scoreMat, AlignConfig<true, true, true, true>{}, options.lower, options.upper);
            stop(mTimer);
            seqan::arrayMoveForward(begin(res, Standard()), end(res, Standard()), begin(options.stats.scores, Standard()));
            break;
        }
    }
    writeScores(options);
}
#endif // ALIGN_BENCH_BANDED
#endif // ALIGN_BENCH_TRACE
#endif // ALIGN_BENCH_SEQAN_HPP
