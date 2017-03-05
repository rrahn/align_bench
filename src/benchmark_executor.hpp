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

#ifndef BENCHMARK_EXECUTOR_HPP_
#define BENCHMARK_EXECUTOR_HPP_

#include <seqan/basic.h>
#include <seqan/align_parallel_2.h>

using namespace seqan;

// ----------------------------------------------------------------------------
// Class Timer
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec = void>
struct Timer
{
    TValue _begin;
    TValue _end;

    unsigned _rep;

    Timer() :
    _begin(0),
    _end(0)
    {};
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setRep()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void setRep(Timer<TValue, TSpec> & timer, unsigned const rep)
{
    timer._rep = rep;
}

// ----------------------------------------------------------------------------
// Function start()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void start(Timer<TValue, TSpec> & timer)
{
    timer._begin = sysTime();
}

// ----------------------------------------------------------------------------
// Function stop()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void stop(Timer<TValue, TSpec> & timer)
{
    timer._end = sysTime();
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline TValue getValue(Timer<TValue, TSpec> const & timer)
{
    return (timer._end - timer._begin) / timer._rep;
}

template <typename TValue, typename TSpec>
inline TValue getValue(Timer<TValue, TSpec> & timer)
{
    return getValue(static_cast<Timer<TValue, TSpec> const &>(timer));
}

// ----------------------------------------------------------------------------
// Function operator<<
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TSpec>
inline TStream & operator<<(TStream & os, Timer<TValue, TSpec> const & timer)
{
    os << getValue(timer) << " sec";
    return os;
}

// ----------------------------------------------------------------------------
// Function printRuler()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void printRuler(TStream & os)
{
    os << std::endl
    << "================================================================================"
    << std::endl << std::endl;
}

// ----------------------------------------------------------------------------
// Class BenchmarkExecutor
// ----------------------------------------------------------------------------

/*!
 * @class BenchmarkExecutor
 * @headerfile benchmark_executor.hpp
 * @brief Configures and executes alingment benchmark.
 *
 * @signature class BenchmarkExecutor;
 *
 * Runnable object, which profiles alignment algorithms.
 */
class BenchmarkExecutor
{
public:

    template <typename TSet1, typename TSet2>
    inline void runGlobalAlignment(TSet1 const &, TSet2 const &, unsigned const, seqan::Parallel const &);

    template <typename TSet1, typename TSet2>
    inline void runGlobalAlignment(TSet1 const &, TSet2 const &, unsigned const, seqan::Sequential const &);

    template <typename TSet1, typename TSet2, typename TScore,
              typename TMethod,
              typename TParSpec, typename TVecSpec>
    inline void runGlobalAlignment(AlignBenchOptions &,
                                   TSet1 const &,
                                   TSet2 const &,
                                   TScore const &,
                                   TMethod const &,
                                   seqan::ExecutionPolicy<TParSpec, TVecSpec> &);

    template <typename TSet1, typename TSet2, typename TScore,
              typename TMethod,
              typename TParSpec, typename TVecSpec>
    inline void runGlobalAlignment(AlignBenchOptions &,
                                   TSet1 const &,
                                   TSet2 const &,
                                   TScore const &,
                                   TMethod const &,
                                   seqan::ExecutionPolicy<seqan::WavefrontAlignment<TParSpec>, TVecSpec> &);

    template <typename TSet1, typename TSet2, typename TScore, typename TMethod>
    inline void runGlobalAlignment(AlignBenchOptions &,
                                   TSet1 const &,
                                   TSet2 const &,
                                   TScore const &,
                                   TMethod const &,
                                   seqan::Sequential const &);

    template <typename TStream>
    inline void
    printProfile(TStream && stream)
    {
        printRuler(stream);
        stream << mTimer << std::endl;
        printRuler(stream);
    }

    inline auto
    getTime()
    {
        return getValue(mTimer);
    }

private:

    Timer<double>   mTimer;
};


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
    std::for_each(std::begin(vec), std::end(vec), [&](auto & sc) { stream << sc << ","; });
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

template <typename TSet1, typename TSet2>
inline void
BenchmarkExecutor::runGlobalAlignment(TSet1 const & set1, TSet2 const & set2, unsigned const rep, seqan::Parallel const & /*tag*/)
{
    // Current old interface!
    SEQAN_ASSERT_EQ(length(set1), length(set2));

    using TAlign = Align<typename Value<TSet1>::Type, ArrayGaps>;
    StringSet<TAlign> alignSet;
    resize(alignSet, length(set1), Exact());

    auto zipView = makeZipView(alignSet, set1, set2);

    for (auto tuple : zipView)
    {
        resize(rows(std::get<0>(tuple)), 2);
        assignSource(row(std::get<0>(tuple), 0), std::get<1>(tuple));
        assignSource(row(std::get<0>(tuple), 1), std::get<2>(tuple));
    }


//    for (auto& align : alignSet)
//        std::cout << align << std::endl;
    setRep(mTimer, rep);
    start(mTimer);
    for (unsigned i = 0; i < rep; ++i)
        auto scores = globalAlignment(alignSet, Blosum62(-2, -8), Gotoh());
    stop(mTimer);

//    std::cout << "Align Profile: " << std::endl;
//    std::cout << "Preprocessing: " << profile.preprTimer << "s" << std::endl;
//    std::cout << "Alignment: " << profile.alignTimer << "s" << std::endl;
//    std::cout << "Traceback: " << profile.traceTimer << "s" << std::endl;

    // Proposed new interface!
//    auto alignObj = createAlignArray(str, str, Score<>, AlignmentTraits);    // pairwise alignment, one-vs-one
//    auto alignObj = createAlign(str, str, Score<>, ArrayGaps);
//    auto alignObj = createAlign(str, str, Score<>, Graph);
//    auto alignObj = createAlign(str, str, Score<>, Fragments);
//    auto alignObj = createAlign(str, str , Score<>, CigarString);
//
//    auto alignObj = createSplitAlign(...);  // get parameters necessary for split alignment.
//    auto alignObj = createExtendAlign(...);  // get parameters necessary for extend alignment.
//    auto alignObj = createBCAlign(...);       // get parameters necessary for banded chain alignment.
//
//    auto alignObj = createMsaAlign(set);     // multiple alignment

    // Global function:
//    align(alignObj, Config<TTraits>(), ExecPolicy());
}

template <typename TSet1, typename TSet2, typename TScore, typename TMethod>
inline void
BenchmarkExecutor::runGlobalAlignment(AlignBenchOptions & options,
                                      TSet1 const & set1,
                                      TSet2 const & set2,
                                      TScore const & scoreMat,
                                      TMethod const &,
                                      seqan::Sequential const & execPolicy)
{
    using TAlign = Align<typename Value<TSet1>::Type, ArrayGaps>;

    // Current old interface!
    SEQAN_ASSERT_EQ(length(set1), length(set2));

    StringSet<TAlign> alignSet;
    resize(alignSet, length(set1), Exact());

    auto zipView = makeZipView(alignSet, set1, set2);

    for (auto tuple : zipView)
    {
        resize(rows(std::get<0>(tuple)), 2);
        assignSource(row(std::get<0>(tuple), 0), std::get<1>(tuple));
        assignSource(row(std::get<0>(tuple), 1), std::get<2>(tuple));
    }

//    unsigned score = MinValue<unsigned>::VALUE;
    setRep(mTimer, options.rep);
    start(mTimer); 
    for (unsigned i = 0; i < options.rep; ++i)
    {
        for (unsigned i = 0; i < length(set1); ++i)
        {
            auto score = globalAlignmentScore(set1[i], set2[i], scoreMat);
            options.stats.scores.push_back(score);
        }
    }
    stop(mTimer);
    writeScores(options);
//    std::cout << "Score: " << score << '\n';
//    writeAlignment(options, alignSet[0]);
}

template <typename TSet1, typename TSet2, typename TScore,
          typename TMethod,
          typename TParSpec, typename TVecSpec>
inline void
BenchmarkExecutor::runGlobalAlignment(AlignBenchOptions & options,
                                      TSet1 const & set1,
                                      TSet2 const & set2,
                                      TScore const & scoreMat,
                                      TMethod const &,
                                      seqan::ExecutionPolicy<TParSpec, TVecSpec> & execPolicy)
{
//    SEQAN_ASSERT_NOT(std::is_same<TMethod, seqan::LocalAlignment_<>>::value);
    setNumThreads(execPolicy, options.threadCount);
    //    setParallelAlignments(execPolicy, options.parallelInstances);
    //    setBlockSize(execPolicy, options.blockSize);

    options.stats.threads = options.threadCount;
    //    options.stats.parallelInstances = options.parallelInstances;
    //    options.stats.blockSize = options.blockSize;

    std::vector<String<int>> alignScores(options.threadCount);

    auto chunkStringSet = [](auto const& infixSet)
    {
        using TSeq = typename std::decay<decltype(infixSet[0])>::type;
        using TSet = seqan::StringSet<TSeq, seqan::Dependent<seqan::Tight>>;

        TSet s;
        for (unsigned i = 0; i < length(infixSet); ++i)
        {
            appendValue(s, infixSet[i]);
        }
        return s;
    };

    setRep(mTimer, options.rep);
    start(mTimer);
    omp_set_num_threads(options.threadCount);
    SEQAN_OMP_PRAGMA(parallel for shared(alignScores))
    for (unsigned thread = 0; thread < options.threadCount; ++thread)
    {
        auto workSetH = infix(set1, thread * (length(set1)/options.threadCount),
                              std::min(length(set1), (thread + 1) * (length(set1)/options.threadCount)));
        auto workSetV = infix(set2, thread * (length(set2)/options.threadCount),
                              std::min(length(set2), (thread + 1) * (length(set2)/options.threadCount)));

        auto wH = chunkStringSet(workSetH);
        auto wV = chunkStringSet(workSetV);

        alignScores[thread] = globalAlignmentScore(wH, wV, scoreMat);
    }
    stop(mTimer);

    for (unsigned i = 0; i < options.threadCount; ++i)
    {
        for (unsigned j = 0; j < length(alignScores[i]); ++j)
            options.stats.scores.push_back(alignScores[i][j]);
    }
    writeScores(options);
}

template <typename TSet1, typename TSet2, typename TScore,
          typename TMethod,
          typename TParSpec, typename TVecSpec>
inline void
BenchmarkExecutor::runGlobalAlignment(AlignBenchOptions & options,
                                      TSet1 const & set1,
                                      TSet2 const & set2,
                                      TScore const & scoreMat,
                                      TMethod const &,
                                      seqan::ExecutionPolicy<seqan::WavefrontAlignment<TParSpec>, TVecSpec> & execPolicy)
{
    // Current old interface!
    SEQAN_ASSERT_EQ(length(set1), length(set2));

    setNumThreads(execPolicy, options.threadCount);
    setParallelAlignments(execPolicy, options.parallelInstances);
    setBlockSize(execPolicy, options.blockSize);

    options.stats.threads = options.threadCount;
    options.stats.parallelInstances = options.parallelInstances;
    options.stats.blockSize = options.blockSize;

    struct DPConfigTraits : seqan::DPTraits::GlobalAffine
    {
        using TAlgorithmType = TMethod;
        using TTracebackType = seqan::TracebackOff;
    };

    using TDPSettings = seqan::DPSettings<TScore, DPConfigTraits>;

    TDPSettings settings;
    settings.mScoringScheme = scoreMat;

    options.stats.scores.resize(length(set1), minValue<int>());

    setRep(mTimer, options.rep);
    start(mTimer);

    seqan::impl::alignExecBatch(execPolicy, set1, set2, settings,
    [&](auto const id, auto const score)
    {
        options.stats.scores[id] = score;
//        std::cout << "Finished alignment " << id << " with score " << score << '\n';
    });
    stop(mTimer);

    writeScores(options);


    //    StringSet<TAlign> alignSet;
    //    resize(alignSet, length(set1), Exact());
    //
    //    auto zipView = makeZipView(alignSet, set1, set2);
    //
    //    for (auto tuple : zipView)
    //    {
    //        resize(rows(std::get<0>(tuple)), 2);
    //        assignSource(row(std::get<0>(tuple), 0), std::get<1>(tuple));
    //        assignSource(row(std::get<0>(tuple), 1), std::get<2>(tuple));
    //    }
    //
    //    setNumThreads(execPolicy, options.threadCount);
    //    options.stats.threads = options.threadCount;
    //    options.stats.blockSize = options.blockSize;
    ////    unsigned score = MinValue<unsigned>::VALUE;
    //    setRep(mTimer, options.rep);
    //    start(mTimer);
    //    for (unsigned i = 0; i < options.rep; ++i)
    //    {
    //        std::cout << "i: " << i << '\n';
    //        auto score = parallelAlign(execPolicy, row(alignSet[0], 0), row(alignSet[0], 1), scoreMat, options.blockSize);
    //        options.stats.scores.push_back(score);
    //    }
    //    stop(mTimer);
    ////    std::cout << "Score: " << score << '\n';
    //    writeAlignment(options, alignSet[0]);
}

#endif  // #ifndef BENCHMARK_EXECUTOR_HPP_
