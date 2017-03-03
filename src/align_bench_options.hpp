// ==========================================================================
//                              Align_Bench
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

#ifndef ALIGN_BENCH_OPTIONS_H_
#define ALIGN_BENCH_OPTIONS_H_

#include <string>

enum class DistributionFunction : uint8_t
{
    UNIFORM_DISTRIBUTION,
    NORMAL_DISTRIBUTION,
};

enum class AlignMethod : uint8_t
{
    GLOBAL,
    SEMIGLOBAL,
    LOCAL
};

enum class ParallelMode : uint8_t
{
    SEQUENTIAL,
    PARALLEL,
    PARALLEL_VEC,
    WAVEFRONT,
    WAVEFRONT_VEC
};

enum class SimdIntegerWidth : uint8_t
{
    BIT_8 = 1,
    BIT_16 = 2,
    BIT_32 = 4,
    BIT_64 = 8
};

enum class ScoreAlphabet : uint8_t
{
    DNA,
    AMINOACID
};

struct AlignBenchStats
{
    std::string             execPolicy;
    std::string             state;
    int                     numSequences;
    size_t                  seqMinLength;
    size_t                  seqMaxLength;
    double                  totalCells;
    std::string             dist;
    std::string             method;
    std::string             scoreValue;
    std::string             scoreAlpha;
    std::vector<size_t>     scores;
    double                  time;

    size_t                  blockSize{0};
    size_t                  threads{0};
    size_t                  parallelInstances{0};
    size_t                  vectorLength{0};

    template <typename TStream>
    void writeHeader(TStream & stream)
    {
        stream << "Policy,";
        stream << "State,";
        stream << "Simulated,";
        stream << "SeqMin,";
        stream << "SeqMax,";
        stream << "#GCells,";
        stream << "Dist,";
        stream << "BitsPerScore,";
        stream << "Alphabet,";
        stream << "Method,";
        stream << "Time,";
        stream << "GCUPS,";
        stream << "BlockSize,";
        stream << "#Threads,";
        stream << "#Instances,";
        stream << "VectorSize,";
        stream << "#SerialBlocks,";
        stream << "#SimdBlocks,";
        stream << "Results\n";
    }

    template <typename TStream>
    void writeStats(TStream & stream)
    {

        stream << execPolicy   << "," <<
                  state        << "," <<
                  numSequences << "," <<
                  seqMinLength << "," <<
                  seqMaxLength << "," <<
                  totalCells   << "," <<
                  dist         << "," <<
                  scoreValue   << "," <<
                  scoreAlpha   << "," <<
                  method       << ",";
        stream << time << ",";
        stream << (static_cast<double>(totalCells) / time) / 1000000000.0 << ",";
        if (blockSize != 0)
            stream << blockSize << ",";
        else
            stream << "n/a,";
        if (threads != 0)
            stream << threads << ",";
        else
            stream << "n/a,";
        if (parallelInstances != 0)
            stream << parallelInstances << ",";
        else
            stream << "n/a,";
        if (vectorLength != 0)
            stream << vectorLength << ",";
        else
            stream << "n/a,";
        #ifdef DP_ALIGN_STATS
            stream << serialCounter.load() << "," ;
            stream << simdCounter.load() << "," ;
        #else
            stream << "n/a,";
            stream << "n/a,";
        #endif

        for (auto val : scores)
            stream << val << " ";
        stream << '\n';
    }
};

struct AlignBenchOptions
{
    std::string queryFile;
    std::string databaseFile;
    std::string alignOut;
    unsigned rep;
    unsigned threadCount;
    unsigned parallelInstances;
    unsigned blockSize;
    int      numSequences;
    unsigned minSize;
    unsigned maxSize;
    AlignMethod  method;
    ParallelMode parMode = ParallelMode::SEQUENTIAL;

    SimdIntegerWidth simdWidth;
    ScoreAlphabet alpha;

    DistributionFunction distFunction;

    AlignBenchStats stats;
};

#endif  // #ifndef ALIGN_BENCH_OPTIONS_H_
