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
#ifndef ALIGN_BENCH_CONFIGURE_HPP
#define ALIGN_BENCH_CONFIGURE_HPP

#include <cxxabi.h>
#include <array>
#include <future>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include "align_bench_seqan.hpp"
#include "sequence_generator.hpp"
#include "align_bench_options.hpp"

using namespace seqan;

template <typename... TArgs>
inline void invoke(AlignBenchOptions & options,
                   TArgs &&... args)
{
    std::cout << "Invoke Alignment...\t" << std::flush;
    BenchmarkExecutor device;
#if defined(ALIGN_BENCH_BANDED)
    if (options.isBanded)
        device.runAlignmentBanded(options, std::forward<TArgs>(args)...);
    else
#endif // ALIGN_BENCH_BANDED
        device.runAlignment(options, std::forward<TArgs>(args)...);
    std::cout << "\t\t\tdone." << std::endl;
    device.printProfile(std::cout);
    options.stats.time = device.getTime();
}

template <typename TScoreValue, typename... TArgs>
inline void
configureScore(Dna5 const & /*tag*/,
               AlignBenchOptions & options,
               TArgs &&... args)
{
    invoke(options, std::forward<TArgs>(args)..., Score<TScoreValue>(6, -4, -1 , -11));
}

template <typename TScoreValue, typename... TArgs>
inline void
configureScore(AminoAcid const & /*tag*/,
               AlignBenchOptions & options,
               TArgs &&... args)
{
    invoke(options, std::forward<TArgs>(args)..., Score<TScoreValue, ScoreMatrix<AminoAcid, ScoreSpecBlosum62> >(-1 , -11));
}

template <typename TAlphabet, typename ...TArgs>
inline void
configureSequences(AlignBenchOptions & options,
                   TArgs && ...args)
{
    // If option for generating is specified.
    using TView = typename Infix<String<TAlphabet> const>::Type;
    StringSet<TView> seqSet1;
    StringSet<TView> seqSet2;

    StringSet<String<TAlphabet>> tmp1;
    StringSet<String<TAlphabet>> tmp2;

    options.stats.totalCells = 0;
    if (options.numSequences != -1)
    {
        std::cout << "Generate sequences ...";
        SequenceGenerator<Dna> gen;
        gen.setNumber(options.numSequences);
        gen.setDistribution(options.distFunction);
        gen.setMinLength(options.minSize);
        gen.setMaxLength(options.maxSize);


        tmp1 = gen.generate();
        tmp2 = gen.generate();

        options.stats.numSequences = options.numSequences;
        options.stats.seqMinLength = options.minSize;
        options.stats.seqMaxLength = options.maxSize;
        options.stats.dist         = (options.distFunction == DistributionFunction::NORMAL_DISTRIBUTION) ? "normal" : "uniform";

        for (unsigned i = 0; i < length(tmp1); ++i)
        {
            appendValue(seqSet1, infix(tmp1[i], 0, length(tmp1[i])));
            appendValue(seqSet2, infix(tmp2[i], 0, length(tmp2[i])));
            options.stats.totalCells += (1+length(tmp1[i]))*(1+length(tmp2[i]));
        }

       SeqFileOut ref_out("ref_seq.fa");
       for (unsigned i = 0; i < length(tmp1); ++i)
           writeRecord(ref_out, "seq1", tmp1[i]);

           SeqFileOut qry_out("qry_seq.fa");
           for (unsigned i = 0; i < length(tmp2); ++i)
               writeRecord(qry_out, "seq2", tmp2[i]);
    } else
    {
        std::cout << "Reading sequences ..." << std::flush;
        StringSet<CharString> meta1;
        StringSet<CharString> meta2;

        try {
            SeqFileIn queryFile{options.queryFile.c_str()};
            readRecords(meta1, tmp1, queryFile);
        } catch(seqan::ParseError & e)
        {
            std::cerr << "Could not read query file" << std::endl;
            std::cerr << e.what() << std::endl;
            return;
        }

//        SeqFileOut fa_out("sd_0001.fa");
//        for (unsigned i = 0; i < 1000; ++i)
//            writeRecord(fa_out, meta1[i], tmp1[i]);

        try {
            SeqFileIn dbFile{options.databaseFile.c_str()};
            readRecords(meta2, tmp2, dbFile);
        } catch(ParseError & e)
        {
            std::cerr << e.what() << std::endl;
            return;
        } catch(...)
        {
            std::cerr << "Database: " << options.databaseFile << "\n";
            std::cerr << "Could not read database file" << std::endl;
            return;
        }

        options.stats.sortSequences = "no";
        if (options.sortSequences)
        {
            options.stats.sortSequences = "yes";
            std::cout << "\t done.\nSorting Sequences ...";

            std::sort(begin(tmp1, Standard()), end(tmp1, Standard()), [](auto const & s1, auto const & s2){ return length(s1) < length(s2); });
            std::sort(begin(tmp2, Standard()), end(tmp2, Standard()), [](auto const & s1, auto const & s2){ return length(s1) < length(s2); });
        }

        std::cout << "\t done.\nGenerating Sequences ..." << std::flush;
        switch (options.mode)
        {
            case AlignmentMode::PAIR:
            {
                options.stats.mode = "pair";
                SEQAN_ASSERT_EQ(length(tmp1), length(tmp2));

                for (unsigned i = 0; i < _min(length(tmp1), length(tmp2)); ++i)
                {
                    appendValue(seqSet1, infix(tmp1[i], 0, length(tmp1[i])), Generous());
                    appendValue(seqSet2, infix(tmp2[i], 0, length(tmp2[i])), Generous());
                    options.stats.totalCells += (1+length(tmp1[i]))*(1+length(tmp2[i]));
                }
                break;
            }
            case AlignmentMode::SEARCH:
            {
                options.stats.mode = "search";
                for (unsigned i = 0; i < length(tmp1); ++i)
                {
                    for (unsigned j = 0; j < length(tmp2); ++j)
                    {
                        appendValue(seqSet1, infix(tmp1[i], 0, length(tmp1[i])), Generous());
                        appendValue(seqSet2, infix(tmp2[j], 0, length(tmp2[j])), Generous());
                        options.stats.totalCells += (1+length(tmp1[i]))*(1+length(tmp2[j]));
                    }
                }
                break;
            }
            case AlignmentMode::OLC:
            {
                options.stats.mode = "olc";
                for (unsigned i = 0; i < length(tmp1); ++i)
                {
                    for (unsigned j = i+1; j < length(tmp1); ++j)
                    {
                        appendValue(seqSet1, infix(tmp1[i], 0, length(tmp1[i])), Generous());
                        appendValue(seqSet2, infix(tmp1[j], 0, length(tmp2[j])), Generous());
                        options.stats.totalCells += (1+length(tmp1[i]))*(1+length(tmp1[j]));
                    }
                }
                break;
            }
            default:
                return;
        }
        options.stats.seqMinLength = length(tmp1);
        options.stats.seqMaxLength = length(tmp2);
        options.stats.numSequences = length(seqSet1);
    }
    options.stats.numAlignments = length(seqSet1);

    std::cout << "\t done.\n";

    switch(options.simdWidth)
    {
        case SimdIntegerWidth::BIT_8:
            options.stats.scoreValue = "int8_t";
            // configureScore<int8_t>(TAlphabet(), options, std::forward<TArgs>(args)..., seqSet1, seqSet2);
        case SimdIntegerWidth::BIT_16:
            options.stats.scoreValue = "int16_t";
            configureScore<int16_t>(TAlphabet(), options, std::forward<TArgs>(args)..., seqSet1, seqSet2);
            break;
        case SimdIntegerWidth::BIT_32:
            options.stats.scoreValue = "int32_t";
            configureScore<int32_t>(TAlphabet(), options, std::forward<TArgs>(args)..., seqSet1, seqSet2);
            break;
        case SimdIntegerWidth::BIT_64:
            options.stats.scoreValue = "int64_t";
            //configureScore<int64_t>(TAlphabet(), options, std::forward<TArgs>(args)..., seqSet1, seqSet2);
            break;
    }
}

template <typename ...TArgs>
inline void
configureAlpha(AlignBenchOptions & options, TArgs && ...args)
{
    if (options.alpha == ScoreAlphabet::DNA)
    {
        options.stats.scoreAlpha = "dna";
        configureSequences<Dna5>(options, std::forward<TArgs>(args)...);
    }
    else
    {
        options.stats.scoreAlpha = "aa";
        configureSequences<AminoAcid>(options, std::forward<TArgs>(args)...);
    }
}

#endif  // ALIGN_BENCH_CONFIGURE_HPP
