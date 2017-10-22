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

#ifndef ALIGN_BENCH_PARSER_HPP_
#define ALIGN_BENCH_PARSER_HPP_

#include <seqan/arg_parse.h>
#include <seqan/stream.h>

#include "align_bench_options.hpp"

using namespace seqan;

auto toString = [](auto elem)
{
    std::stringstream tmp;
    tmp << elem;
    return tmp.str();
};

template <typename TParser>
void setup_parser(TParser & parser)
{
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "QUERY"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "DATABASE"));

    addOption(parser, seqan::ArgParseOption("o", "output", "Output file to write the alignment to.", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "align_bench_res.csv");

    addOption(parser, seqan::ArgParseOption("s", "simulate", "Number of sequences to be simulated.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "s", "-1");

    addOption(parser, seqan::ArgParseOption("m", "min-length", "Minimal simulated length.", seqan::ArgParseArgument::INTEGER, "INT"));
//    setMinValue(parser, "ml", "1000");
    setDefaultValue(parser, "m", "1000");

    addOption(parser, seqan::ArgParseOption("x", "max-length", "Maximal simulated length.", seqan::ArgParseArgument::INTEGER, "INT"));
//    setMinValue(parser, "xl", "1000");
    setDefaultValue(parser, "x", "1000");

    addOption(parser, seqan::ArgParseOption("", "pdf", "The probability distribution function to use.", seqan::ArgParseArgument::STRING, "PDF"));
    setValidValues(parser, "pdf", "uniform normal");
    setDefaultValue(parser, "pdf", "uniform");

    addOption(parser, seqan::ArgParseOption("r", "repetition", "Number of repeated runs.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "r", "1");
    setDefaultValue(parser, "r", "1");

    addOption(parser, seqan::ArgParseOption("i", "integer-width", "Width of integers in bits used for score", seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "i", "8 16 32 64");
    setDefaultValue(parser, "i", "32");

    addOption(parser, seqan::ArgParseOption("a", "alphabet", "Size of blocks", seqan::ArgParseArgument::STRING, "STRING"));
    setValidValues(parser, "a", "dna aa");
    setDefaultValue(parser, "a", "dna");

    addOption(parser, seqan::ArgParseOption("d", "dp-algorithm", "Alignment method", seqan::ArgParseArgument::STRING, "STRING"));
    setValidValues(parser, "d", "global semi overlap local");
    setDefaultValue(parser, "d", "global");

#if defined(ALIGN_BENCH_BANDED)
    addOption(parser, seqan::ArgParseOption("", "lower-diagonal", "Lower diagonal of band.", seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("", "upper-diagonal", "Upper diagonal of band.", seqan::ArgParseArgument::INTEGER, "INT"));
#endif // ALIGN_BENCH_BANDED

    addOption(parser, seqan::ArgParseOption("", "alignment-mode", "How the input sequences should be aligned", seqan::ArgParseArgument::STRING, "STRING"));
    setValidValues(parser, "alignment-mode", "pair search olc");
    setDefaultValue(parser, "alignment-mode", "pair");

    addOption(parser, seqan::ArgParseOption("", "sort-sequences", "Whether the input sequences should be sorted by their lengths"));

    addOption(parser, seqan::ArgParseOption("v", "vectorization", "If set, executes vectorized alignment code."));
}

template <typename TOptions, typename TParser>
inline void get_arguments(TOptions & options, TParser & parser)
{
    getArgumentValue(options.queryFile, parser, 0);
    getArgumentValue(options.databaseFile, parser, 1);

    getOptionValue(options.alignOut, parser, "o");
    getOptionValue(options.rep, parser, "r");
    getOptionValue(options.numSequences, parser, "s");
    getOptionValue(options.minSize, parser, "m");
    getOptionValue(options.maxSize, parser, "x");

    std::string tmp;

    clear(tmp);
    if (getOptionValue(tmp, parser, "pdf"))
    {
        if (tmp == "normal")
            options.distFunction = DistributionFunction::NORMAL_DISTRIBUTION;
        else if (tmp == "uniform")
            options.distFunction = DistributionFunction::UNIFORM_DISTRIBUTION;
    }

    std::string bitWidth;
    getOptionValue(bitWidth, parser, "i");
    if (bitWidth == "8")
        options.simdWidth = SimdIntegerWidth::BIT_8;
    else if (bitWidth == "16")
        options.simdWidth = SimdIntegerWidth::BIT_16;
    else if (bitWidth == "32")
        options.simdWidth = SimdIntegerWidth::BIT_32;
    else if (bitWidth == "64")
        options.simdWidth = SimdIntegerWidth::BIT_64;

    std::string alpha;
    getOptionValue(alpha, parser, "a");
    if (alpha == "dna")
        options.alpha = ScoreAlphabet::DNA;
    else
        options.alpha = ScoreAlphabet::AMINOACID;

    if (getOptionValue(options.stats.method, parser, "d"))
    {
        if (options.stats.method == "global")
            options.method = AlignMethod::GLOBAL;
        else if (options.stats.method == "semiglobal")
            options.method = AlignMethod::SEMIGLOBAL;
        else if (options.stats.method == "overlap")
            options.method = AlignMethod::OVERLAP;
        else if (options.stats.method == "local")
            options.method = AlignMethod::LOCAL;
        }

    clear(tmp);
    if (getOptionValue(tmp, parser, "alignment-mode"))
    {
        if (tmp == "pair")
            options.mode = AlignmentMode::PAIR;
        else if (tmp == "search")
            options.mode = AlignmentMode::SEARCH;
        else if (tmp == "olc")
            options.mode = AlignmentMode::OLC;
    }

    options.sortSequences = isSet(parser, "sort-sequences");
    options.simd = isSet(parser, "v");

#if defined(ALIGN_BENCH_BANDED)
    if (isSet(parser, "lower-diagonal") && isSet(parser, "upper-diagonal"))
    {
        options.isBanded = true;
        getOptionValue(options.lower, parser, "lower-diagonal");
        getOptionValue(options.upper, parser, "upper-diagonal");
    }
#endif // ALIGN_BENCH_BANDED
}

#endif // ALIGN_BENCH_PARSER_HPP_
