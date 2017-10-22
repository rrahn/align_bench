#include <iostream>
#include <type_traits>
#include <fstream>


#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/simd.h>

#include "align_bench_options.hpp"
#include "align_bench_parser.hpp"

inline ArgumentParser::ParseResult
parseCommandLine(AlignBenchOptions & options, int const argc, char* argv[])
{
    ArgumentParser parser("align_bench_par");

    setShortDescription(parser, "Alignment Benchmark Tool");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    setup_parser(parser);

    addOption(parser, seqan::ArgParseOption("t", "threads", "Number of threads", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "t", toString(std::thread::hardware_concurrency()));

    addOption(parser, seqan::ArgParseOption("", "lower-diagonal", "Lower diagonal of band.", seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("", "upper-diagonal", "Upper diagonal of band.", seqan::ArgParseArgument::INTEGER, "INT"));

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

    get_arguments(options, parser);

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

    getOptionValue(options.threadCount, parser, "t");
    if (isSet(parser, "lower-diagonal") && isSet(parser, "upper-diagonal"))
    {
        options.isBanded = true;
        getOptionValue(options.lower, parser, "lower-diagonal");
        getOptionValue(options.upper, parser, "upper-diagonal");
    }

    return ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function _createSimdRepImpl()
// ----------------------------------------------------------------------------

template <typename TBuffer, typename TSequences>
void _transform(TBuffer & buffer, TSequences const & sequences)
{
    using TSimdVector = std::decay_t<decltype(buffer[0][0])>;

    buffer.resize((length(sequences) + LENGTH<TSimdVector>::VALUE - 1)/LENGTH<TSimdVector>::VALUE);

    for (unsigned i = 0; i < length(buffer); ++i)
    {
        auto inf = infix(sequences, i*LENGTH<TSimdVector>::VALUE, (i+1)*LENGTH<TSimdVector>::VALUE);
        resize(buffer[i], length(inf[0]));
        _createSimdRepImpl(buffer[i], inf);
    }
}

int main(int argc, char* argv[])
{
    AlignBenchOptions options;

    if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK)
        return EXIT_FAILURE;

    using TSimdVector = typename SimdVector<int16_t>::Type;

    /*
     * Read the sequences.
     */

    std::cout << "Reading sequences ..." << std::flush;
    StringSet<CharString> meta1;
    StringSet<Dna5String> set1;
    StringSet<CharString> meta2;
    StringSet<Dna5String> set2;
    try {
        SeqFileIn queryFile{options.queryFile.c_str()};
        readRecords(meta1, set1, queryFile);
    } catch(seqan::ParseError & e)
    {
        std::cerr << "Could not read query file" << std::endl;
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    try {
        SeqFileIn dbFile{options.databaseFile.c_str()};
        readRecords(meta2, set2, dbFile);
    } catch(ParseError & e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    } catch(...)
    {
        std::cerr << "Database: " << options.databaseFile << "\n";
        std::cerr << "Could not read database file" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "\t\tdone.\n";

    double start = sysTime();
    /* prepare simd vector representation of sequences. */
    std::vector<String<TSimdVector, Alloc<OverAligned>>> blockBufferH;
    std::vector<String<TSimdVector, Alloc<OverAligned>>> blockBufferV;
    String<TSimdVector, Alloc<OverAligned>> results;

    _transform(blockBufferH, set1);
    _transform(blockBufferV, set2);
    std::cout << "Prepare time: " << sysTime() - start << "s.\n";

    start = sysTime();
    /*
     * Define score values.
     */
    TSimdVector gapOpen = createVector<TSimdVector>(-11);
    TSimdVector gapExtend = createVector<TSimdVector>(-1);
    TSimdVector match = createVector<TSimdVector>(6);
    TSimdVector mismatch = createVector<TSimdVector>(-4);
    TSimdVector infinity = createVector<TSimdVector>(minValue<int16_t>() / 2);

    // TSimdVector tmpD;
    // TSimdVector tmp;
    // TSimdVector colVert;

    String<TSimdVector, Alloc<OverAligned>> colDiag;
    String<TSimdVector, Alloc<OverAligned>> colHori;

    resize(colDiag, length(blockBufferV[0]) + 1, Exact());
    resize(colHori, length(blockBufferV[0]) + 1, Exact());

    /* The algorithm */
    for (unsigned i = 0; i < length(blockBufferH); ++i)
    {
        /* initialize the algorithm */
        auto & bufferH = blockBufferH[i];
        auto & bufferV = blockBufferV[i];

        colDiag[0] = createVector<TSimdVector>(0);
        colHori[0] = gapOpen;
        TSimdVector colVert    = gapOpen;

        for (unsigned i = 1; i < length(colDiag); ++i)
        {
            colDiag[i] = colVert;
            colHori[i] = colDiag[i] + gapOpen;
            colVert    += gapExtend;
        }

        /* start the recursion */
        for (unsigned col = 1; col < length(bufferH) + 1; ++col)
        {
            TSimdVector valueH = bufferH[col - 1];

            TSimdVector H = colDiag[0];
            colDiag[0] = colHori[0];
            colHori[0] += gapExtend;
            colVert = colDiag[0] + gapOpen;

            for (unsigned row = 1; row < length(bufferV) + 1; ++row)
            {
                // H += blend(mismatch, match, cmpEq(valueH, bufferV[row - 1]));
                H = max(H + blend(mismatch, match, cmpEq(valueH, bufferV[row - 1])), max(colHori[row],colVert));
                TSimdVector N = H;
                // colVert -= gapExtend;
                // colHori[row] -= gapExtend;
                // H -= gapOpen;
                colVert = max(colVert - gapExtend, H - gapOpen);
                colHori[row] = max(colHori[row] - gapExtend, H - gapOpen);
                // cache prevDiag
                H = colDiag[row];
                colDiag[row] = N;
                // colHori[row]    = max(colHori[row] + gapExtend, colDiag[row] + gapOpen);
                // colVert         = max(colVert + gapExtend, H + gapOpen);
                // TSimdVector tmp = H + blend(mismatch, match, cmpEq(valueH, bufferV[row - 1]));
                // H            = colDiag[row];
                // colDiag[row]    = max(tmp, max(colHori[row], colVert));
            }
        }
        appendValue(results, colDiag[length(bufferH)]);
    }
    std::cout << "Runtime: " << sysTime() - start << "s.\n";

    std::ofstream alignOut;
    alignOut.open(options.alignOut.c_str());
    if (!alignOut.good())
    {
        std::cerr << "Could not open file << " << options.alignOut.c_str() << ">>!" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::for_each(begin(results), end(results), [&](auto & sc)
    {
        for (unsigned i = 0; i < LENGTH<TSimdVector>::VALUE; ++i)
            alignOut << sc[i] << ",\n";
    });
}
