//#define SEQAN_ALIGN_SIMD_PROFILE 1
//#define SEQAN_ENABLE_TESTING 0
//#define SEQAN_ENABLE_DEBUG 0

//#define DP_PARALLEL_SHOW_PROGRESS
//#define DP_ALIGN_STATS

#include <cxxabi.h>
#include <array>
#include <future>

#ifdef DP_ALIGN_STATS
std::atomic<uint32_t> simdCounter;
std::atomic<uint32_t> serialCounter;
#endif

#include <iostream>

#include <seqan/basic.h>
#include <seqan/align_parallel_2.h>

#include <seqan/arg_parse.h>
#include <seqan/stream.h>

#include <seqan/seq_io.h>

#include "align_bench_options.hpp"
#include "sequence_generator.hpp"
#include "benchmark_executor.hpp"

using namespace seqan;

auto toString = [](auto elem)
{
    std::stringstream tmp;
    tmp << elem;
    return tmp.str();
};
/*
 * @fn parsCommandLine
 *
 * @brief Parses the command line arguments and options.
 *
 * @signature ParseResult parseCommandLine(options, argc, argv)
 * @param   options The @link AlignBenchOptions @endlink to be created.
 * @param   argc    The number of input arguments. Of type <tt>int</tt>
 * @param   argv    The argument values. Of type <tt>char **</tt>.
 *
 * @return ParseResult PARSE_OK on success, otherwise PARSE_ERROR.
 */
inline ArgumentParser::ParseResult
parseCommandLine(AlignBenchOptions & options, int const argc, char* argv[])
{
    ArgumentParser parser("align_bench");

    setShortDescription(parser, "Alignment Benchmark Tool");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "QUERY"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "DATABASE"));

    addOption(parser, seqan::ArgParseOption("o", "output", "Output file to write the alignment to.", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "align_bench_res.csv");

    addOption(parser, seqan::ArgParseOption("s", "simulate", "Number of sequences to be simulated.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "s", "-1");

    addOption(parser, seqan::ArgParseOption("ml", "min-length", "Minimal length of sequence.", seqan::ArgParseArgument::INTEGER, "INT"));
//    setMinValue(parser, "ml", "1000");
    setDefaultValue(parser, "ml", "1000");

    addOption(parser, seqan::ArgParseOption("xl", "max-lenght", "Maximal length of sequence.", seqan::ArgParseArgument::INTEGER, "INT"));
//    setMinValue(parser, "xl", "1000");
    setDefaultValue(parser, "xl", "1000");

    addOption(parser, seqan::ArgParseOption("d", "distribution-function", "The distribution functio to use.", seqan::ArgParseArgument::STRING, "PDF"));
    setValidValues(parser, "d", "uniform normal");
    setDefaultValue(parser, "d", "uniform");

    addOption(parser, seqan::ArgParseOption("r", "repetition", "Number of repeated runs.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "r", "1");
    setDefaultValue(parser, "r", "1");

    addOption(parser, seqan::ArgParseOption("e", "execution-mode", "Parallelization strategy. If not set, execution is forced to serial mode", seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "e", "seq par par_vec wave wave_vec");

    addOption(parser, seqan::ArgParseOption("t", "threads", "Number of threads", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "t", toString(std::thread::hardware_concurrency()));

    addOption(parser, seqan::ArgParseOption("p", "parallel-instances", "Number of parallel alignment instances", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "p", "1");
    setDefaultValue(parser, "p", toString(std::thread::hardware_concurrency() << 1));

    addOption(parser, seqan::ArgParseOption("sw", "score-width", "Width of integers in bits used for score", seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "sw", "8 16 32 64");
    setDefaultValue(parser, "sw", "32");

    addOption(parser, seqan::ArgParseOption("bs", "block-size", "Size of blocks", seqan::ArgParseArgument::INTEGER, "INTEGER"));
    setDefaultValue(parser, "bs", "100");

    addOption(parser, seqan::ArgParseOption("sa", "score-alphabet", "Size of blocks", seqan::ArgParseArgument::STRING, "STRING"));
    setValidValues(parser, "sa", "dna aminoacid");
    setDefaultValue(parser, "sa", "dna");

    addOption(parser, seqan::ArgParseOption("m", "method", "Alignment method", seqan::ArgParseArgument::STRING, "STRING"));
    setValidValues(parser, "m", "global semi local");
    setDefaultValue(parser, "m", "global");

    addOption(parser, seqan::ArgParseOption("", "alignment-mode", "How the input sequences should be aligned", seqan::ArgParseArgument::STRING, "STRING"));
    setValidValues(parser, "alignment-mode", "pair search olc");
    setDefaultValue(parser, "alignment-mode", "pair");

    addOption(parser, seqan::ArgParseOption("", "sort-sequences", "Whether the sequences should be sorted"));

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

    getArgumentValue(options.queryFile, parser, 0);
    getArgumentValue(options.databaseFile, parser, 1);

    getOptionValue(options.alignOut, parser, "o");
    getOptionValue(options.rep, parser, "r");
    getOptionValue(options.numSequences, parser, "s");
    getOptionValue(options.minSize, parser, "ml");
    getOptionValue(options.maxSize, parser, "xl");

    std::string tmp;
    if (getOptionValue(tmp, parser, "e"))
    {
        if (tmp == "seq")
            options.parMode = ParallelMode::SEQUENTIAL;
        else if (tmp == "par")
            options.parMode = ParallelMode::PARALLEL;
        else if (tmp == "par_vec")
            options.parMode = ParallelMode::PARALLEL_VEC;
        else if (tmp == "wave")
            options.parMode = ParallelMode::WAVEFRONT;
        else if (tmp == "wave_vec")
            options.parMode = ParallelMode::WAVEFRONT_VEC;
    }

    clear(tmp);
    if (getOptionValue(tmp, parser, "d"))
    {
        if (tmp == "normal")
            options.distFunction = DistributionFunction::NORMAL_DISTRIBUTION;
        else if (tmp == "uniform")
            options.distFunction = DistributionFunction::UNIFORM_DISTRIBUTION;
    }

    std::string bitWidth;
    getOptionValue(bitWidth, parser, "sw");
    if (bitWidth == "8")
        options.simdWidth = SimdIntegerWidth::BIT_8;
    else if (bitWidth == "16")
        options.simdWidth = SimdIntegerWidth::BIT_16;
    else if (bitWidth == "32")
        options.simdWidth = SimdIntegerWidth::BIT_32;
    else if (bitWidth == "64")
        options.simdWidth = SimdIntegerWidth::BIT_64;

    std::string alpha;
    getOptionValue(alpha, parser, "sa");
    if (alpha == "dna")
        options.alpha = ScoreAlphabet::DNA;
    else
        options.alpha = ScoreAlphabet::AMINOACID;

    if (getOptionValue(options.stats.method, parser, "m"))
    {
        if (options.stats.method == "global")
            options.method = AlignMethod::GLOBAL;
        else if (options.stats.method == "semiglobal")
            options.method = AlignMethod::SEMIGLOBAL;
        else if (options.stats.method == "local")
            options.method = AlignMethod::LOCAL;
        }

    getOptionValue(options.threadCount, parser, "t");
    getOptionValue(options.parallelInstances, parser, "p");

    // Read block size.
    getOptionValue(options.blockSize, parser, "bs");

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

    return ArgumentParser::PARSE_OK;
}

template <typename... TArgs>
inline void invoke(AlignBenchOptions & options,
                   TArgs &&... args)
{
    std::cout << "Invoke Alignment...\t" << std::flush;
    BenchmarkExecutor device;
    device.runGlobalAlignment(options, std::forward<TArgs>(args)...);
    std::cout << "\t\t\tdone." << std::endl;
    device.printProfile(std::cout);
    options.stats.time = device.getTime();
}

template <typename... TArgs>
inline void configureExec(AlignBenchOptions & options,
                          TArgs &&... args)
{
    switch (options.parMode)
    {
        case ParallelMode::PARALLEL:
        {
            options.stats.execPolicy = "parallel";
            options.stats.threads = options.threadCount;
            seqan::ExecutionPolicy<seqan::Parallel, seqan::Serial> exec;
            setNumThreads(exec, options.threadCount);
            invoke(options, std::forward<TArgs>(args)..., exec);
            break;
        }
        case ParallelMode::PARALLEL_VEC:
        {
            #ifdef SEQAN_SIMD_ENABLED
            options.stats.execPolicy = "parallel_vec";
            options.stats.vectorLength = SEQAN_SIZEOF_MAX_VECTOR / static_cast<unsigned>(options.simdWidth);
            options.stats.threads = options.threadCount;
            seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial> exec;
            setNumThreads(exec, options.threadCount);
            invoke(options, std::forward<TArgs>(args)..., exec);
            #endif
            break;
        }
        case ParallelMode::WAVEFRONT:
        {
                options.stats.execPolicy = "wavefront";
                seqan::ExecutionPolicy<seqan::WavefrontAlignment<>, seqan::Serial> exec;
                options.stats.threads = options.threadCount;
                options.stats.parallelInstances = options.parallelInstances;
                options.stats.blockSize = options.blockSize;
                setNumThreads(exec, options.threadCount);
                setParallelAlignments(exec, options.parallelInstances);
                setBlockSize(exec, options.blockSize);
                invoke(options, std::forward<TArgs>(args)..., exec);
            break;
        }
        case ParallelMode::WAVEFRONT_VEC:
        {
            #ifdef SEQAN_SIMD_ENABLED
            options.stats.execPolicy = "wavefront_vec";
            options.stats.vectorLength = SEQAN_SIZEOF_MAX_VECTOR / static_cast<unsigned>(options.simdWidth);
            options.stats.threads = options.threadCount;
            options.stats.parallelInstances = options.parallelInstances;
            options.stats.blockSize = options.blockSize;
            #if defined(SEQAN_BLOCK_OFFSET_OPTIMIZATION)
            seqan::ExecutionPolicy<seqan::WavefrontAlignment<BlockOffsetOptimization>, seqan::Vectorial> exec;
            #else
            seqan::ExecutionPolicy<seqan::WavefrontAlignment<>, seqan::Vectorial> exec;
            #endif // defined(SEQAN_BLOCK_OFFSET_OPTIMIZATION)

            setNumThreads(exec, options.threadCount);
            setParallelAlignments(exec, options.parallelInstances);
            setBlockSize(exec, options.blockSize);
            invoke(options, std::forward<TArgs>(args)..., exec);
            #endif
            break;
        }
        default:
        {
            SEQAN_ASSERT(options.parMode == ParallelMode::SEQUENTIAL);
            options.stats.execPolicy = "sequential";
            options.stats.threads = 1;
            seqan::ExecutionPolicy<seqan::Serial, seqan::Serial> exec;
            invoke(options, std::forward<TArgs>(args)..., exec);
        }
    }
}

template <typename TScoreValue, typename... TArgs>
inline void
configureScore(Dna5 const & /*tag*/,
               AlignBenchOptions & options,
               TArgs &&... args)
{
    configureExec(options, std::forward<TArgs>(args)..., Score<TScoreValue>(6, -4, -1 , -11));
}

template <typename TScoreValue, typename... TArgs>
inline void
configureScore(AminoAcid const & /*tag*/,
               AlignBenchOptions & options,
               TArgs &&... args)
{
    configureExec(options, std::forward<TArgs>(args)..., Score<TScoreValue, ScoreMatrix<AminoAcid, ScoreSpecBlosum62> >(-1 , -11));
}

template <typename TAlphabet>
inline void
configureAlpha(AlignBenchOptions & options)
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
        SequenceGenerator<TAlphabet> gen;
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
    } else
    {
        std::cout << "Reading sequences ...";
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
            std::cerr << "8 bit values not supported. Continue with 16 bit values.\n";
        case SimdIntegerWidth::BIT_16:
            options.stats.scoreValue = "int16_t";
            configureScore<int16_t>(TAlphabet(), options, seqSet1, seqSet2);
            break;
        case SimdIntegerWidth::BIT_32:
            options.stats.scoreValue = "int32_t";
            configureScore<int32_t>(TAlphabet(), options, seqSet1, seqSet2);
            break;
        case SimdIntegerWidth::BIT_64:
            options.stats.scoreValue = "int64_t";
            //configureScore<int64_t>(TAlphabet(), options, seqSet1, seqSet2);
            break;
    }
}

inline void
configure(AlignBenchOptions & options)
{
    if (options.alpha == ScoreAlphabet::DNA)
    {
        options.stats.scoreAlpha = "dna";
        configureAlpha<Dna5>(options);
    }
    else
    {
        options.stats.scoreAlpha = "amino acid";
        configureAlpha<AminoAcid>(options);
    }
}

int main(int argc, char* argv[])
{
    AlignBenchOptions options;

    if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK)
        return EXIT_FAILURE;

// TODO(rrahn): Make object configurable.
#ifdef DP_ALIGN_STATS
    simdCounter = 0;
    serialCounter = 0;
#endif

    configure(options);
    options.stats.state = "done";
    options.stats.writeHeader(std::cout);
    options.stats.writeStats(std::cout);

    return EXIT_SUCCESS;
}
