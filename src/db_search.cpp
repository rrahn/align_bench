//#define SEQAN_ALIGN_SIMD_PROFILE 1
//#define SEQAN_ENABLE_TESTING 0
//#define SEQAN_ENABLE_DEBUG 0

#define DP_ALIGN_STATS

#include <cxxabi.h>
#include <array>
#include <future>

#ifdef DP_ALIGN_STATS
std::atomic<uint32_t> simdCounter;
std::atomic<uint32_t> serialCounter;
#endif

#include <iostream>


#include <seqan/basic.h>
#include <seqan/align.h>
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
    ArgumentParser parser("db_search");

    setShortDescription(parser, "Database search using naive search with dynamic programming.");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "DATABASE"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "QUERY"));

    addOption(parser, seqan::ArgParseOption("o", "output", "Output file to write the alignment to.", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "align_bench_res.csv");

    addOption(parser, seqan::ArgParseOption("e", "execution-mode", "Parallelization strategy. If not set, execution is forced to serial mode", seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "e", "seq par par_vec wave wave_vec");

    addOption(parser, seqan::ArgParseOption("t", "threads", "Number of threads", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "t", toString(std::thread::hardware_concurrency()));

    addOption(parser, seqan::ArgParseOption("p", "parallel-instances", "Number of parallel alignment instances", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "p", "1");
    setDefaultValue(parser, "p", toString(std::thread::hardware_concurrency() << 1));

    addOption(parser, seqan::ArgParseOption("s", "score-width", "Width of integers in bits used for score", seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "s", "16 32 64");
    setDefaultValue(parser, "s", "32");

    addOption(parser, seqan::ArgParseOption("b", "block-size", "Size of blocks", seqan::ArgParseArgument::INTEGER, "INTEGER"));
    setDefaultValue(parser, "b", "100");

    addOption(parser, seqan::ArgParseOption("a", "score-alphabet", "Size of blocks", seqan::ArgParseArgument::STRING, "STRING"));
    setValidValues(parser, "a", "dna aa");
    setDefaultValue(parser, "a", "dna");

    addOption(parser, seqan::ArgParseOption("m", "method", "Alignment method", seqan::ArgParseArgument::STRING, "STRING"));
    setValidValues(parser, "m", "global semi local");
    setDefaultValue(parser, "m", "global");

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

    getArgumentValue(options.databaseFile, parser, 0);
    getArgumentValue(options.queryFile, parser, 1);

    getOptionValue(options.alignOut, parser, "o");

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

    std::string bitWidth;
    getOptionValue(bitWidth, parser, "s");
    if (bitWidth == "16")
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
    getOptionValue(options.blockSize, parser, "b");

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
            invoke(options, std::forward<TArgs>(args)..., seq);
        }
    }
}

template <typename TScoreValue, typename... TArgs>
inline void
configureScore(Dna const & /*tag*/,
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
    configureExec(options, std::forward<TArgs>(args)..., Score<TScoreValue,
                  ScoreMatrix<AminoAcid, ScoreSpecBlosum62> >(-1 , -11));
}

template <typename TAlphabet>
inline void
configureAlpha(AlignBenchOptions & options)
{
    // If option for generating is specified.
    StringSet<String<TAlphabet>> seqSet1;
    StringSet<String<TAlphabet>> seqSet2;

    options.stats.totalCells = 0;

    StringSet<CharString> meta1;
    StringSet<CharString> meta2;

    // read query.
    try {
        std::cout << "Reading query";
        SeqFileIn queryFile{options.queryFile.c_str()};
        readRecords(meta1, seqSet2, queryFile);
    } catch(...)
    {
        std::cerr << "\nCould not open or read  query file" << std::endl;
        return;
    }
    std::cout << " ... done\n";
    // read database.
    try {
        std::cout << "Reading database";
        SeqFileIn dbFile{options.databaseFile.c_str()};
        readRecords(meta2, seqSet1, dbFile);
    } catch(ParseError & e)
    {
        std::cerr << e.what() << std::endl;
        return;
    } catch(...)
    {
        std::cerr << "\nCould not open or read database file" << std::endl;
        return;
    }
    std::cout << " ... done\n";
    // TODO(rmaerker): Add sort option, but do it with modified string.
    // We might want to have a sorted modifid string, with only views to the underlying container?
    // std::sort(begin(tmp1, Standard()), end(tmp1, Standard()), [](auto const & s1, auto const & s2){ return length(s1) < length(s2); });
    // std::sort(begin(tmp2, Standard()), end(tmp2, Standard()), [](auto const & s1, auto const & s2){ return length(s1) < length(s2); });
    //
    options.stats.seqMinLength = [](auto & set1, auto & set2)
    {
        decltype(length(set1)) m{0};
        for (auto && seq : set1)
        {
            m = std::min(length(seq), m);
        }

        for (auto && seq : set2)
        {
            m = std::min(length(seq), m);
        }
        return m;
    }(seqSet1, seqSet2);
    options.stats.seqMaxLength = [](auto & set1, auto & set2)
    {
        decltype(length(set1)) m{0};
        for (auto && seq : set1)
        {
            m = std::max(length(seq), m);
        }

        for (auto && seq : set2)
        {
            m = std::max(length(seq), m);
        }
        return m;
    }(seqSet1, seqSet2);

    options.stats.numAlignments = length(seqSet1) * length(seqSet2);

    StringSet<String<TAlphabet>, Dependent<Tight>> setH;
    StringSet<String<TAlphabet>, Dependent<Tight>> setV;

    std::cout << "Calculating total size";
    for (unsigned i = 0; i < length(seqSet1); ++i)
    {
        for (unsigned j = 0; j < length(seqSet2); ++j)
        {
            options.stats.totalCells += (1+length(seqSet1[i]))*(1+length(seqSet2[j]));
            appendValue(setH, seqSet1[i]);
            appendValue(setV, seqSet2[j]);
        }
    }

    std::cout << "\t ...done.\n";

    switch(options.simdWidth)
    {
        case SimdIntegerWidth::BIT_16:
            options.stats.scoreValue = "int16_t";
            configureScore<int16_t>(TAlphabet(), options, setH, setV);
            break;
        case SimdIntegerWidth::BIT_32:
            options.stats.scoreValue = "int32_t";
            configureScore<int32_t>(TAlphabet(), options, setH, setV);
            break;
        case SimdIntegerWidth::BIT_64:
            options.stats.scoreValue = "int64_t";
            configureScore<int64_t>(TAlphabet(), options, setH, setV);
            break;
    }
}

inline void
configure(AlignBenchOptions & options)
{
    if (options.alpha == ScoreAlphabet::DNA)
    {
        options.stats.scoreAlpha = "dna";
        configureAlpha<Dna>(options);
    }
    else
    {
        options.stats.scoreAlpha = "aa";
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
