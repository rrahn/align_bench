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
    ArgumentParser parser("align_bench");

    setShortDescription(parser, "Alignment Benchmark Tool");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

//    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "SEQUENCE_FILE"));
//    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "SEQUENCE_FILE"));
    addOption(parser, seqan::ArgParseOption("o", "output", "Output file to write the alignment to.", seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "stdout");

    addOption(parser, seqan::ArgParseOption("s", "sequences", "Number of sequences to be generated.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "s", "1");
    setDefaultValue(parser, "s", "1");

    addOption(parser, seqan::ArgParseOption("ml", "min-length", "Minimal length of sequence.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "ml", "1000");
    setDefaultValue(parser, "ml", "1000");

    addOption(parser, seqan::ArgParseOption("xl", "max-lenght", "Maximal length of sequence.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "xl", "1000");
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

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

//    getArgumentValue(options.inputFileOne, parser, 0);
//    getArgumentValue(options.inputFileTwo, parser, 1);

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
            invoke(options, std::forward<TArgs>(args)..., par);
            break;
        }
        case ParallelMode::PARALLEL_VEC:
        {
            options.stats.execPolicy = "parallel_vec";
            options.stats.vectorLength = SEQAN_SIZEOF_MAX_VECTOR / static_cast<unsigned>(options.simdWidth);
            invoke(options, std::forward<TArgs>(args)..., parVec);
            break;
        }
        case ParallelMode::WAVEFRONT:
        {
                options.stats.execPolicy = "wavefront";
                seqan::ExecutionPolicy<seqan::WavefrontAlignment<>, seqan::Serial> wavePolicy;
                invoke(options, std::forward<TArgs>(args)..., wavePolicy);
            break;
        }
        case ParallelMode::WAVEFRONT_VEC:
        {
            options.stats.execPolicy = "wavefront_vec";
            options.stats.vectorLength = SEQAN_SIZEOF_MAX_VECTOR / static_cast<unsigned>(options.simdWidth);
            SEQAN_ASSERT_FAIL("Not implemented yet");
            seqan::ExecutionPolicy<seqan::WavefrontAlignment<>, seqan::Vectorial> wavePolicy;
            invoke(options, std::forward<TArgs>(args)..., wavePolicy);
            break;
        }
        default:
        {
            SEQAN_ASSERT(options.parMode == ParallelMode::SEQUENTIAL);
            options.stats.execPolicy = "sequential";
            invoke(options, std::forward<TArgs>(args)..., seq);
        }
    }
}

template <typename... TArgs>
inline void
configureMethod(AlignBenchOptions & options,
                TArgs &&... args)
{
    switch (options.method)
    {
        case AlignMethod::GLOBAL:
            configureExec(options, std::forward<TArgs>(args)..., GlobalAlignment_<>());
            break;
        case AlignMethod::LOCAL:
            configureExec(options, std::forward<TArgs>(args)..., LocalAlignment_<>());
            break;
        case AlignMethod::SEMIGLOBAL:
            configureExec(options, std::forward<TArgs>(args)..., GlobalAlignment_<seqan::FreeEndGaps_<True, False, True, False>>());
            break;
    }
}

template <typename TAlphabet, typename TScoreValue, typename... TArgs>
inline void
configureScore(AlignBenchOptions & options,
               TArgs &&... args)
{
    if (IsSameType<TAlphabet, Dna>::VALUE)
        configureMethod(options, std::forward<TArgs>(args)..., Score<TScoreValue>(6, -4, -1 , -11));
    else
        configureMethod(options, std::forward<TArgs>(args)..., Score<TScoreValue, ScoreMatrix<AminoAcid, ScoreSpecBlosum62> >(-1 , -11));
}

template <typename TAlphabet>
inline void
configureAlpha(AlignBenchOptions & options)
{
    std::cout << "Generate sequences ...";
    SequenceGenerator<TAlphabet> gen;
    gen.setNumber(options.numSequences);
    gen.setDistribution(options.distFunction);
    gen.setMinLength(options.minSize);
    gen.setMaxLength(options.maxSize);

    auto seqSet1 = gen.generate();
    auto seqSet2 = gen.generate();

    std::cout << "\t done.\n";

    options.stats.numSequences = options.numSequences;
    options.stats.seqMinLength = options.minSize;
    options.stats.seqMaxLength = options.maxSize;
    options.stats.dist         = (options.distFunction == DistributionFunction::NORMAL_DISTRIBUTION) ? "normal" : "uniform";

    options.stats.totalCells = 0;
    for (unsigned i = 0; i < length(seqSet1); ++i)
    {
        options.stats.totalCells += (1+length(seqSet1[i]))*(1+length(seqSet2[i]));
    }

    switch(options.simdWidth)
    {
        case SimdIntegerWidth::BIT_8:
            std::cerr << "8 bit values not supported. Continue with 16 bit values.\n";
            break;
        case SimdIntegerWidth::BIT_16:
            options.stats.scoreValue = "int16_t";
            configureScore<TAlphabet, int16_t>(options, seqSet1, seqSet2);
            break;
        case SimdIntegerWidth::BIT_32:
            options.stats.scoreValue = "int32_t";
            configureScore<TAlphabet, int32_t>(options, seqSet1, seqSet2);
            break;
        case SimdIntegerWidth::BIT_64:
            options.stats.scoreValue = "int64_t";
//            configureScore<TAlphabet, int64_t>(options, seqSet1, seqSet2);
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
//    else
//    {
//        options.stats.scoreAlpha = "amino acid";
//        configureAlpha<AminoAcid>(options);
//    }
}

int main(int argc, char* argv[])
{
    AlignBenchOptions options;

    if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK) 
        return EXIT_FAILURE;

// NOTE(rrahn): Later we might read from input files.
//    SeqFileIn seqFileOne;
//    if (!open(seqFileOne, options.inputFileOne.c_str()))
//    {
//        std::cerr << "Could not open file: " <<
//    }
//    loadSequence(setLeft, options.inputFileOne.c_str());
//    loadSequence(setRight, options.inputFileTwo.c_str());

//    seqSet1[0] = "ATGTGTTACTGGGAGTAGGTTCTCCACTCTCTTCCAGTTAGGCTTGTAAGCGCTAATCGTTCTTGGGAAAGCCGGCTAAACCTTTGGACCAGCTGCAGCGGTATGATGTTTCTCAGAATCTATCGGGAAACAAGACACTCGGCATTTTATTGGTACGACCATTAAGGGTGCTTTGGATTTGGTCAGGACGTCAGTGTAACGACGGATCGCCCACAGCGATCTTGTATCTCGGGGACTCGAAACCCAGAACAATTCATCTATTACCGGCAAACGACGAGCGGGCCAAGCTGGCATTAGCCCGATAACAAAGACCTCGATC";
//    seqSet2[0] = "TAACAGGTATTGTAGCGTTTAAAGGCTCGAGTCACAGACATGATATAAGATTGGTCTGTAACTGGGCATTCTGAAAGCAAAACAATCCAGTATAATGGTGGTAGGGTGTTTGATTGAAGGTAGAGTAGTTACATCCGTCACGCACGCCTACCATTGGACAAAGGGCGCTCGTACCTTCTGTTTGTCGTTATCACGAGATGGGACTCGCAAATCAGACCTTCTCTTAGGGTCTTTACGTCTGTCAGAACAGATCTGTTCCACTTTTGGTCGCATTTGGCAAGCTGGAAGCGAATTGGGGACCATAATAAATGCGATTGGTG";
//    std::cout << length(seqSet1[0]) << " \t";
//    std::cout << length(seqSet2[0]);
//    std::cout << "\t\t\tdone\n" << std::flush;

//    for (auto& str : seqSet)
//        std::cout << str << "\n";
//    std::cout << std::flush;

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
