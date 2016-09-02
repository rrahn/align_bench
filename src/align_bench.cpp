//#define SEQAN_ALIGN_SIMD_PROFILE 1
//#define SEQAN_ENABLE_TESTING 0
//#define SEQAN_ENABLE_DEBUG 0

#include <cxxabi.h>
#include <array>
#include <future>

#include <iostream>
#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/stream.h>
//#include <seqan/align.h>
#include <seqan/seq_io.h>

#include "align_bench_options.hpp"
#include "sequence_generator.hpp"
#include "benchmark_executor.hpp"

//#include <tbb/task.h>
//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range.h>
//#include <tbb/task_scheduler_init.h>
//#include <tbb/tick_count.h>
//
//long SerialFib( long n ) {
//    if( n<2 )
//        return n;
//    else
//        return SerialFib(n-1)+SerialFib(n-2);
//}
//
//class FibTask: public tbb::task {
//public:
//
//    const long n;
//    long* const sum;
//
//    FibTask( long n_, long* sum_ ) :
//        n(n_), sum(sum_)
//    {}
//    tbb::task* execute()
//    {      // Overrides virtual function task::execute
//        if( n < 4 ) {
//            *sum = SerialFib(n);
//        } else {
//            long x, y;
//            FibTask& a = *new( allocate_child() ) FibTask(n-1,&x);
//            FibTask& b = *new( allocate_child() ) FibTask(n-2,&y);
//            // Set ref_count to 'two children plus one for the wait".
//            set_ref_count(3);
//            // Start b running.
//            spawn( b );
//            // Start a running and wait for all children (a and b).
//            spawn_and_wait_for_all(a);
//            // Do the sum
//            *sum = x+y;
//        }
//        return NULL;
//    }
//};
//
//long ParallelFib( long n ) {
//    long sum;
//    FibTask& a = *new(tbb::task::allocate_root()) FibTask(n,&sum);
//    tbb::task::spawn_root_and_wait(a);
//    return sum;
//}

using namespace seqan;

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

    addOption(parser, seqan::ArgParseOption("ml", "min-length", "Minimal length of sequence.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "ml", "1000");
    setDefaultValue(parser, "ml", "1000");

    addOption(parser, seqan::ArgParseOption("xl", "max-lenght", "Maximal length of sequence.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "xl", "1000");
    setDefaultValue(parser, "xl", "1000");

    addOption(parser, seqan::ArgParseOption("r", "repetition", "Number of repeated runs.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "r", "1");
    setDefaultValue(parser, "r", "1");

    addOption(parser, seqan::ArgParseOption("p", "parallel-execution", "Parallelization strategy. If not set, execution is forced to serial mode", seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "p", "native native_vec omp omp_vec tbb tbb_vec");

    addOption(parser, seqan::ArgParseOption("t", "threads", "Number of threads", seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption("sw", "score-width", "Width of integers in bits used for score", seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "sw", "8 16 32 64");
    setDefaultValue(parser, "sw", "32");

    addOption(parser, seqan::ArgParseOption("bs", "block-size", "Size of blocks", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "bs", "100");

    addOption(parser, seqan::ArgParseOption("sa", "score-alphabet", "Size of blocks", seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "sa", "dna aminoacid");
    setDefaultValue(parser, "sa", "dna");

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

//    getArgumentValue(options.inputFileOne, parser, 0);
//    getArgumentValue(options.inputFileTwo, parser, 1);

    getOptionValue(options.alignOut, parser, "o");
    getOptionValue(options.rep, parser, "r");
    getOptionValue(options.minSize, parser, "ml");
    getOptionValue(options.maxSize, parser, "xl");

    std::string tmp;
    if (getOptionValue(tmp, parser, "p"))
    {
        if (tmp == "native")
            options.parMode = ParallelMode::NATIVE;
        else if (tmp == "native_vec")
            options.parMode = ParallelMode::NATIVE_VEC;
        else if (tmp == "omp")
            options.parMode = ParallelMode::OMP;
        else if (tmp == "omp_vec")
            options.parMode = ParallelMode::OMP_VEC;
        else if (tmp == "tbb")
            options.parMode = ParallelMode::TBB;
        else if (tmp == "tbb_vec")
            options.parMode = ParallelMode::TBB_VEC;
        else
            options.parMode = ParallelMode::SERIAL;
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

    if (!getOptionValue(options.threadCount, parser, "t"))
        options.threadCount = std::thread::hardware_concurrency();

    // Read block size.
    getOptionValue(options.blockSize, parser, "bs");

    return ArgumentParser::PARSE_OK;
}

template <typename... TArgs>
inline void invoke(AlignBenchOptions & options,
                   TArgs &&... args)
{
//    std::cout << "Invoke Alignment...\t" << std::flush;
    BenchmarkExecutor device;
    device.runGlobalAlignment(options, std::forward<TArgs>(args)...);
//    std::cout << "\t\t\tdone." << std::endl;
//    device.printProfile(std::cout);
    options.stats.time = device.getTime();
}

template <typename... TArgs>
inline void configureExec(AlignBenchOptions & options,
                          TArgs &&... args)
{
    switch (options.parMode)
    {
        case ParallelMode::NATIVE:
        {
            options.stats.execPolicy = "par_native";
            invoke(options, std::forward<TArgs>(args)..., parNative);
            break;
        }
        case ParallelMode::NATIVE_VEC:
        {
            options.stats.execPolicy = "par_native_vec";
            options.stats.vectorLength = SEQAN_SIZEOF_MAX_VECTOR / static_cast<unsigned>(options.simdWidth);
            invoke(options, std::forward<TArgs>(args)..., parNativeVec);
            break;
        }
        case ParallelMode::OMP:
        {
            #if !defined(_OPENMP)
                SEQAN_ASSERT_FAIL("OpenMP not enabled.");
            #else
                options.stats.execPolicy = "par_omp";
                invoke(options, std::forward<TArgs>(args)..., parOmp);
            #endif
            break;
        }
        case ParallelMode::OMP_VEC:
        {
            #if !defined(_OPENMP)
                SEQAN_ASSERT_FAIL("OpenMP not enabled.");
            #else
                options.stats.execPolicy = "par_omp_vec";
                options.stats.vectorLength = SEQAN_SIZEOF_MAX_VECTOR / static_cast<unsigned>(options.simdWidth);
                invoke(options, std::forward<TArgs>(args)..., parOmpVec);
            #endif
            break;
        }
        case ParallelMode::TBB:
        {
            #if !defined(SEQAN_TBB)
                SEQAN_ASSERT_FAIL("TBB not enabled.");
            #else
                options.stats.execPolicy = "par_tbb";
                invoke(options, std::forward<TArgs>(args)..., parTbb);
            #endif
            break;
        }
        case ParallelMode::TBB_VEC:
        {
            #if !defined(SEQAN_TBB)
                SEQAN_ASSERT_FAIL("TBB not enabled.");
            #else
                options.stats.execPolicy = "par_tbb_vec";
                options.stats.vectorLength = SEQAN_SIZEOF_MAX_VECTOR / static_cast<unsigned>(options.simdWidth);
                invoke(options, std::forward<TArgs>(args)..., parTbbVec);
            #endif
            break;
        }
        default:
        {
            SEQAN_ASSERT(options.parMode == ParallelMode::SERIAL);
            options.stats.execPolicy = "serial";
            invoke(options, std::forward<TArgs>(args)..., ser);
        }
    }
}

template <typename TAlphabet, typename TScoreValue, typename... TArgs>
inline void
configureScore(AlignBenchOptions & options,
               TArgs &&... args)
{
    if (IsSameType<TAlphabet, Dna>::VALUE)
        configureExec(options, std::forward<TArgs>(args)..., Score<TScoreValue>(6, -4, -1 , -11));
    else
        configureExec(options, std::forward<TArgs>(args)..., Score<TScoreValue, ScoreMatrix<AminoAcid, ScoreSpecBlosum62> >(-1 , -11));
}

template <typename TAlphabet>
inline void
configureAlpha(AlignBenchOptions & options)
{
    SequenceGenerator<TAlphabet> gen;
    gen.setNumber(1);
    gen.setMinLength(options.minSize);
    gen.setMaxLength(options.maxSize);

    auto seqSet1 = gen.generate();
    auto seqSet2 = gen.generate();

    options.stats.seqHLength = length(seqSet1[0]);
    options.stats.seqVLength = length(seqSet2[0]);

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
            configureScore<TAlphabet, int64_t>(options, seqSet1, seqSet2);
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
        options.stats.scoreAlpha = "amino acid";
        configureAlpha<AminoAcid>(options);
    }
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
    auto handle = std::async(std::launch::async, configure, std::ref(options));

    auto res = handle.wait_for(std::chrono::minutes(40));
    if (res == std::future_status::ready)
        options.stats.state = "done";
    else
        options.stats.state = "timeout";
    options.stats.writeStats(std::cout);

    return EXIT_SUCCESS;
}