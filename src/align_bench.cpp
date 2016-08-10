//#define SEQAN_ALIGN_SIMD_PROFILE 1
//#define SEQAN_ENABLE_TESTING 0
//#define SEQAN_ENABLE_DEBUG 0

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

    addOption(parser, seqan::ArgParseOption("r", "repetition", "Number of repeated runs.", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "r", "1");
    setDefaultValue(parser, "r", "1");

    addOption(parser, seqan::ArgParseOption("p", "parallel", "Parallelization strategy. If not set, execution is forced to serial mode", seqan::ArgParseArgument::STRING, "STR"));
    setValidValues(parser, "p", "native omp tbb");

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

//    getArgumentValue(options.inputFileOne, parser, 0);
//    getArgumentValue(options.inputFileTwo, parser, 1);

    getOptionValue(options.alignOut, parser, "o");
    getOptionValue(options.rep, parser, "r");

    std::string tmp;
    if (getOptionValue(tmp, parser, "p"))
        options.parMode = ((tmp == "native") ? ParallelMode::NATIVE : ((tmp == "omp") ? ParallelMode::OMP : ParallelMode::TBB));

    return ArgumentParser::PARSE_OK;
}

template <typename... TArgs>
inline void invoke(AlignBenchOptions const & options,
                   TArgs &&... args)
{
    std::cout << "Invoke Alignment...\t" << std::flush;
    BenchmarkExecutor device;

    device.runGlobalAlignment(options, std::forward<TArgs>(args)...);
    std::cout << "\t\t\tdone." << std::endl;
    device.printProfile(std::cout);
}

template <typename... TArgs>
inline void configureExec(AlignBenchOptions const & options,
                          TArgs &&... args)
{
    switch (options.parMode)
    {
        case ParallelMode::NATIVE:
        {
            invoke(options, std::forward<TArgs>(args)..., parNative);
            break;
        }
        case ParallelMode::OMP:
        {
            #if !defined(_OPENMP)
                SEQAN_ASSERT_FAIL("OpenMP not enabled.");
            #else
                invoke(options, std::forward<TArgs>(args)..., parOmp);
            #endif
            break;
        }
        case ParallelMode::TBB:
        {
            #if !defined(SEQAN_TBB)
                SEQAN_ASSERT_FAIL("TBB not enabled.");
            #else
                invoke(options, std::forward<TArgs>(args)..., parTbb);
            #endif
            break;
        }
        default:
        {
            SEQAN_ASSERT_EQ(options.parMode, ParallelMode::SERIAL);
            invoke(options, std::forward<TArgs>(args)..., ser);
        }
    }
}

int main(int argc, char* argv[])
{
    AlignBenchOptions options;

    if (parseCommandLine(options, argc, argv) != ArgumentParser::PARSE_OK)
        return EXIT_FAILURE;


    // TODO(rrahn): We need to read the scoring model from command line.
    // TODO(rrahn): We need to read the execution model from command line.
    // TODO(rrahn): We need to read the Alignment Configuration from command line

// NOTE(rrahn): Later we might read from input files.
//    SeqFileIn seqFileOne;
//    if (!open(seqFileOne, options.inputFileOne.c_str()))
//    {
//        std::cerr << "Could not open file: " <<
//    }
//    loadSequence(setLeft, options.inputFileOne.c_str());
//    loadSequence(setRight, options.inputFileTwo.c_str());

    SequenceGenerator<Dna> gen;
    gen.setNumber(1);
    gen.setMinLength(32000);
    gen.setMaxLength(35000);

    std::cout << "Start generating sequences...\t" << std::flush;
    auto seqSet1 = gen.generate();
    auto seqSet2 = gen.generate();
//    seqSet1[0] = "ATGTGTTACTGGGAGTAGGTTCTCCACTCTCTTCCAGTTAGGCTTGTAAGCGCTAATCGTTCTTGGGAAAGCCGGCTAAACCTTTGGACCAGCTGCAGCGGTATGATGTTTCTCAGAATCTATCGGGAAACAAGACACTCGGCATTTTATTGGTACGACCATTAAGGGTGCTTTGGATTTGGTCAGGACGTCAGTGTAACGACGGATCGCCCACAGCGATCTTGTATCTCGGGGACTCGAAACCCAGAACAATTCATCTATTACCGGCAAACGACGAGCGGGCCAAGCTGGCATTAGCCCGATAACAAAGACCTCGATC";
//    seqSet2[0] = "TAACAGGTATTGTAGCGTTTAAAGGCTCGAGTCACAGACATGATATAAGATTGGTCTGTAACTGGGCATTCTGAAAGCAAAACAATCCAGTATAATGGTGGTAGGGTGTTTGATTGAAGGTAGAGTAGTTACATCCGTCACGCACGCCTACCATTGGACAAAGGGCGCTCGTACCTTCTGTTTGTCGTTATCACGAGATGGGACTCGCAAATCAGACCTTCTCTTAGGGTCTTTACGTCTGTCAGAACAGATCTGTTCCACTTTTGGTCGCATTTGGCAAGCTGGAAGCGAATTGGGGACCATAATAAATGCGATTGGTG";
//    std::cout << length(seqSet1[0]) << "\t";
//    std::cout << length(seqSet2[0]);
    std::cout << "\t\t\tdone\n" << std::flush;

//    for (auto& str : seqSet)
//        std::cout << str << "\n";
//    std::cout << std::flush;

// TODO(rrahn): Make object configurable.

    configureExec(options, seqSet1, seqSet2);
    return EXIT_SUCCESS;
}