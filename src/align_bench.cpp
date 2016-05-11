#include <iostream>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/stream.h>
#include <seqan/align.h>
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

//    addOption(parser, seqan::ArgParseOption(
//                                            "i", "period", "Period to use for the index.",
//                                            seqan::ArgParseArgument::INTEGER, "INT"));
//    addOption(parser, seqan::ArgParseOption(
//                                            "U", "uppercase", "Select to-uppercase as operation."));

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

//    getArgumentValue(options.inputFileOne, parser, 0);
//    getArgumentValue(options.inputFileTwo, parser, 1);

    return ArgumentParser::PARSE_OK;
}

//template <typename TString, typename TSetSpec>
//inline void loadSequence(StringSet<TString, TSetSpec> & set,
//                         char const * fileName)
//{
//    SeqFileIn
//}

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

    SequenceGenerator<AminoAcid> gen;
    gen.setNumber(200000);
    gen.setMinLength(250);
    gen.setMaxLength(250);

    std::cout << "Start generating sequences..." << std::flush;
    auto seqSet1 = gen.generate();
    auto seqSet2 = gen.generate();
    std::cout << "\t done\n" << std::flush;

//    for (auto& str : seqSet)
//        std::cout << str << "\n";
//    std::cout << std::flush;

// TODO(rrahn): Make object configurable.
    BenchmarkExecutor device;

    device.runGlobalAlignment(seqSet1, seqSet2, 4, seqan::Parallel());
    device.printProfile(std::cout);

//    device.runGlobalAlignment(seqSet1, seqSet2, seqan::Serial());
//    device.printProfile(std::cout);

    return EXIT_SUCCESS;
}