//#define SEQAN_ALIGN_SIMD_PROFILE 1
//#define SEQAN_ENABLE_TESTING 0
//#define SEQAN_ENABLE_DEBUG 0

//#define DP_PARALLEL_SHOW_PROGRESS
//#define DP_ALIGN_STATS

#ifdef DP_ALIGN_STATS
std::atomic<uint32_t> simdCounter;
std::atomic<uint32_t> serialCounter;
#endif

#include "align_bench_parser.hpp"
#include "align_bench_configure.hpp"
#include "align_bench_options.hpp"

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
    ArgumentParser parser("align_bench_par");

    setShortDescription(parser, "Alignment Benchmark Tool");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    setup_parser(parser);

    addOption(parser, seqan::ArgParseOption("t", "threads", "Number of threads", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "t", toString(std::thread::hardware_concurrency()));

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

    get_arguments(options, parser);

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

    getOptionValue(options.threadCount, parser, "t");

    return ArgumentParser::PARSE_OK;
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
    options.stats.threads = options.threadCount;
    if (options.simd)
    {
        options.stats.execPolicy = "parallel_vec";
        options.stats.vectorLength = SEQAN_SIZEOF_MAX_VECTOR / static_cast<unsigned>(options.simdWidth);
        seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial> exec_policy;
        setNumThreads(exec_policy, options.threadCount);
        configureAlpha(options, exec_policy);
    }
    else
    {
        options.stats.execPolicy = "parallel";
        seqan::ExecutionPolicy<seqan::Parallel, seqan::Serial> exec_policy;
        setNumThreads(exec_policy, options.threadCount);
        configureAlpha(options, exec_policy);
    }

    options.stats.state = "done";
    options.stats.writeHeader(std::cout);
    options.stats.writeStats(std::cout);

    return EXIT_SUCCESS;
}
