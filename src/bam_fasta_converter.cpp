//#define SEQAN_ALIGN_SIMD_PROFILE 1
//#define SEQAN_ENABLE_TESTING 0
//#define SEQAN_ENABLE_DEBUG 0

//#define DP_PARALLEL_SHOW_PROGRESS
//#define DP_ALIGN_STATS

#include <cxxabi.h>
#include <array>
#include <string>
#include <future>

#ifdef DP_ALIGN_STATS
std::atomic<uint32_t> simdCounter;
std::atomic<uint32_t> serialCounter;
#endif

#include <iostream>

#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>


using namespace seqan;

struct options
{
    std::string bam_file;
    std::string ref_file;

    std::string reads_out;
    std::string subjects_out;
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
parseCommandLine(options & opt, int const argc, char* argv[])
{
    ArgumentParser parser("bam_fasta_converter");

    setShortDescription(parser, "Converts bam files into fasta files with a read file and query file");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);

    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM_FILE"));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "REF_FILE"));

    addOption(parser, seqan::ArgParseOption("r", "reads", "Output file containing the extracted reads in fasta format.", seqan::ArgParseArgument::OUTPUT_FILE, "READS"));
    setValidValues(parser, "r", "fa fasta");
    setDefaultValue(parser, "r", "reads.fa");

    addOption(parser, seqan::ArgParseOption("s", "subjects", "Output file containing the extracted subregions of the reference sequence in fasta format.", seqan::ArgParseArgument::OUTPUT_FILE, "SUBJECTS"));
    setValidValues(parser, "s", "fa fasta");
    setDefaultValue(parser, "s", "subject.fa");

    // Parse command line.
    if (parse(parser, argc, argv) != ArgumentParser::PARSE_OK)
        return ArgumentParser::PARSE_ERROR;

    getArgumentValue(opt.bam_file, parser, 0);
    getArgumentValue(opt.ref_file, parser, 1);

    getOptionValue(opt.reads_out, parser, "r");
    getOptionValue(opt.subjects_out, parser, "s");

    return ArgumentParser::PARSE_OK;
}

int main(int argc, char* argv[])
{
    options opt;

    if (parseCommandLine(opt, argc, argv) != ArgumentParser::PARSE_OK)
        return EXIT_FAILURE;

    try {
        SeqFileOut read_file{opt.reads_out.c_str()};  // Open read output file.
        SeqFileOut sbj_file{opt.subjects_out.c_str()};  // Open subject output file.

        SeqFileIn ref_file{opt.ref_file.c_str()};  // Open reference file.
        BamFileIn bam_file{opt.bam_file.c_str()};  // Open bam file.

        StringSet<CharString> ref_ids;
        StringSet<Dna5String> ref_seqs;

        BamHeader header;
        readHeader(header, bam_file);

        size_t numRecords = 100000;

        readRecords(ref_ids, ref_seqs, ref_file);
        while(!atEnd(bam_file))
        {
            std::cout << "Load next chunk\n" << std::flush;

            StringSet<BamAlignmentRecord> records;
            readRecords(records, bam_file, numRecords);

            auto extract_seq = [&](auto & align_row)
            {
                Dna5String seq;
                reserve(seq, length(align_row), Exact());
                for (auto it = begin(align_row); it != end(align_row); ++it)
                {
                    if (!isGap(it))
                        appendValue(seq, *it);
                }
                return seq;
            };

            size_t threadNum = std::thread::hardware_concurrency();

            std::vector<String<CharString>> reads_id_local(threadNum);
            std::vector<String<Dna5String>> reads_local(threadNum);
            std::vector<String<CharString>> subjects_id_local(threadNum);
            std::vector<String<Dna5String>> subjects_local(threadNum);

            for (unsigned  i = 0; i < threadNum; ++i)
            {
                reserve(reads_id_local[i],    numRecords);
                reserve(reads_local[i],       numRecords);
                reserve(subjects_id_local[i], numRecords);
                reserve(subjects_local[i],    numRecords);
            }

            std::cout << "Start extracting\n" << std::flush;
            SEQAN_OMP_PRAGMA(parallel for schedule(dynamic) num_threads(threadNum))
            for (int i = 0; i < length(records); ++i)
            // for (int job = 0; job < static_cast<int>(length(splitter)); ++job)
            {
                auto it = begin(records, Standard()) + i;
                size_t thread_id = omp_get_thread_num();
                if ((*it).rID == -1)
                {
                    continue;
                }
                Align<Dna5String, AnchorGaps<>> align;
                bamRecordToAlignment(align, ref_seqs[(*it).rID], *it);

//                std::cout << "Alignment: " << std::endl;
//                std::cout << align << std::endl;

                appendValue(reads_id_local[thread_id], (*it).qName, Generous());
                appendValue(reads_local[thread_id], extract_seq(row(align, 1)), Generous());
//                std::cout << "length(reads_local[thread_id]) = " << length(back(reads_local[thread_id])) << std::endl;

                CharString sbj_id{ref_ids[(*it).rID]};
                append(sbj_id, "_");
                append(sbj_id, (*it).qName);
                appendValue(subjects_id_local[thread_id], sbj_id, Generous());
                appendValue(subjects_local[thread_id], extract_seq(row(align, 0)), Generous());
//                printf("t_%d: read = %d - ref = %d\n", omp_get_thread_num(), length(back(reads_local[thread_id])),
//                length(back(subjects_local[thread_id])));
                //printf();
                //std::cout << "length(subjects_local[thread_id]) = " << length(back(subjects_local[thread_id])) << std::endl;
                // }
            }
            std::cout << "Dump to file\n\n" << std::flush;
            // This section must be serialized.
            for (unsigned  i = 0; i < threadNum; ++i)
            {
                writeRecords(read_file, reads_id_local[i], reads_local[i]);
                writeRecords(sbj_file, subjects_id_local[i], subjects_local[i]);
            }
        }
    } catch(seqan::ParseError & e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch(seqan::IOError & e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
