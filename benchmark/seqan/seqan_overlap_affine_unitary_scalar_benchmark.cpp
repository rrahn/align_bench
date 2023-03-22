#include <benchmark/benchmark.h>

#include <algorithm>
#include <utility>

#include <seqan/align.h>
#include <seqan/align_parallel.h>

#include <pairalign/benchmark/units.hpp>
#include "fastq_reader.hpp"

static void bench(benchmark::State& state) {
    // ----------------------------------------------------------------------------
    // Prepare sequences
    // ----------------------------------------------------------------------------
    using sequence_type = seqan::String<seqan::Dna5>;
    using sequence_collection_type = seqan::StringSet<sequence_type>;
    sequence_collection_type seq1 = pairalign::input::load_fastq(DATADIR"sim_reads_n1K_rl150.fq");
    sequence_collection_type seq2{seq1};
    std::ranges::rotate(seq2, std::ranges::next(std::ranges::begin(seq2)));

    // ----------------------------------------------------------------------------
    // Configure alignment
    // ----------------------------------------------------------------------------
    using exec_policy_t = seqan::ExecutionPolicy<seqan::Serial, seqan::Serial>;
    using end_gap_config_t = seqan::AlignConfig<true, true, true, true>;
    seqan::Score<int, seqan::Simple> scoringScheme(5, -4, -1, -11);

    // ----------------------------------------------------------------------------
    // Run alignment
    // ----------------------------------------------------------------------------
    seqan::String<int32_t> res{};
    for (auto _ : state) {
        res = seqan::globalAlignmentScore(exec_policy_t{}, seq1, seq2, scoringScheme, end_gap_config_t{});
    }

    // ----------------------------------------------------------------------------
    // Prepare output
    // ----------------------------------------------------------------------------
    using alignment_instance_t = std::pair<sequence_type, sequence_type>;
    std::vector<alignment_instance_t> sequences{};
    sequences.resize(std::ranges::size(seq1));
    for (int32_t i = 0; i < std::ranges::ssize(sequences); ++i) {
        sequences[i] = std::pair{seq1[i], seq2[i]};
    }

    benchmark::DoNotOptimize(res);
    state.counters["CUPS"] = pairalign::units::cups(sequences);
}
// Register the function as a benchmark
BENCHMARK(bench);
// Run the benchmark
BENCHMARK_MAIN();
