#include <benchmark/benchmark.h>

#include <utility>

#include <seqan/align.h>

#include <pairalign/benchmark/units.hpp>

static void bench(benchmark::State& state) {

    seqan::AlignConfig<true, true, true, true> alignConfig;

    seqan::Dna5String strH = "ATGT";
    seqan::DnaString strV = "ATAGAT";

    using alignment_instance_t = std::pair<seqan::Dna5String, seqan::Dna5String>;

    seqan::StringSet<alignment_instance_t> strings;
    seqan::appendValue(strings, std::pair{strH, strV});

    seqan::Score<int, seqan::Simple> scoringScheme(2, -1, -1);

    int32_t res = 0;

    for (auto _ : state) {
        res = seqan::globalAlignmentScore(strH, strV, scoringScheme, alignConfig, seqan::NeedlemanWunsch());
    }

    benchmark::DoNotOptimize(res);
    state.counters["CUPS"] = pairalign::units::cups(strings);
}
// Register the function as a benchmark
BENCHMARK(bench);
// Run the benchmark
BENCHMARK_MAIN();
