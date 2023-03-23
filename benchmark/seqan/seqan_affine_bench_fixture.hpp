#include <string_view>

#include <seqan/align.h>

#include <pairalign/benchmark/data_sources.hpp>

#include "seqan_bench_fixture.hpp"

template <auto * data, typename alphabet_t, typename score_t = int32_t>
class seqan_affine_bench_fixture : public seqan_bench_fixture<alphabet_t, data> {

    using scoring_scheme_t = seqan::Score<score_t, seqan::Simple>;
    using matrix_scheme_t = seqan::Score<score_t, seqan::ScoreMatrix<seqan::AminoAcid, seqan::ScoreSpecBlosum62>>;

    scoring_scheme_t _scheme{5, -4, -1, -11};

public:
    using score_type = score_t;
    scoring_scheme_t scoring_scheme() const noexcept {
        return _scheme;
    }

    matrix_scheme_t matrix_scheme() const noexcept {
        return matrix_scheme_t{-1, -11};
    }
};
