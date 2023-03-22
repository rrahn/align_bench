
#include <seqan/align.h>

#include "seqan_bench_fixture.hpp"
template <typename alphabet_t, typename score_t = int32_t>
class seqan_affine_bench_fixture : public seqan_bench_fixture<alphabet_t> {

    using scoring_scheme_t = seqan::Score<score_t, seqan::Simple>;

    scoring_scheme_t _scheme{5, -4, -1, -11};

public:
    using score_type = score_t;
    scoring_scheme_t scoring_scheme() const noexcept {
        return _scheme;
    }
};
