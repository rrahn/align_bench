#include "bench_fixture.hpp"

#include <seqan/align.h>

template <typename alphabet_t, typename score_t = int32_t>
class affine_bench_fixture : public bench_fixture<seqan::String<alphabet_t>> {

    using scoring_scheme_t = seqan::Score<score_t, seqan::Simple>;

    scoring_scheme_t _scheme{5, -4, -1, -11};

public:
    using score_type = score_t;
    scoring_scheme_t scoring_scheme() const noexcept {
        return _scheme;
    }
};
