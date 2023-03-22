#include "pa_bench_fixture.hpp"

#include <pairwise_aligner/configuration/gap_model_affine.hpp>

template <typename alphabet_t, typename score_t = int32_t>
class pa_affine_bench_fixture : public pa_bench_fixture<alphabet_t> {

public:
    using score_type = score_t;

    score_type match_cost() const noexcept {
        return 5;
    }

    score_type mismatch_cost() const noexcept {
        return -4;
    }

    constexpr auto gap_cost() const noexcept {
        return seqan::pairwise_aligner::cfg::gap_model_affine(-10, -1);
    }
};
