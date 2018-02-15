// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef SEQUENCE_GENERATOR_HPP_
#define SEQUENCE_GENERATOR_HPP_

#include <random>

namespace seqan
{

template <typename TAlphabet_ = Dna>
class SequenceGenerator
{
public:

    using TAlphabet = TAlphabet_;

    // ----------------------------------------------------------------------------
    // Constructor.
    // ----------------------------------------------------------------------------

    // Default c'tor. - Needed?
    SequenceGenerator() : mSeed(-1), mMinLength(100), mMaxLength(1000), mNum(1),
                          mDist(DistributionFunction::UNIFORM_DISTRIBUTION)
    {
        mRng.seed(mSeed);
    }

    // ----------------------------------------------------------------------------
    // Member functions.
    // ----------------------------------------------------------------------------

    inline void setMinLength(unsigned const size)
    {
        mMinLength = size;
    }

    inline void setMaxLength(unsigned const size)
    {
        mMaxLength = size;
    }

    inline void setNumber(unsigned const size)
    {

        mNum = size;
    }

    inline void setDistribution(DistributionFunction dist)
    {

        mDist = dist;
    }

    inline StringSet<String<TAlphabet> > generate();

private:

    template <typename TPdf>
    inline auto generateImpl(TPdf & pdf);


    // ----------------------------------------------------------------------------
    // Internal member variables.
    // ----------------------------------------------------------------------------

    unsigned                mSeed;
    unsigned                mMinLength;
    unsigned                mMaxLength;
    unsigned                mNum;
    DistributionFunction    mDist;
    std::mt19937            mRng;
};

// ----------------------------------------------------------------------------
// Function generate()
// ----------------------------------------------------------------------------

template <typename TAlphabet>
template <typename TPdf>
inline auto
SequenceGenerator<TAlphabet>::generateImpl(TPdf & pdf)
{
    StringSet<String<TAlphabet> > set;
    resize(set, mNum, Exact());

    std::uniform_int_distribution<> charDis(0, ValueSize<TAlphabet>::VALUE - 1);

    for (auto& str : set)
    {
        resize(str, static_cast<size_t>(pdf(mRng)), Exact());
        for (auto& val : str)
            val = static_cast<TAlphabet>(charDis(this->mRng));
    }
    return set;
}

template <typename TAlphabet>
inline StringSet<String<TAlphabet> >
SequenceGenerator<TAlphabet>::generate()
{
    // Sanity check.
    if (this->mMaxLength < this->mMinLength)
        this->mMaxLength = this->mMinLength;

    if (mDist == DistributionFunction::UNIFORM_DISTRIBUTION)
    {
        std::uniform_int_distribution<> lenDis(this->mMinLength, this->mMaxLength);
        return generateImpl(lenDis);
    }
    else
    {
        std::normal_distribution<> lenDis((mMinLength + mMaxLength)/2, 1.0);
        return generateImpl(lenDis);
    }
}

}

#endif  // #ifndef SEQUENCE_GENERATOR_HPP_
