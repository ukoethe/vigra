/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by  Ullrich Koethe and Rahul Nair         */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_INDEX_SAMPLING_HXX
#define VIGRA_INDEX_SAMPLING_HXX

#include "array_vector.hxx"
#include "random.hxx"
#include <map>
#include <memory>
#include <cmath>

namespace vigra
{

class SamplingOptions
{
  public:

    double sample_proportion;
    unsigned int sample_count;
    bool   sample_with_replacement;
    bool   stratified_sampling;
    
    SamplingOptions()
    : sample_proportion(1.0),
      sample_count(0), 
      sample_with_replacement(true),
      stratified_sampling(false)
    {}

        /**\brief Sample from training population with replacement.
         *
         * <br> Default: true
         */
    SamplingOptions& withReplacement(bool in = true)
    {
        sample_with_replacement = in;
        return *this;
    }

        /**\brief Sample from training population without replacement.
         *
         * <br> Default (if you don't call this function): false
         */
    SamplingOptions& withoutReplacement(bool in = true)
    {
        sample_with_replacement = !in;
        return *this;
    }

        /**\brief Draw the given number of samples.
         * If stratifiedSampling is true, the count is equally distributed
         * accross all strata (e.g. <tt>count / strataCount</tt> samples are taken 
         * from each stratum).
         *
         * <br> Default: 0 (i.e. determine the count by means of sampleProportion())
         */
    SamplingOptions& sampleCount(unsigned int count)
    {
        sample_count = count;
        return *this;
    }


        /**\brief Determine the number of samples to draw as a proportion of the total
         * number. That is, we draw <tt>count = totalCount * proportion</tt> samples. 
         * This option is overridden when an absolute count is specified by sampleCount().
         * 
         * If stratifiedSampling is true, the count is equally distributed
         * accross all strata (e.g. <tt>totalCount * proportion / strataCount</tt> samples are taken 
         * from each stratum).
         *
         * <br> Default: 1.0
         */
    SamplingOptions& sampleProportion(double proportion)
    {
        vigra_precondition(proportion >= 0.0,
               "SamplingOptions::sampleProportion(): argument must not be negative.");
        sample_proportion = proportion;
        return *this;
    }

        /**\brief Draw equally many samples from each "stratum". 
         *  A stratum is a group of like entities, e.g. pixels belonging 
         *  to the same object class. This is useful to create balanced samples  
         *  when the class probabilities are very unbalanced (e.g.
         *  when there are many background and few foreground pixels).
         *  Stratified sampling thus avoids that a trained classifier is biased 
         *  towards the majority class. 
         *
         * <br> Default (if you don't call this function): false
         */
    SamplingOptions& stratified(bool in = true)
    {
        stratified_sampling = in;
        return *this;
    }
};

template<class Random = MersenneTwister >
class Sampler
{
  public:
    typedef Int32                               IndexType;
    typedef ArrayVector     <IndexType>  IndexArrayType;
    typedef ArrayVectorView <IndexType>  IndexArrayViewType;

  private:
    typedef std::map<IndexType, IndexArrayType> StrataIndicesType;
    typedef std::map<IndexType, int> StrataSizesType;
    typedef ArrayVector     <bool>       IsUsedArrayType;
    typedef ArrayVectorView <bool>       IsUsedArrayViewType;
    
    static const int oobInvalid = -1;

    int total_count_, sample_count_;
	mutable int current_oob_count_;
    StrataIndicesType     strata_indices_;
    StrataSizesType       strata_sample_count_;
    IndexArrayType        current_sample_;
    mutable IndexArrayType        current_oob_sample_;
    IsUsedArrayType       is_used_;
    Random random_;
    SamplingOptions options_;

    void initStrataCount()
    {
        // compute how many samples to take from each stratum
        // (may be unequal if sample_count_ is not a multiple of strataCount())
        int strata_sample_count = (int)std::ceil(double(sample_count_) / strataCount());
        int strata_total_count = strata_sample_count * strataCount();

        for(StrataIndicesType::iterator i = strata_indices_.begin(); 
             i != strata_indices_.end(); ++i)
        {
            if(strata_total_count > sample_count_)
            {
                strata_sample_count_[i->first] = strata_sample_count - 1;
                --strata_total_count;
            }
            else
            {
                strata_sample_count_[i->first] = strata_sample_count;
            }
        }
    }

  public:
    
        /** Create a sampler for \a totalCount data objects.
            
            In each invocation of <tt>sample()</tt> below, it will sample
            indices according to the options passed. If no options are given, 
            <tt>totalCount</tt> indices will be drawn with replacement.
        */
    Sampler(UInt32 totalCount, SamplingOptions const & opt = SamplingOptions(), 
            Random const & rnd = Random(RandomSeed))
    : total_count_(totalCount),
      sample_count_(opt.sample_count == 0
                         ? (int)(std::ceil(total_count_ * opt.sample_proportion))
                         : opt.sample_count),
      current_oob_count_(oobInvalid),
      current_sample_(sample_count_),
      current_oob_sample_(total_count_),
      is_used_(total_count_),
      random_(rnd),
      options_(opt)
    {
        vigra_precondition(opt.sample_with_replacement || sample_count_ <= total_count_,
          "Sampler(): Cannot draw without replacement when data size is smaller than sample count.");
          
        vigra_precondition(!opt.stratified_sampling,
          "Sampler(): Stratified sampling requested, but no strata given.");
          
        // initialize a single stratum containing all data
        strata_indices_[0].resize(total_count_);
        for(int i=0; i<total_count_; ++i)
            strata_indices_[0][i] = i;

        initStrataCount();
        sample();
    }
    
        /** Create a sampler for stratified sampling.
            
            <tt>strataBegin</tt> and <tt>strataEnd</tt> must refer to a sequence 
            which specifies for each sample the stratum it belongs to. The
            total number of data objects will be set to <tt>strataEnd - strataBegin</tt>.
            Equally many samples (+- 1) will be drawn from each stratum, unless the 
            option object explicitly requests unstratified sampling, in which case the 
            strata are ignored.
        */
    template <class Iterator>
    Sampler(Iterator strataBegin, Iterator strataEnd, SamplingOptions const & opt = SamplingOptions(), 
            Random const & rnd = Random(RandomSeed))
    : total_count_(strataEnd - strataBegin),
      sample_count_(opt.sample_count == 0
                         ? (int)(std::ceil(total_count_ * opt.sample_proportion))
                         : opt.sample_count),
      current_oob_count_(oobInvalid),
      current_sample_(sample_count_),
      current_oob_sample_(total_count_),
      is_used_(total_count_),
      random_(rnd),
      options_(opt)
    {
        vigra_precondition(opt.sample_with_replacement || sample_count_ <= total_count_,
          "Sampler(): Cannot draw without replacement when data size is smaller than sample count.");
          
        // copy the strata indices
        for(int i = 0; strataBegin != strataEnd; ++i, ++strataBegin)
        {
            strata_indices_[*strataBegin].push_back(i);
        }
        vigra_precondition(sample_count_ >= (int)strata_indices_.size(),
            "Sampler(): Requested sample count must be at least as large as the number of strata.");

        initStrataCount();
        sample();
    }

        /** Return the k-th sampled index.
         */
    IndexType const & operator[](int k)
    {
        return current_sample_[k];
    }

        /** Create a new sample.
         */
    void sample()
    {
        current_oob_count_ = oobInvalid;
        is_used_.init(false);
        
        if(options_.sample_with_replacement)
        {
            //Go thru all strata
            int j = 0;
            StrataIndicesType::iterator iter;
            for(iter = strata_indices_.begin(); iter != strata_indices_.end(); ++iter)
            {
                // do sampling with replacement in each strata and copy data.
                int stratum_size = iter->second.size();
                for(int i = 0; i < (int)strata_sample_count_[iter->first]; ++i, ++j)
                {
                    current_sample_[j] = iter->second[random_.uniformInt(stratum_size)];
                    is_used_[current_sample_[j]] = true;
                }
            }
        }
        else
        {
            //Go thru all strata
            int j = 0;
            StrataIndicesType::iterator iter;
            for(iter = strata_indices_.begin(); iter != strata_indices_.end(); ++iter)
            {
                // do sampling without replacement in each strata and copy data.
                int stratum_size = iter->second.size();
                for(int i = 0; i < (int)strata_sample_count_[iter->first]; ++i, ++j)
                {
                    std::swap(iter->second[i], iter->second[i+ random_.uniformInt(stratum_size - i)]);
                    current_sample_[j] = iter->second[i];
                    is_used_[current_sample_[j]] = true;
                }
            }
        }
    }

        /** The total number of data elements.
         */
    int totalCount() const
    {
        return total_count_;
    }

        /** The number of data elements to be sampled.
         */
    int sampleCount() const
    {
        return sample_count_;
    }

        /** The number of strata to be used.
            will be 1 if no strata are given.
         */
    int strataCount() const
    {
        return strata_indices_.size();
    }

        /** Whether to use stratified sampling.
            (If this is false, strata will be ignored even if present.)
         */
    bool stratifiedSampling() const
    {
        return options_.stratified_sampling;
    }
    
        /** Whether sampling should be performed with replacement.
         */
    bool withReplacement() const
    {
        return options_.sample_with_replacement;
    }
    
        /** Return an array view containing the indices in the current sample.
         */
    IndexArrayViewType sampledIndices() const
    {
        return current_sample_;
    }
    
        /** Return an array view containing the out-of-bag indices.
            (i.e. the indices that are not in the current sample)
         */
    IndexArrayViewType oobIndices() const
    {
        if(current_oob_count_ == oobInvalid)
        {
            current_oob_count_ = 0;
            for(int i = 0; i<total_count_; ++i)
            {
                if(!is_used_[i])
                {
                    current_oob_sample_[current_oob_count_] = i;
                    ++current_oob_count_;
                }
            }
        }
        return current_oob_sample_.subarray(0, current_oob_count_);
    }
};

} // namespace vigra

#endif /*VIGRA_INDEX_SAMPLING_HXX*/
