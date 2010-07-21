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

/** \addtogroup MachineLearning Machine Learning
**/
//@{


/**\brief Options object for the Sampler class.
 
  <b>usage:</b>
 
  \code
  SamplerOptions opt =  SamplerOptions()
                               .withReplacement()
                               .sampleProportion(0.5);
  \endcode
 
  Note that the return value of all methods is <tt>*this</tt> which makes
  concatenating of options as above possible.
*/
class SamplerOptions
{
  public:

    double sample_proportion;
    unsigned int sample_size;
    bool   sample_with_replacement;
    bool   stratified_sampling;
    
    SamplerOptions()
    : sample_proportion(1.0),
      sample_size(0), 
      sample_with_replacement(true),
      stratified_sampling(false)
    {}

        /**\brief Sample from training population with replacement.
         *
         * <br> Default: true
         */
    SamplerOptions& withReplacement(bool in = true)
    {
        sample_with_replacement = in;
        return *this;
    }

        /**\brief Sample from training population without replacement.
         *
         * <br> Default (if you don't call this function): false
         */
    SamplerOptions& withoutReplacement(bool in = true)
    {
        sample_with_replacement = !in;
        return *this;
    }

        /**\brief Draw the given number of samples.
         * If stratifiedSampling is true, the \a size is equally distributed
         * accross all strata (e.g. <tt>size / strataCount</tt> samples are taken 
         * from each stratum, subject to rounding).
         *
         * <br> Default: 0 (i.e. determine the count by means of sampleProportion())
         */
    SamplerOptions& sampleSize(unsigned int size)
    {
        sample_size = size;
        return *this;
    }


        /**\brief Determine the number of samples to draw as a proportion of the total
         * number. That is, we draw <tt>count = totalCount * proportion</tt> samples. 
         * This option is overridden when an absolute count is specified by sampleSize().
         * 
         * If stratifiedSampling is true, the count is equally distributed
         * accross all strata (e.g. <tt>totalCount * proportion / strataCount</tt> samples are taken 
         * from each stratum).
         *
         * <br> Default: 1.0
         */
    SamplerOptions& sampleProportion(double proportion)
    {
        vigra_precondition(proportion >= 0.0,
               "SamplerOptions::sampleProportion(): argument must not be negative.");
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
    SamplerOptions& stratified(bool in = true)
    {
        stratified_sampling = in;
        return *this;
    }
};

/************************************************************/
/*                                                          */
/*                        Sampler                           */
/*                                                          */
/************************************************************/

/** \brief Create random samples from a sequence of indices.

    Selecting data items at random is a basic task of machine learning,
    for example in boostrapping, RandomForest training, and cross validation.
    This class implements various ways to select random samples via their indices. 
    Indices are assumed to be consecutive in
    the range <tt>0 &lt;= index &lt; total_sample_count</tt>.
    
    The class always contains a current sample which can be accessed by 
    the index operator or by the function sampledIndices(). The indices
    that are not in the current sample (out-of-bag indices) can be accessed
    via the function oobIndices().
    
    The sampling method (with/without replacement, stratified or not) and the
    number of samples to draw are determined by the option object 
    SamplerOptions.
    
    <b>Usage:</b>
    
    <b>\#include</b> \<<a href="index__sampling_8hxx-source.html">vigra/index_sampling.hxx</a>\><br>
    Namespace: vigra
    
    Create a Sampler with default options, i.e. sample as many indices as there 
    are data elements, with replacement. On average, the sample will contain 
    <tt>0.63*totalCount</tt> distinct indices.
    
    \code
    int totalCount = 10000;   // total number of data elements
    int numberOfSamples = 20; // repeat experiment 20 times 
    Sampler<> sampler(totalCount);
    for(int k=0; k<numberOfSamples; ++k)
    {
        // process current sample
        for(int i=0; i<sampler.sampleSize(); ++i)
        {
            int currentIndex = sampler[i];
            processData(data[currentIndex]);
        }
        // create next sample
        sampler.sample();
    }
    \endcode
    
    Create a Sampler for stratified sampling, without replacement.
    
    \code
    // prepare the strata (i.e. specify which stratum each element belongs to)
    int stratumSize1 = 2000, stratumSize2 = 8000,
        totalCount = stratumSize1 + stratumSize2;
    ArrayVerctor<int> strata(totalCount);
    for(int i=0; i<stratumSize1; ++i)
        strata[i] = 1;
    for(int i=stratumSize1; i<stratumSize2; ++i)
        strata[i] = 2;
        
    int sampleSize = 200; // i.e. sample 100 elements from each of the two strata
    int numberOfSamples = 20; // repeat experiment 20 times 
    Sampler<> stratifiedSampler(strata.begin(), strata.end(),
                     SamplerOptions().withoutReplacement().stratified().sampleSize(sampleSize));

    for(int k=0; k<numberOfSamples; ++k)
    {
        // process current sample
        for(int i=0; i<sampler.sampleSize(); ++i)
        {
            int currentIndex = sampler[i];
            processData(data[currentIndex]);
        }
        // create next sample
        sampler.sample();
    }
    \endcode
*/
template<class Random = MersenneTwister >
class Sampler
{
  public:
        /** Internal type of the indices.
            Currently, 64-bit indices are not supported because this
            requires extension of the random number generator classes.
        */
    typedef Int32                               IndexType;
    
    typedef ArrayVector     <IndexType>  IndexArrayType;
    
        /** Type of the array view object that is returned by 
            sampledIndices() and oobIndices().
        */
    typedef ArrayVectorView <IndexType>  IndexArrayViewType;

  private:
    typedef std::map<IndexType, IndexArrayType> StrataIndicesType;
    typedef std::map<IndexType, int> StrataSizesType;
    typedef ArrayVector     <bool>       IsUsedArrayType;
    typedef ArrayVectorView <bool>       IsUsedArrayViewType;
    
    static const int oobInvalid = -1;

    int total_count_, sample_size_;
	mutable int current_oob_count_;
    StrataIndicesType     strata_indices_;
    StrataSizesType       strata_sample_size_;
    IndexArrayType        current_sample_;
    mutable IndexArrayType        current_oob_sample_;
    IsUsedArrayType       is_used_;
    Random random_;
    SamplerOptions options_;

    void initStrataCount()
    {
        // compute how many samples to take from each stratum
        // (may be unequal if sample_size_ is not a multiple of strataCount())
        int strata_sample_size = (int)std::ceil(double(sample_size_) / strataCount());
        int strata_total_count = strata_sample_size * strataCount();

        for(StrataIndicesType::iterator i = strata_indices_.begin(); 
             i != strata_indices_.end(); ++i)
        {
            if(strata_total_count > sample_size_)
            {
                strata_sample_size_[i->first] = strata_sample_size - 1;
                --strata_total_count;
            }
            else
            {
                strata_sample_size_[i->first] = strata_sample_size;
            }
        }
    }

  public:
    
        /** Create a sampler for \a totalCount data objects.
            
            In each invocation of <tt>sample()</tt> below, it will sample
            indices according to the options passed. If no options are given, 
            <tt>totalCount</tt> indices will be drawn with replacement.
        */
    Sampler(UInt32 totalCount, SamplerOptions const & opt = SamplerOptions(), 
            Random const & rnd = Random(RandomSeed))
    : total_count_(totalCount),
      sample_size_(opt.sample_size == 0
                         ? (int)(std::ceil(total_count_ * opt.sample_proportion))
                         : opt.sample_size),
      current_oob_count_(oobInvalid),
      current_sample_(sample_size_),
      current_oob_sample_(total_count_),
      is_used_(total_count_),
      random_(rnd),
      options_(opt)
    {
        vigra_precondition(opt.sample_with_replacement || sample_size_ <= total_count_,
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
            Equally many samples (subject to rounding) will be drawn from each stratum, 
            unless the option object explicitly requests unstratified sampling, 
            in which case the strata are ignored.
        */
    template <class Iterator>
    Sampler(Iterator strataBegin, Iterator strataEnd, SamplerOptions const & opt = SamplerOptions(), 
            Random const & rnd = Random(RandomSeed))
    : total_count_(strataEnd - strataBegin),
      sample_size_(opt.sample_size == 0
                         ? (int)(std::ceil(total_count_ * opt.sample_proportion))
                         : opt.sample_size),
      current_oob_count_(oobInvalid),
      current_sample_(sample_size_),
      current_oob_sample_(total_count_),
      is_used_(total_count_),
      random_(rnd),
      options_(opt)
    {
        vigra_precondition(opt.sample_with_replacement || sample_size_ <= total_count_,
          "Sampler(): Cannot draw without replacement when data size is smaller than sample count.");
          
        // copy the strata indices
        for(int i = 0; strataBegin != strataEnd; ++i, ++strataBegin)
        {
            strata_indices_[*strataBegin].push_back(i);
        }
        vigra_precondition(sample_size_ >= (int)strata_indices_.size(),
            "Sampler(): Requested sample count must be at least as large as the number of strata.");

        initStrataCount();
        sample();
    }

        /** Return the k-th index in the current sample.
         */
    IndexType operator[](int k) const
    {
        return current_sample_[k];
    }

        /** Create a new sample.
         */
    void sample();

        /** The total number of data elements.
         */
    int totalCount() const
    {
        return total_count_;
    }

        /** The number of data elements that have been sampled.
         */
    int sampleSize() const
    {
        return sample_size_;
    }

        /** Same as sampleSize().
         */
    int size() const
    {
        return sample_size_;
    }

        /** The number of strata to be used.
            Will be 1 if no strata are given. Will be ognored when
            stratifiedSampling() is false.
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


template<class Random>
void Sampler<Random>::sample()
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
            for(int i = 0; i < (int)strata_sample_size_[iter->first]; ++i, ++j)
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
            for(int i = 0; i < (int)strata_sample_size_[iter->first]; ++i, ++j)
            {
                std::swap(iter->second[i], iter->second[i+ random_.uniformInt(stratum_size - i)]);
                current_sample_[j] = iter->second[i];
                is_used_[current_sample_[j]] = true;
            }
        }
    }
}


//@}

} // namespace vigra

#endif /*VIGRA_INDEX_SAMPLING_HXX*/
