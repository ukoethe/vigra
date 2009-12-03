/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by  Ullrich Koethe and Rahul Nair         */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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
#ifndef RN_SAMPLER_HXX
#define RN_SAMPLER_HXX

#include <vigra/array_vector.hxx>
#include <vigra/random.hxx>
#include <map>
#include <math.h>
#include <memory>
namespace vigra
{

class SamplingOptions
{
    public:
    typedef std::auto_ptr<vigra::ArrayVectorView<int> > strata_ptr_type;
    strata_ptr_type strata;
    bool sample_with_replacement;
    bool sample_stratified;
    bool sample_classes_individually;
    bool use_internal_mem;
    ArrayVectorView<Int32> mem;


    SamplingOptions(): sample_with_replacement(true),
                        sample_stratified(false), sample_classes_individually(false),
                        use_internal_mem(true)
    { this->strata = strata_ptr_type(new vigra::ArrayVectorView<int>());
    }

    SamplingOptions( SamplingOptions const& rhs): sample_with_replacement(rhs.sample_with_replacement),
                        sample_stratified(rhs.sample_stratified), sample_classes_individually(rhs.sample_classes_individually)
    {
        this->strata = strata_ptr_type(new vigra::ArrayVectorView<int>(*(rhs.strata)));
    }

    void operator=(SamplingOptions const& rhs)
    {
        this->sample_with_replacement = rhs.sample_with_replacement;
        this->sample_stratified = rhs.sample_stratified;
        this->sample_classes_individually = rhs.sample_classes_individually;
        this->strata = strata_ptr_type(new vigra::ArrayVectorView<int>(*(rhs.strata)));
        this->use_internal_mem = rhs.use_internal_mem;
        this->mem = rhs.mem;
    }
    SamplingOptions& sampleWithReplacement(bool in =true)
    {
        sample_with_replacement = in;
        return *this;
    }

    SamplingOptions& useExternalMemory(ArrayVectorView<Int32> & mem_)
    {
        this->use_internal_mem = false;
        this->mem = mem_;
        return *this;
    }


    SamplingOptions& sampleWithoutReplacement()
    {
        sample_with_replacement = false;
        return *this;
    }


    SamplingOptions& sampleStratified(vigra::ArrayVectorView<int> in)
    {

        strata = strata_ptr_type(new vigra::ArrayVectorView<int>(in));
        sample_stratified = true;
        sample_classes_individually = false;
        return *this;
    }

    SamplingOptions& sampleClassesIndividually(vigra::ArrayVectorView<int> in)
    {
        strata = strata_ptr_type(new vigra::ArrayVectorView<int>(in));
        sample_stratified = false;
        sample_classes_individually = true;
        return *this;
    }

    SamplingOptions& resetStrata()
    {
        strata.reset(new vigra::ArrayVectorView<int>());
        sample_stratified = false;
        sample_classes_individually = false;
        return *this;
    }
};








template<class Random =UniformIntRandomFunctor<RandomTT800> >
class Sampler
{
    public:
    typedef Int32                               IndexType;
    typedef vigra::ArrayVector     <IndexType>  IndexArrayType;
    typedef vigra::ArrayVectorView <IndexType>  IndexArrayViewType;
    typedef vigra::ArrayVector     <bool>       IsUsedArrayType;
    typedef vigra::ArrayVectorView <bool>       IsUsedArrayViewType;

    Random  randint;

    private:
        bool unused_indices_set;
        typedef std::map<IndexType, vigra::ArrayVector< IndexType> > StrataIndicesType;
        typedef std::map<IndexType, int> StrataSizesType;
        StrataIndicesType     strata_indices_;
        StrataSizesType       strata_sizes_;
        IndexArrayType        internal_used_indices_;
        IndexArrayViewType    used_indices_;
        IndexArrayType        unused_indices_;
        IsUsedArrayType         is_used_;

        void (Sampler::*sampling_func_ptr_)(void);

        void sample_stratified_w_rep()
        {

            is_used_.init(false);

            //Go thru each strata
            StrataIndicesType::iterator iter;
            int jj = 0;
            for(iter = strata_indices_.begin(); iter != strata_indices_.end(); ++iter)
            {
                // do sampling with replacement in each strata and copy
                // data.
                int sze = iter->second.size();
                for(int ii = 0; ii < (int)strata_sizes_[iter->first]; ++ii)
                {
                    used_indices_[jj] = iter->second[randint(sze)];
                    is_used_[used_indices_[jj] ] = 1;
                    jj++;
                }

            }
        }
        void sample_stratified_wo_rep()
        {


            // reset is_used
            is_used_.init(false);

            //Go thru each strata
            StrataIndicesType::iterator iter;
            int jj = 0;
            for(iter = strata_indices_.begin(); iter != strata_indices_.end(); ++iter)
            {
                // do sampling without replacement in each strata and copy
                // data.
                for(int ii = 0; ii < (int)strata_sizes_[iter->first]; ++ii)
                {


                    std::swap(iter->second[ii], iter->second[ii+ randint(iter->second.size() - ii)]);
                    used_indices_[jj] = iter->second[ii];

                    is_used_[used_indices_[jj] ] = 1;
                    jj++;
                }


            }
        }
        void sample_w_rep()
        {

            is_used_.init(false);
            for(int ii = 0; ii < num_of_samples; ++ii)
            {
                used_indices_[ii] = randint(max_index);
                is_used_[used_indices_[ii] ] = true;
            }
        }
        void sample_wo_rep()
        {

            is_used_.init(false);
            for(int ii = 0; ii < num_of_samples; ++ii)
            {
                std::swap(used_indices_[ii], used_indices_[ii+ randint(max_index - ii)]);
                is_used_[used_indices_[ii] ] = true;
            }
        }

    public:

    SamplingOptions options_;

    IndexType const & operator[](int in)
    {
        return used_indices_[in];
    }

    IndexArrayViewType used_indices()
    {
        return used_indices_.subarray(0, num_of_samples);
    }

    IndexArrayViewType unused_indices()
    {
        if(unused_indices_set)
        {
            return unused_indices_;
        }
        else
        {
            unused_indices_.clear();
            for(int ii = 0; ii < (int)is_used().size(); ++ii)
            {
                if(is_used_[ii])
                    unused_indices_.push_back(ii);
            }
            unused_indices_set = true;
        }
        return unused_indices_;
    }

    IndexArrayViewType used_indices_volatile()
    {
        return IndexArrayViewType(used_indices_);
    }

    IsUsedArrayViewType is_used()
    {
        return IsUsedArrayViewType(is_used_);
    }

    int num_of_samples;
    int max_index;

    inline void sample(  )
    {
        unused_indices_set = false;
        (this->*sampling_func_ptr_)();
    }

    int numOfSamples()
    {
        return num_of_samples;
    }

    int numOfSamples(int n)
    {
        num_of_samples = n;
        if(!options_.sample_with_replacement && ((*options_.strata).data() != 0) && options_.use_internal_mem)
        {
            internal_used_indices_.resize(num_of_samples);
            used_indices_ = internal_used_indices_;
        }
        else if (!options_.use_internal_mem)
            vigra_fail("Sampler::numOfSamples() : Sampler is using external memory - no resize possible");
        return num_of_samples;
    }

    void init(  int numOfSamples ,
                int maxIndex     ,
                SamplingOptions const & opt)
    {
        unused_indices_set = false;
        strata_indices_.clear();
        strata_sizes_.clear();
        max_index = maxIndex;
        num_of_samples = numOfSamples;
        options_ = opt;
        if(!options_.use_internal_mem)
        {
            used_indices_ = options_.mem;
        }
        //Allocate memory for the used/unused vector.

        is_used_.resize(max_index);



        if((*options_.strata).data() != 0)
        {


            // Set up the data struct used.
            for(int ii = 0; ii < maxIndex; ++ii)
            {

                strata_indices_[(*(options_.strata))[ii]].push_back(ii);
            }

            // Set up the strata size
            if(options_.sample_classes_individually)
            {
                StrataIndicesType::difference_type total_size = 0;
                StrataIndicesType::iterator iter;
                for(iter = strata_indices_.begin(); iter != strata_indices_.end(); ++iter)
                {
                    // Set the size of the strata to be sampled to the fixed proportion of num_of_samples
                    // This value is the same for all strata
                    strata_sizes_[iter->first] = (int)ceil(double(num_of_samples) / strata_indices_.size());
                    total_size += strata_sizes_[iter->first];
                }
                int cut_off = std::abs(int(total_size - num_of_samples));
                if(cut_off != 0)
                {
                    StrataIndicesType::iterator end_iter = strata_indices_.begin();
                    for(int ii = 0; ii < cut_off; ++ii)
                    {
                        ++end_iter;
                    }
                    for(iter = strata_indices_.begin(); iter != end_iter; ++iter)
                    {
                        strata_sizes_[iter->first] -= 1;

                    }
                }

            }

            else
            {
                StrataIndicesType::iterator iter;
                unsigned int total_size = 0;
                for(iter = strata_indices_.begin(); iter != strata_indices_.end(); ++iter)
                {
                    // Set the size of the strata to be sampled to be proportional of the size of the strata
                    // This value is different for each strata
                    strata_sizes_[iter->first] = (int)ceil((double(iter->second.size())/double(max_index)) * num_of_samples );
                    total_size += strata_sizes_[iter->first];
                }
                int cut_off = std::abs(int(total_size - num_of_samples));
                if(cut_off != 0)
                {
                    for(int ii = 0; ii < cut_off; ++ii)
                    {
                        StrataIndicesType::iterator curmax = strata_indices_.begin();
                        for(iter = strata_indices_.begin(); iter != strata_indices_.end(); ++iter)
                        {
                            if(strata_sizes_[iter->first] > strata_sizes_[curmax->first])
                                curmax = iter;
                        }
                        strata_sizes_[curmax->first] -= 1;
                    }
                }
            }


              // Set sampling function
            if(options_.sample_with_replacement)
            {
                sampling_func_ptr_ = &Sampler::sample_stratified_w_rep;
            }
            else
            {
                // Check whether there are enough samples to sample the needed strata size from the population
                StrataSizesType::iterator iter;
                for(iter = strata_sizes_.begin(); iter != strata_sizes_.end(); ++iter)
                {
                    vigra_precondition(iter->second <= (int)strata_indices_[iter->first].size(),
                                                "Not enough samples to sample classes individually //stratified and\
                                                    without replacement");
                }
                sampling_func_ptr_ = &Sampler::sample_stratified_wo_rep;
            }
            // Allocate memory for output
            if(options_.use_internal_mem)
            {
                internal_used_indices_.resize(num_of_samples);
                used_indices_ = internal_used_indices_;
            }
        }
        else // unstratified sampling
        {

              // Set sampling function
            if(options_.sample_with_replacement)
            {
                if(options_.use_internal_mem)
                {
                    internal_used_indices_.resize(num_of_samples);
                    used_indices_ = internal_used_indices_;
                }
                sampling_func_ptr_ = &Sampler::sample_w_rep;
            }
            else
            {
                vigra_precondition(max_index >= num_of_samples,
                                        "Not enough samples to sample without replacement");
                if(options_.use_internal_mem)
                {
                    internal_used_indices_.resize(max_index);
                    used_indices_ = internal_used_indices_;
                    for(int ii = 0; ii < max_index; ++ii)
                    {
                        used_indices_[ii] = ii;
                    }
                }

                sampling_func_ptr_ = &Sampler::sample_wo_rep;

            }
        }
        if(options_.use_internal_mem)
        {
            used_indices_ = internal_used_indices_;
        }
    }

    inline Sampler(int numOfSamples,int maxIndex, SamplingOptions const & opt, Random & rnd)
    :
        randint(rnd)
    {
        init(numOfSamples,maxIndex, opt);
    }

    inline Sampler(int numOfSamples, int maxIndex, SamplingOptions const & opt)
    {
         init(numOfSamples,maxIndex, opt);
    }

};

template<class Random =RandomTT800 >
class PoissonSampler
{
    Random  randfloat;
    typedef Int32                               IndexType;
    typedef vigra::ArrayVector     <IndexType>  IndexArrayType;
    IndexArrayType        used_indices_;
    double lambda;
    int minIndex;
    int maxIndex;
    inline PoissonSampler(double lambda,IndexType minIndex,IndexType maxIndex,Random & rnd)
        : randfloat(rnd)
    {
        this->lambda=lambda;
        this->minIndex=minIndex;
        this->maxIndex=maxIndex;
    }

    inline void sample(  )
    {
        IndexType i;
        for(i=minIndex;i<maxIndex;++i)
        {
            //from http://en.wikipedia.org/wiki/Poisson_distribution
            int k=0;
            double p=1;
            double L=exp(-lambda);
            do
            {
                ++k;
                p*=randfloat.uniform53();

            }while(p>L);
            --k;
            //Insert i this many time
            while(k>0)
            {
                used_indices_.push_back(i);
                --k;
            }
        }
    }

    IndexType const & operator[](int in)
    {
        return used_indices_[in];
    }
    int numOfSamples()
    {
        return used_indices_.size();
    }

};

} /*end of namespace rn*/
#endif /*RN_SAMPLER_HXX*/
