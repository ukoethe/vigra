/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
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

#define VIGRA_CHECK_BOUNDS
#include "unittest.hxx"
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <vigra/mathutil.hxx>
#include <vigra/index_sampling.hxx>
#include <vigra/random_forest/rf_sampling.hxx>
#include <map>

using namespace vigra;

struct wierdOrderingFunctor
{
    bool operator()(vigra::ArrayVectorView<unsigned int> a, vigra::ArrayVectorView<unsigned int> b)
    {
        size_t ii = 0;
        while(ii < a.size() && a[ii] == b[ii])
        {
            ++ii;
        }
        if(ii == a.size())return false;
        else return a[ii] < b[ii];
    }
};

struct RFSamplerTests
{

    vigra::ArrayVector<int > strata;
    vigra::ArrayVector<int > strata_ungleich;

    RFSamplerTests()
    {
        for(int ii = 0; ii < 10; ++ii)
        {
            strata.push_back(1);
            strata.push_back(2);
        }
    }

    void testSamplingWithReplacement();
    void testSamplingWithoutReplacement();
    void testStratifiedSamplingWithReplacement();
    void testStratifiedSamplingWithoutReplacement();

    void testVolitality();

    void testSamplingWithReplacementChi2();
    void testSamplingWithoutReplacementChi2();

};



void RFSamplerTests::testSamplingWithoutReplacement()
{

    size_t SizeOfIndex = 20;

    size_t numOfSamples = int(SizeOfIndex/2);

    for(int ii = 0; ii <300; ++ii)
    {
        //Check whether samples are doubled or samplersize is != size given
        RFSampler<> sampler( numOfSamples,SizeOfIndex, RFSamplingOptions().sampleWithReplacement(false));
        sampler.sample();
        shouldEqualMessage(sampler.used_indices().size(), numOfSamples,
                            "number of samples in RFSampler != number of samples requested");
        {
            RFSampler<>::IndexArrayType usedIndices(sampler.used_indices());
            std::sort(usedIndices.begin(), usedIndices.end());
            for(unsigned int ii = 0; ii < usedIndices.size()-1; ++ii)
            {
                shouldMsg(sampler.used_indices()[ii] != sampler.used_indices()[ii+1],
                                "multiple entries while sampling without replacement");
            }

            should(*(usedIndices.end()-1) < int(SizeOfIndex));
            should(*(usedIndices.begin()) >= 0 );
        }


        // This test is probabilistic
        {
            sampler.sample();
            RFSampler<>::IndexArrayType lastSampledIndices(sampler.used_indices());
            sampler.sample();
            shouldMsg(lastSampledIndices != sampler.used_indices(),
                                    "Results of consecutive sampling is equal - Probability of this ist 1/10 ~ 2.8e-7!");
        }



        // Sample all all indices and check whether all indices are present
        {
            RFSampler<>::IndexArrayType orderedIndex;
            for(unsigned int ii = 0; ii < SizeOfIndex; ++ii)
            {
                orderedIndex.push_back(ii);
            }
            sampler.init(SizeOfIndex,SizeOfIndex, RFSamplingOptions().sampleWithoutReplacement());
            sampler.sample();
            RFSampler<>::IndexArrayType usedIndices(sampler.used_indices());

            std::sort(usedIndices.begin(), usedIndices.end());
            shouldEqual(usedIndices, orderedIndex);
        }
    }
}

void RFSamplerTests::testVolitality()
{
    RFSampler<> sampler(10,10, RFSamplingOptions());
    sampler.sample();
    //Check for constness of the used_indices method and unconstness of used_indices_volatile
    try
    {
        sampler.used_indices()[1] = 20;
        shouldMsg(1 != 1, "Indices are not const");
    }
    catch(...)
    {
        sampler.used_indices_volatile()[1] = 30;
    }

    shouldEqual(sampler.used_indices()[1], (30));
}


void RFSamplerTests::testSamplingWithReplacement()
{

    size_t SizeOfIndex = 20;
    size_t numOfSamples = int(SizeOfIndex/2);

    for(int ii = 0; ii <300; ++ii)
    {
        //Check whether to big sample indices are created or samplersize is != size given
        RFSampler<> sampler(numOfSamples,SizeOfIndex, RFSamplingOptions().sampleWithReplacement());
        sampler.sample();
        shouldEqualMessage(sampler.used_indices().size(), numOfSamples,
                            "number of samples in RFSampler != number of samples requested");
        {
            RFSampler<>::IndexArrayType usedIndices(sampler.used_indices());
            std::sort(usedIndices.begin(), usedIndices.end());
            should(*(usedIndices.end()-1) < int(SizeOfIndex));
            should(*(usedIndices.begin()) >= 0 );
        }


        // This test is probabilistic
        {
            sampler.sample();
            RFSampler<>::IndexArrayType lastSampledIndices(sampler.used_indices());
            sampler.sample();
            shouldMsg(lastSampledIndices != sampler.used_indices(),
                                    "Results of consecutive sampling is equal - Probability of this ist 10^10 ~ 1e-10!");
        }
    }
}


void RFSamplerTests::testStratifiedSamplingWithReplacement()
{
    //gleiche Daten;
    vigra::ArrayVector<int> strata;
    for(int ii = 0; ii < 10; ++ii)
    {
        strata.push_back(1);
        strata.push_back(2);
    }

    RFSampler<> sampler( 10, 20, RFSamplingOptions().sampleWithReplacement().sampleStratified(strata));
    sampler.sample();
    RFSampler<> sampler2( 10, 20, RFSamplingOptions().sampleWithReplacement().sampleClassesIndividually(strata));
    sampler2.sample();

    for(int ii = 0; ii < 5; ++ii)
    {
        shouldEqual(strata[sampler.used_indices()[ii]], 1);
        shouldEqual(strata[sampler.used_indices()[5+ii]], 2);
        shouldEqual(strata[sampler2.used_indices()[ii]], 1);
        shouldEqual(strata[sampler2.used_indices()[5+ii]], 2);
    }

    vigra::ArrayVector<int> strata1_2;
    for(int ii = 0; ii < 10; ++ii)
    {
        strata1_2.push_back(1);
        strata1_2.push_back(2);
        strata1_2.push_back(2);
    }

    //ungleiche daten;
    sampler.init( 9,30, RFSamplingOptions().sampleWithReplacement().sampleStratified(strata1_2));
    sampler.sample();
    sampler2.init( 9,30, RFSamplingOptions().sampleWithReplacement().sampleClassesIndividually(strata1_2));
    sampler2.sample();
    for(int ii = 0; ii < 3; ++ii)
    {
        shouldEqual(strata1_2[sampler.used_indices()[ii]], 1);
        shouldEqual(strata1_2[sampler.used_indices()[3+ii]], 2);
        shouldEqual(strata1_2[sampler.used_indices()[6+ii]], 2);
    }
    for(int ii = 0; ii < 4; ++ii)
    {
        shouldEqual(strata1_2[sampler2.used_indices()[ii]], 1);
        shouldEqual(strata1_2[sampler2.used_indices()[5+ii]], 2);
    }
    shouldEqual(strata1_2[sampler2.used_indices()[4]], 2);

}

void RFSamplerTests::testStratifiedSamplingWithoutReplacement()
{
    //gleiche Daten;
    vigra::ArrayVector<int> strata;
    for(int ii = 0; ii < 10; ++ii)
    {
        strata.push_back(1);
        strata.push_back(2);
    }

    RFSampler<> sampler( 10, 20, RFSamplingOptions().sampleWithoutReplacement().sampleStratified(strata));
    sampler.sample();
    RFSampler<> sampler2( 10, 20, RFSamplingOptions().sampleWithoutReplacement().sampleClassesIndividually(strata));
    sampler2.sample();

    for(int ii = 0; ii < 5; ++ii)
    {
        shouldEqual(strata[sampler.used_indices()[ii]], 1);
        shouldEqual(strata[sampler.used_indices()[5+ii]], 2);
        shouldEqual(strata[sampler2.used_indices()[ii]], 1);
        shouldEqual(strata[sampler2.used_indices()[5+ii]], 2);
    }

    vigra::ArrayVector<int> strata1_2;
    for(int ii = 0; ii < 10; ++ii)
    {
        strata1_2.push_back(1);
        strata1_2.push_back(2);
        strata1_2.push_back(2);
    }

    //ungleiche daten;

    sampler.init(9, 30, RFSamplingOptions().sampleWithoutReplacement().sampleStratified(strata1_2));

    sampler.sample();
    sampler2.init( 9, 30, RFSamplingOptions().sampleWithoutReplacement().sampleClassesIndividually(strata1_2));
    sampler2.sample();

    for(int ii = 0; ii < 3; ++ii)
    {

        shouldEqual(strata1_2[sampler.used_indices()[ii]], 1);
        shouldEqual(strata1_2[sampler.used_indices()[3+ii]], 2);
        shouldEqual(strata1_2[sampler.used_indices()[6+ii]], 2);
    }
    for(int ii = 0; ii < 4; ++ii)
    {
        shouldEqual(strata1_2[sampler2.used_indices()[ii]], 1);
        shouldEqual(strata1_2[sampler2.used_indices()[5+ii]], 2);
    }
    shouldEqual(strata1_2[sampler2.used_indices()[4]], 2);


    //nicht genug daten in einer Klasse
    try
    {
        sampler2.init( 24, 30, RFSamplingOptions().sampleWithReplacement().sampleClassesIndividually(strata1_2));
        sampler2.sample();
        shouldMsg(1!= 1, "No Exceptionthrown while Sampling Classes Individually if not enough Samples are available");
    }
    catch(...)
    {

    }

}

void RFSamplerTests::testSamplingWithReplacementChi2()
{
    vigra::ArrayVector<int> observed(10);
    observed.init(0);
    size_t numOfSamples = 1000;
    size_t maxIndex = 10;
    double chi_squared = 0;

    {
    RFSampler<> sampler( numOfSamples, maxIndex, RFSamplingOptions().sampleWithReplacement(true));
    sampler.sample();
        RFSampler<>::IndexArrayType::iterator iter;

        for(iter = sampler.used_indices().begin(); iter!= sampler.used_indices().end(); ++iter)
        {
            observed[*iter]++;
        }
        for(int ii = 0; ii < 10; ++ii)
        {
            chi_squared += double(((observed[ii] - 100.0)*(observed[ii] - 100.0))/100.0);
        }
        should(chi_squared < 21.67);
    }
    maxIndex = 10000;
    {
    RFSampler<> sampler( maxIndex, maxIndex, RFSamplingOptions().sampleWithReplacement(true));
    sampler.sample();
        double numPositives = 0;
        RFSampler<>::IsUsedArrayType::iterator iter;
        for(iter = sampler.is_used().begin(); iter != sampler.is_used().end(); ++iter)
        {
            if(*iter)numPositives++;
        }
        numPositives = numPositives/double(maxIndex);

        /* when sampling k times without replacement
        the probability p of a sample not being pickes is ((k-1)/k)^k
        The Distribution of the number of samples not being choosen is a
        Binomial distribution with n = k and p = ((k-1)/k)^k.
        The expectation value of a Binomial Distribution is n*p ==>
        The percentage of samples not used is (n*p)/k = p
        For large n p converges to 0.366 => numpositives should be
        around 0.63 +/- 0.01 */
        shouldEqualTolerance (numPositives, 0.63, 0.01);
    }
}


void RFSamplerTests::testSamplingWithoutReplacementChi2()
{
    double nsamples = 10000;
    double nclasses = 120;
    RFSampler<> sampler( 5, 5, RFSamplingOptions().sampleWithoutReplacement());
    std::map<unsigned int, int> wierdmap;
    std::map<unsigned int , int>::iterator iter;
    for(int ii = 0; ii < 1000; ++ii)
    {
        sampler.sample();
        int dec = 1;
        unsigned int hash = 0;
        for(size_t ii = 0; ii < sampler.used_indices().size(); ++ii)
        {
            hash += dec* sampler.used_indices()[ii];
            dec = dec *10;
        }
        wierdmap[hash] = 0;
    }
    shouldEqual(wierdmap.size(), nclasses);
    for(int ii = 0; ii < nsamples; ++ii)
    {
        sampler.sample();
        int dec = 1;
        unsigned int hash = 0;
        for(size_t ii = 0; ii < sampler.used_indices().size(); ++ii)
        {
            hash += dec* sampler.used_indices()[ii];
            dec = dec *10;
        }
        wierdmap[hash] += 1;
    }
    double chi_squared = 0;
    for(iter = wierdmap.begin(); iter != wierdmap.end(); ++iter)
    {
        //std::cerr << iter->second<< "  ";
            chi_squared += double(((iter->second - (nsamples/nclasses))*(iter->second - (nsamples/nclasses)))/(nsamples/nclasses));
    }
    //std::cerr << "chi squared while sampling wo replacement :" << std::endl << chi_squared << std::endl;
    should(chi_squared < 117.7995);

}

/***********************************************************************/

struct SamplerTests
{
    SamplerTests()
    {}

    void testSamplingWithoutReplacement();
    void testStratifiedSamplingWithoutReplacement();
    void testSamplingWithReplacement();
    void testStratifiedSamplingWithReplacement();
    void testSamplingWithoutReplacementChi2();
    void testSamplingWithReplacementChi2();
    
    void testSamplingImpl(bool withReplacement);
    void testStratifiedSamplingImpl(bool withReplacement);
};

void SamplerTests::testSamplingWithoutReplacement()
{
    testSamplingImpl(false);
}

void SamplerTests::testSamplingWithReplacement()
{
    testSamplingImpl(true);
}

void SamplerTests::testSamplingImpl(bool withReplacement)
{
    int totalDataCount = 20;
    int numOfSamples = int(totalDataCount/2);

    for(int ii = 0; ii <300; ++ii)
    {
        // check basic attributes
        Sampler<> sampler( totalDataCount, 
            SamplingOptions().withReplacement(withReplacement).sampleCount(numOfSamples));
        shouldEqual(sampler.totalCount(), totalDataCount);
        shouldEqual(sampler.sampleCount(), numOfSamples);
        shouldEqual(sampler.strataCount(), 1);
        shouldEqual(sampler.stratifiedSampling(), false);
        shouldEqual(sampler.withReplacement(), withReplacement);

        Sampler<> samplerProportional( totalDataCount, 
            SamplingOptions().withReplacement(withReplacement).sampleProportion(0.5));
        shouldEqual(samplerProportional.totalCount(), totalDataCount);
        shouldEqual(samplerProportional.sampleCount(), numOfSamples);
        shouldEqual(samplerProportional.strataCount(), 1);
        shouldEqual(samplerProportional.withReplacement(), withReplacement);
        
        // check that indices are either sampled or out-of-bag
        {
            ArrayVector<bool> wasPicked(totalDataCount, false);
            Sampler<>::IndexArrayType usedIndices(sampler.sampledIndices());
            Sampler<>::IndexArrayType unusedIndices(sampler.oobIndices());
            
            shouldEqual(usedIndices.size(), numOfSamples);
            if(withReplacement)
                should(usedIndices.size()+unusedIndices.size() >= (unsigned int)totalDataCount);
            else
                shouldEqual(usedIndices.size()+unusedIndices.size(), totalDataCount);
                
            for(unsigned int ii = 0; ii < usedIndices.size(); ++ii)
            {
                should(usedIndices[ii] >= 0 && usedIndices[ii] < int(totalDataCount));
                wasPicked[usedIndices[ii]] = true;
            }
            for(unsigned int ii = 0; ii < unusedIndices.size(); ++ii)
            {
                should(unusedIndices[ii] >= 0 && unusedIndices[ii] < int(totalDataCount));
                wasPicked[unusedIndices[ii]] = true;
            }
            for(int ii = 0; ii < totalDataCount; ++ii)
            {
                should(wasPicked[ii]);
            }
        }

        // check that consecutive samples differ
        {
            Sampler<>::IndexArrayType lastSampledIndices(sampler.sampledIndices());
            sampler.sample();
            shouldMsg(lastSampledIndices != sampler.sampledIndices(),
                       "Consecutive samples are equal - probability of this is less than 10^-5.");
        }
    }

    if(!withReplacement)
    {
        // check exception when more samples than data
        try
        {
            Sampler<> sampler( 2, 
                 SamplingOptions().withoutReplacement().sampleCount(9));
            failTest("No exception thrown when there are too few data for sampling.");
        }
        catch(PreconditionViolation &)
        {}
    }
}



void SamplerTests::testStratifiedSamplingWithoutReplacement()
{
    testStratifiedSamplingImpl(false);
}

void SamplerTests::testStratifiedSamplingWithReplacement()
{
    testStratifiedSamplingImpl(true);
}

void SamplerTests::testStratifiedSamplingImpl(bool withReplacement)
{
    vigra::ArrayVector<int> strata;
    for(int ii = 0; ii < 10; ++ii)
    {
        strata.push_back(1);
        strata.push_back(2);
    }
    
    {
        int  totalDataCount = strata.size();
        Sampler<> sampler( strata.begin(), strata.end(), 
             SamplingOptions().withReplacement(withReplacement).sampleCount(10).stratified());
        shouldEqual(sampler.totalCount(), totalDataCount);
        shouldEqual(sampler.sampleCount(), 10);
        shouldEqual(sampler.strataCount(), 2);
        shouldEqual(sampler.stratifiedSampling(), true);
        shouldEqual(sampler.withReplacement(), withReplacement);
        if(withReplacement)
            should(sampler.sampledIndices().size()+sampler.oobIndices().size() >= (unsigned int)totalDataCount);
        else
            shouldEqual(sampler.sampledIndices().size()+sampler.oobIndices().size(), totalDataCount);
            
        ArrayVector<bool> wasPicked(totalDataCount, false);

        for(int ii = 0; ii < 5; ++ii)
        {
            int index = sampler.sampledIndices()[ii];
            should(index >= 0 && index < int(totalDataCount));
            shouldEqual(strata[index], 1);
            wasPicked[index] = true;
        }
        for(int ii = 5; ii < 10; ++ii)
        {
            int index = sampler.sampledIndices()[ii];
            should(index >= 0 && index < int(totalDataCount));
            shouldEqual(strata[index], 2);
            wasPicked[index] = true;
        }
        for(unsigned int ii = 0; ii < sampler.oobIndices().size(); ++ii)
        {
            int index = sampler.oobIndices()[ii];
            should(index >= 0 && index < int(totalDataCount));
            wasPicked[index] = true;
        }
        for(int ii = 0; ii < totalDataCount; ++ii)
        {
            should(wasPicked[ii]);
        }
        Sampler<>::IndexArrayType lastSampledIndices(sampler.sampledIndices());
        sampler.sample();
        should(lastSampledIndices != sampler.sampledIndices());
    }

    {
        int  totalDataCount = strata.size();
        Sampler<> sampler( strata.begin(), strata.end(), 
             SamplingOptions().withReplacement(withReplacement).sampleCount(9).stratified());
        shouldEqual(sampler.sampleCount(), 9);

        ArrayVector<bool> wasPicked(totalDataCount, false);
        for(int ii = 0; ii < 4; ++ii)
        {
            int index = sampler.sampledIndices()[ii];
            should(index >= 0 && index < int(totalDataCount));
            shouldEqual(strata[index], 1);
            wasPicked[index] = true;
        }
        for(int ii = 4; ii < 9; ++ii)
        {
            int index = sampler.sampledIndices()[ii];
            should(index >= 0 && index < int(totalDataCount));
            shouldEqual(strata[index], 2);
            wasPicked[index] = true;
        }
        for(unsigned int ii = 0; ii < sampler.oobIndices().size(); ++ii)
        {
            int index = sampler.oobIndices()[ii];
            should(index >= 0 && index < int(totalDataCount));
            wasPicked[index] = true;
        }
        for(int ii = 0; ii < totalDataCount; ++ii)
        {
            should(wasPicked[ii]);
        }
    }

    for(int ii = 0; ii < 10; ++ii)
    {
        strata.push_back(1);
    }

    {
        int  totalDataCount = strata.size();
        Sampler<> sampler( strata.begin(), strata.end(), 
             SamplingOptions().withReplacement(withReplacement).sampleCount(10).stratified());

        for(int ii = 0; ii < 5; ++ii)
        {
            shouldEqual(strata[sampler.sampledIndices()[ii]], 1);
        }
        for(int ii = 5; ii < 10; ++ii)
        {
            shouldEqual(strata[sampler.sampledIndices()[ii]], 2);
        }
    }

    for(int ii = 0; ii < 10; ++ii)
    {
        strata.push_back(3);
    }

    // need at most one sample per stratum
    try
    {
        Sampler<> sampler( strata.begin(), strata.end(), 
             SamplingOptions().withReplacement(withReplacement).sampleCount(2).stratified());
        failTest("No exception thrown when there are too few data for stratified sampling.");
    }
    catch(PreconditionViolation &)
    {}
}

void SamplerTests::testSamplingWithoutReplacementChi2()
{
    // check that all permutations of the indices occur with about equal frequency
    // use fixed random numbers so that the sampling is reproducible 
    int nsamples = 120000;
    int nclasses = 120;
    Sampler<> sampler( 5, 
                SamplingOptions().withoutReplacement().sampleCount(5),
                MersenneTwister(987928));
    std::map<unsigned int, int> wierdmap;
    std::map<unsigned int , int>::iterator iter;
    for(int ii = 0; ii < 1000; ++ii)
    {
        sampler.sample();
        int dec = 1;
        unsigned int hash = 0;
        for(size_t ii = 0; ii < sampler.sampledIndices().size(); ++ii)
        {
            hash += dec* sampler.sampledIndices()[ii];
            dec = dec *10;
        }
        wierdmap[hash] = 0;
    }
    shouldEqual(wierdmap.size(), nclasses);
    for(int ii = 0; ii < nsamples; ++ii)
    {
        sampler.sample();
        int dec = 1;
        unsigned int hash = 0;
        for(size_t ii = 0; ii < sampler.sampledIndices().size(); ++ii)
        {
            hash += dec* sampler.sampledIndices()[ii];
            dec = dec *10;
        }
        wierdmap[hash] += 1;
    }
    double chi_squared = 0;
    double ratio = nsamples/nclasses;
    for(iter = wierdmap.begin(); iter != wierdmap.end(); ++iter)
    {
        chi_squared += sq(iter->second - ratio)/ratio;
    }
    should(chi2CDF(119, chi_squared) < 0.1);
}

void SamplerTests::testSamplingWithReplacementChi2()
{
    vigra::ArrayVector<int> observed(10);
    observed.init(0);
    int totalDataCount = 10;
    int numOfSamples = 1000000;
    double chi_squared = 0,
           ratio = double(numOfSamples) / totalDataCount;

    {
        Sampler<> sampler(totalDataCount, 
             SamplingOptions().withReplacement().sampleCount(numOfSamples),
             MersenneTwister());

        for(int ii = 0; ii < numOfSamples; ++ii)
        {
            observed[sampler.sampledIndices()[ii]]++;
        }
        for(int ii = 0; ii < totalDataCount; ++ii)
        {
            chi_squared += sq(observed[ii] - ratio)/ratio;
        }
        should(chi2CDF(9, chi_squared) < 0.1);
    }
    totalDataCount = 10000;
    {
        Sampler<> sampler( totalDataCount, 
             SamplingOptions().withReplacement().sampleCount(totalDataCount));
        double numPositives = double(totalDataCount - sampler.oobIndices().size()) / totalDataCount;

        /* when sampling k times without replacement
        the probability p of a sample not being pickes is ((k-1)/k)^k
        The Distribution of the number of samples not being choosen is a
        Binomial distribution with n = k and p = ((k-1)/k)^k.
        The expectation value of a Binomial Distribution is n*p ==>
        The percentage of samples not used is (n*p)/k = p
        For large n p converges to 0.366 => numpositives should be
        around 0.63 +/- 0.01 */
        shouldEqualTolerance (0, numPositives-0.63, 0.01);
    }
}

struct SamplerTestSuite
: public vigra::test_suite
{
    SamplerTestSuite()
    : vigra::test_suite("Sampler Test")
    {
        add(testCase(&RFSamplerTests::testSamplingWithReplacement));
        add(testCase(&RFSamplerTests::testSamplingWithReplacementChi2));
        add(testCase(&RFSamplerTests::testSamplingWithoutReplacement));
        add(testCase(&RFSamplerTests::testSamplingWithoutReplacementChi2));
        add(testCase(&RFSamplerTests::testStratifiedSamplingWithReplacement));
        add(testCase(&RFSamplerTests::testStratifiedSamplingWithoutReplacement));
        add(testCase(&RFSamplerTests::testVolitality));

        add(testCase(&SamplerTests::testSamplingWithoutReplacement));
        add(testCase(&SamplerTests::testStratifiedSamplingWithoutReplacement));
        add(testCase(&SamplerTests::testSamplingWithReplacement));
        add(testCase(&SamplerTests::testStratifiedSamplingWithReplacement));
        add(testCase(&SamplerTests::testSamplingWithoutReplacementChi2));
        add(testCase(&SamplerTests::testSamplingWithReplacementChi2));
    }
};

int main(int argc, char **argv)
{

    SamplerTestSuite samplerTest;

    int failed = samplerTest.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << samplerTest.report() << std::endl;


    return (failed != 0);
}
