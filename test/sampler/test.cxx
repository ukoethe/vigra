/************************************************************************/
/*                                                                      */
/*                 Copyright 2013 by Ullrich Koethe                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
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

#define VIGRA_CHECK_BOUNDS
#include "vigra/unittest.hxx"
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <vigra/mathutil.hxx>
#include <vigra/sampling.hxx>
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
            SamplerOptions().withReplacement(withReplacement).sampleSize(numOfSamples));
        shouldEqual(sampler.totalCount(), totalDataCount);
        shouldEqual(sampler.sampleSize(), numOfSamples);
        shouldEqual(sampler.strataCount(), 1);
        shouldEqual(sampler.stratifiedSampling(), false);
        shouldEqual(sampler.withReplacement(), withReplacement);

        Sampler<> samplerProportional( totalDataCount, 
            SamplerOptions().withReplacement(withReplacement).sampleProportion(0.5));
        shouldEqual(samplerProportional.totalCount(), totalDataCount);
        shouldEqual(samplerProportional.sampleSize(), numOfSamples);
        shouldEqual(samplerProportional.strataCount(), 1);
        shouldEqual(samplerProportional.withReplacement(), withReplacement);
        
        sampler.sample();
        samplerProportional.sample();
        // check that indices are either sampled or out-of-bag
        {
            ArrayVector<bool> wasPicked(totalDataCount, false);
            Sampler<>::IndexArrayType usedIndices(sampler.sampledIndices());
            Sampler<>::IndexArrayType unusedIndices(sampler.oobIndices());
            
            shouldEqual((int)usedIndices.size(), numOfSamples);
            if(withReplacement)
                should(usedIndices.size()+unusedIndices.size() >= (unsigned int)totalDataCount);
            else
                shouldEqual(usedIndices.size()+unusedIndices.size(), (unsigned int)totalDataCount);
                
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
                 SamplerOptions().withoutReplacement().sampleSize(9));
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
             SamplerOptions().withReplacement(withReplacement).sampleSize(10).stratified());
        shouldEqual(sampler.totalCount(), totalDataCount);
        shouldEqual(sampler.sampleSize(), 10);
        shouldEqual(sampler.strataCount(), 2);
        shouldEqual(sampler.stratifiedSampling(), true);
        shouldEqual(sampler.withReplacement(), withReplacement);
        sampler.sample();
        if(withReplacement)
            should(int(sampler.sampledIndices().size()+sampler.oobIndices().size()) >= totalDataCount);
        else
            shouldEqual(int(sampler.sampledIndices().size()+sampler.oobIndices().size()), totalDataCount);
            
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
            should(index >= 0 && index < totalDataCount);
            shouldEqual(strata[index], 2);
            wasPicked[index] = true;
        }
        for(int ii = 0; ii < (int)sampler.oobIndices().size(); ++ii)
        {
            int index = sampler.oobIndices()[ii];
            should(index >= 0 && index < totalDataCount);
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
        int totalDataCount = strata.size();
        Sampler<> sampler( strata.begin(), strata.end(), 
             SamplerOptions().withReplacement(withReplacement).sampleSize(9).stratified());
        sampler.sample();
        shouldEqual(sampler.sampleSize(), 9);

        ArrayVector<bool> wasPicked(totalDataCount, false);
        for(int ii = 0; ii < 4; ++ii)
        {
            int index = sampler.sampledIndices()[ii];
            should(index >= 0 && index < totalDataCount);
            shouldEqual(strata[index], 1);
            wasPicked[index] = true;
        }
        for(int ii = 4; ii < 9; ++ii)
        {
            int index = sampler.sampledIndices()[ii];
            should(index >= 0 && index < totalDataCount);
            shouldEqual(strata[index], 2);
            wasPicked[index] = true;
        }
        for(int ii = 0; ii < (int)sampler.oobIndices().size(); ++ii)
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
        Sampler<> sampler( strata.begin(), strata.end(), 
             SamplerOptions().withReplacement(withReplacement).sampleSize(10).stratified());
        sampler.sample();

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
             SamplerOptions().withReplacement(withReplacement).sampleSize(2).stratified());
        failTest("No exception thrown when there are too few data for stratified sampling.");
    }
    catch(PreconditionViolation &)
    {}
}

void SamplerTests::testSamplingWithoutReplacementChi2()
{
    // Check that all permutations of the indices occur with about equal frequency.
    // (Residuals are Poisson distributed, which is approximated by a Gaussian
    //  distribution with data-dependent variance, so that conformance can be checked
    //  by a chi-square test.)
    // Use fixed random numbers so that the sampling is reproducible.
    int nsamples = 120000;
    int nclasses = 120;
    MersenneTwister randomGenerator;
    Sampler<> sampler( 5, 
                SamplerOptions().withoutReplacement().sampleSize(5),
                &randomGenerator);
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
    
    // check that all 120 permutations occured after 1000 trials
    shouldEqual((int)wierdmap.size(), nclasses);
    
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
    
    // check that we are in the 80% quantile of the expected distribution
    shouldEqualTolerance (0, chi2CDF(119, chi_squared)-0.5, 0.4);
}

void SamplerTests::testSamplingWithReplacementChi2()
{
    // Check that samples are selected with uniform probability
    // (Residuals are Poisson distributed, which is approximated by a Gaussian
    //  distribution with data-dependent variance, so that conformance can be checked
    //  by a chi-square test.)
    // Use fixed random numbers so that the sampling is reproducible.
    vigra::ArrayVector<int> observed(10);
    observed.init(0);
    int totalDataCount = 10;
    int numOfSamples = 100000;
    double chi_squared = 0,
           ratio = double(numOfSamples) / totalDataCount;

    {
        MersenneTwister randomGenerator;
        Sampler<> sampler(totalDataCount, 
             SamplerOptions().withReplacement().sampleSize(numOfSamples),
             &randomGenerator);

        sampler.sample();
        for(int ii = 0; ii < numOfSamples; ++ii)
        {
            observed[sampler.sampledIndices()[ii]]++;
        }
        for(int ii = 0; ii < totalDataCount; ++ii)
        {
            chi_squared += sq(observed[ii] - ratio)/ratio;
        }
        // check that we are in the 80% quantile of the expected distribution
        shouldEqualTolerance (0, chi2CDF(9, chi_squared)-0.5, 0.4);
    }

    /* when sampling k times without replacement
    the probability p of a sample not being picked is ((k-1)/k)^k
    The distribution of the number of samples not being chosen is a
    binomial distribution with n = k and p = ((k-1)/k)^k.
    The expectation value of a binomial distribution is n*p ==>
    The percentage of samples not used is (n*p)/k = p
    For large n p converges to 0.366 => numpositives should be
    around 0.63 +/- 0.015 */
    totalDataCount = 10000;
    {
        Sampler<> sampler( totalDataCount, 
             SamplerOptions().withReplacement().sampleSize(totalDataCount));
        sampler.sample();
        double numPositives = double(totalDataCount - sampler.oobIndices().size()) / totalDataCount;

        shouldEqualTolerance (0, numPositives-0.63, 0.015);
    }
}

struct SamplerTestSuite
: public vigra::test_suite
{
    SamplerTestSuite()
    : vigra::test_suite("Sampler Test")
    {
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
