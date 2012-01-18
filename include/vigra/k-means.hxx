/************************************************************************/
/*                                                                      */
/*               Copyright 2011-2012 by Ullrich Koethe                  */
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

#ifndef VIGRA_KMEANS_HXX
#define VIGRA_KMEANS_HXX

#include "multi_array.hxx"
#include "multi_math.hxx"
#include "matrix.hxx"
#include "sampling.hxx"
#include <vector>
#include <algorithm>

namespace vigra {

    /** \brief Randomly initialize cluster centers for k-means.
    
        Array \a data contains the data to be clustered, and the initial 
        cluster centers will be written into \a clusterCenters. The number of
        features corresponds to the length of the first dimension of \a data and
        \a clusterCenters, which must match. The desired number of clusters
        corresponds to the length of the second dimension in \a clusterCenters
        and must be at least 2. This function chooses centers 
        by sampling random points (columns) from the data array without replacement.
        
        See \ref kMeans() for a usage example.
    */
template<class T1, class Stride1, class T2, class Stride2>
void 
kMeansRandomInitialization(MultiArrayView<2, T1, Stride1> const & data, 
                           MultiArrayView<2, T2, Stride2> & clusterCenters)
{
    int numFeatures = rowCount(data);
    int numSamples = columnCount(data);
    int clusterCount = columnCount(clusterCenters);
    vigra_precondition(clusterCount >= 2,
      "kMeansRandomInitialization(): At least 2 clusters are required.");
    vigra_precondition(numFeatures == rowCount(clusterCenters),
      "kMeansRandomInitialization(): Number of features mismatch between data and clusterCenter arrays.");

    Sampler<> sampler(numSamples, 
                      SamplerOptions().sampleSize(clusterCount).withoutReplacement());
    sampler.sample();
    for(int k=0; k<sampler.sampleSize(); ++k)
    {
        columnVector(clusterCenters,k) = columnVector(data,sampler[k]);
    }
}

    /** \brief Initialize cluster centers according to the k-means++ algorithm.
    
        Array \a data contains the data to be clustered, and the initial 
        cluster centers will be written into \a clusterCenters. The number of
        features corresponds to the length of the first dimension of \a data and
        \a clusterCenters, which must match. The desired number of clusters
        corresponds to the length of the second dimension in \a clusterCenters
        and must be at least 2. This function chooses centers 
        by sampling points (columns) from \a data with a probability that increases 
        when the point is far away from the centers selected so far,
        as descibed in:
        
        D. Arthur, D. and S. Vassilvitskii: <em>k-means++: the advantages of careful seeding</em>, SIAM Symposium on Discrete Algorithms, 2007, pp. 1027-1035.
        
        See \ref kMeans() for a usage example.
    */
template<class T1, class Stride1, class T2, class Stride2>
void 
kMeansPlusPlusInitialization(MultiArrayView<2, T1, Stride1> const & data, 
                             MultiArrayView<2, T2, Stride2> & clusterCenters)
{
    using namespace multi_math;

    MultiArrayIndex numFeatures = rowCount(data);
    MultiArrayIndex numSamples = columnCount(data);
    MultiArrayIndex clusterCount = columnCount(clusterCenters);
    vigra_precondition(clusterCount >= 2,
      "kMeansPlusPlusInitialization(): At least 2 clusters are required.");
    vigra_precondition(numFeatures == rowCount(clusterCenters),
      "kMeansPlusPlusInitialization(): Number of features mismatch between data and clusterCenter arrays.");
    
    typedef typename NormTraits<MultiArrayView<2, T2, Stride2> >::SquaredNormType DistType;
    std::vector<DistType> distances(numSamples, NumericTraits<DistType>::max());
    
    MersenneTwister random(RandomSeed);
    
    columnVector(clusterCenters, 0) = columnVector(data, random.uniformInt(data.size()));
    for(MultiArrayIndex k=1; k<clusterCount; ++k)
    {
        DistType distSum = NumericTraits<DistType>::zero();
        MultiArrayIndex s=0;
        for(s=0; s < numSamples; ++s)
        {
            DistType dist = sum(sq(columnVector(clusterCenters,k-1) - columnVector(data, s)), DistType());
            if(dist < distances[s])
                distances[s] = dist;
            distSum += distances[s];
        }
        
        DistType threshold = distSum*random.uniform();
        distSum = NumericTraits<DistType>::zero();
        for(s=0; s < numSamples; ++s)
        {
            distSum += distances[s];
            if(distSum >= threshold)
                break;
        }
        if(s < numSamples)
            columnVector(clusterCenters,k) = columnVector(data,s);
    }
}

    /** \brief Perform k-means clustering of an array.
    
        Array \a data holds the data to be clustered, and the array shape must be 
        <tt>(numFeatures * numSamples)</tt>. Cluster assignments
        (i.e. the index of the nearest cluster to each point) are returned
        in array \a assignments, whose length must be equal to <tt>numSamples</tt>. 
        The cluster centers (i.e. the mean of the points assigned to each cluster) 
        is returned in array \a clusterCenters whose shape must be
        <tt>(numFeatures * numClusters)</tt>.
        
        Array \a clusterCenters also serves as input: the length of its second dimension
        determines the number of clusters to be created, and its contents are the initial
        cluster centers needed by k-means.
        
        <b> Usage:</b>

        <b>\#include</b> \<vigra/k-means.hxx\><br>
        Namespace: vigra
        
        \code
        int numFeatures = 3, numSamples = 1000;
        MultiArray<2, double> data(Shape2(numFeatures, numSamples)); 
        ... // fill the data array
        
        MultiArray<1, unsigned int> assignments(numSamples);
        
        int clusterCount = 200;
        MultiArray<2, double> clusterCenters(Shape2(numFeatures, clusterCount)); 
        
        // initialize the cluster centers
        // alternatively, you may use kMeansPlusPlusInitialization() or
        // your own initialization code
        kMeansRandomInitialization(data, clusterCenters);
        
        // perform k-means clustering
        kMeans(data, assignments, clusterCenters);
        
        // replace the initial data values with their cluster representatives
        for(int k=0; k<numSamples; ++k)
            columnVector(data, k) = columnVector(clusterCenters, assignments(k));
        \endcode
        
        See also: \ref kMeansRandomInitialization(), \ref kMeansPlusPlusInitialization()
    */
template<class T1, class Stride1, class T2, class Stride2, class T3, class Stride3>
void 
kMeans(MultiArrayView<2, T1, Stride1> const & data, 
       MultiArrayView<1, T2, Stride2> & assignments,
       MultiArrayView<2, T3, Stride3> & clusterCenters)
{
    using namespace multi_math;

    MultiArrayIndex numFeatures = rowCount(data);
    MultiArrayIndex numSamples = columnCount(data);
    MultiArrayIndex clusterCount = columnCount(clusterCenters);
    vigra_precondition(clusterCount >= 2,
      "kMeans(): At least 2 clusters are required.");
    vigra_precondition(numFeatures == rowCount(clusterCenters),
      "kMeans(): Number of features mismatch between data and clusterCenter arrays.");
    vigra_precondition(numSamples == assignments.size(),
      "kMeans(): Number of samples mismatch between data and assignments arrays.");
    
    typedef typename NormTraits<MultiArrayView<2, T3, Stride3> >::SquaredNormType DistType;
    
    std::vector<MultiArrayIndex> clusterSizes(clusterCount);
    assignments.init(clusterCount+1);
    
    int maxIter = 100, iter=0;
    for(; iter<maxIter; ++iter)
    {
        bool assignmentsChanged = false;
        
        // assign all points to the nearest cluster
        for(MultiArrayIndex s=0; s < numSamples; ++s)
        {
            MultiArrayIndex nearestCluster = 0;
            DistType minDist = sum(sq(columnVector(clusterCenters,0) - columnVector(data, s)), DistType());
            
            for(MultiArrayIndex k=1; k<clusterCount; ++k)
            {
                DistType dist = sum(sq(columnVector(clusterCenters,k) - columnVector(data, s)), DistType());
                if(dist < minDist)
                {
                    minDist = dist;
                    nearestCluster = k;
                }
            }
            
            if(assignments(s) != nearestCluster)
            {
                assignments(s) = nearestCluster;
                assignmentsChanged = true;
            }
        }
        
        // stop if assignments didn't change in the last iteration
        if(!assignmentsChanged)
            break;
        
        // update cluster centers
        clusterCenters.init(T3());
        std::fill(clusterSizes.begin(), clusterSizes.end(), 0);
        
        for(MultiArrayIndex s=0; s < numSamples; ++s)
        {
            columnVector(clusterCenters, (MultiArrayIndex)assignments(s)) += columnVector(data, s);
            clusterSizes[(MultiArrayIndex)assignments(s)]++;
        }
        for(MultiArrayIndex k=0; k<clusterCount; ++k)
        {
            if(clusterSizes[k])
                columnVector(clusterCenters, k) *= 1.0/clusterSizes(k);
        }
    }
}

}

#endif // VIGRA_KMEANS_HXX

