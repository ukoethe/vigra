/************************************************************************/
/*                                                                      */
/*    Copyright 2008-2011 by Michael Hanselmann and Ullrich Koethe      */
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


#ifndef VIGRA_UNSUPERVISED_DECOMPOSITION_HXX
#define VIGRA_UNSUPERVISED_DECOMPOSITION_HXX

#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <list>
#include <numeric>
#include "mathutil.hxx"
#include "array_vector.hxx"
#include "sized_int.hxx"
#include "matrix.hxx"
#include "random.hxx"
#include "functorexpression.hxx"
#include "time.h"

namespace vigra
{

/** \addtogroup Unsupervised_Decomposition Unsupervised Decomposition

    This module provides unsupervised decomposition 
**/
//@{


/*****************************************************************/
/*                                                               */
/*         probabilistic latent semantic analysis (pLSA)         */
/*         see T Hofmann, Probabilistic Latent Semantic          */
/*         Indexing for details                                  */
/*                                                               */
/*****************************************************************/

class PLSAOptions
{
  public:
        /** Initialize all options with default values.
        */
    PLSAOptions()
    : numComponents(1),
      max_iterations(50),
      min_rel_gain(1e-4)
    {}

        /** Number of components.

            The number of components into which the data matrix is decomposed.<br>
        */
    PLSAOptions & numberOfComponents(unsigned int n)
    {
        vigra_precondition(n >= 1,
            "PLSAOptions::numberOfComponents(): number must be a positive integer."); // we could also check that the number of components is lower than the dimensions of the data matrix (maybe still unkown!)
        numComponents = n;
        return *this;
    }

        /** Maximum number of iterations which is performed by the pLSA algorithm.
        */
    PLSAOptions & maximumNumberOfIterations(unsigned int n)
    {
        vigra_precondition(n >= 1,
            "PLSAOptions::maximumNumberOfIterations(): number must be a positive integer.");
        max_iterations = n;
        return *this;
    }

        /** Minimum relative gain which is required for the algorithm to continue the iterations.
        */
    PLSAOptions & minimumRelativeGain(double g)
    {
        vigra_precondition(g >= 0.0,
            "PLSAOptions::minimumRelativeGain(): number must be positive or zero.");
        min_rel_gain = g;
        return *this;
    }
	
    double min_rel_gain;
    int max_iterations, numComponents;
};

class PLSA
{
  public:
    PLSAOptions options_;

  public:

    PLSA(PLSAOptions const & options = PLSAOptions())
    : options_(options)
    {
    }

		/** each component represents a hidden/latent "topic" in the data
			features is a matrix of dimension NUMFEATURESxNUMVOXELS
		    FZ is a matrix of dimension NUMFEATURESxNUMCOMPONENTS - the features/components matrix of typical features for each component type
		    ZV is a matrix of dimension NUMCOMPONENTSxNUMVOXELS - the components/voxels matrix of component type mixtures per voxel
		*/
    template <class U, class C, class Random>
    void decompose(MultiArrayView<2, U, C> const & features, MultiArrayView<2, U, C> & FZ, MultiArrayView<2, U, C> & ZV, Random const& random);

    template <class U, class C>
    void decompose(MultiArrayView<2, U, C> const & features, MultiArrayView<2, U, C> & FZ, MultiArrayView<2, U, C> & ZV)
    {
		RandomNumberGenerator<> generator(RandomSeed);
        return decompose(features, FZ, ZV, generator);
    }

		/** normalize the columns such that the sum of each column equals 1.
		*/
	template <class U>
	inline Matrix<U> normalizeColumns(Matrix<U> const & features)
	{
		double eps = 1/NumericTraits<U>::max();
		int rows = rowCount(features);
		int cols = columnCount(features);
		Matrix<U> result(rows, cols);
		Matrix<U> colSums = features.sum(0); //columnSums(features);
	    for (MultiArrayIndex j = 0; j < cols; j++)
	    {
		    for (MultiArrayIndex i = 0; i < rows; i++)
		    {
				result(i, j) = features(i, j)/(colSums(0, j) + (U)eps);
		    }
		}
        return result;
		//Matrix<U> ones(rows, 1, 1);
		//return pdiv(features, divisor + eps); //pdiv(features, ones * colSums + eps);
	}

};

template <class U, class C, class Random>
void
PLSA::decompose(MultiArrayView<2, U, C> const & features,
                                             MultiArrayView<2, U, C> & FZ, 
											 MultiArrayView<2, U, C> & ZV, 
											 Random const& random)
{
    Matrix<U> fz = FZ; // since we need the matrix multiplication operator later on, we have to work with the matrix class. ToDo: avoid copying!
    Matrix<U> zv = ZV;
	Matrix<U> feats = features;

	int numFeatures = rowCount(feats);
    int numVoxels = columnCount(feats);
	int numComponents = options_.numComponents;
    vigra_precondition((int)(numFeatures) >= (int)numComponents,
      "PLSA::decompose(): The number of features has to be larger or equal to the number of components in which the feature matrix is decomposed.");
    vigra_precondition(((int)(columnCount(FZ)) == (int)numComponents) && ((int)(rowCount(FZ)) == (int)numFeatures),
      "PLSA::decompose(): The output matrix FZ has to be of dimension NUMFEATURESxNUMCOMPONENTS.");
    vigra_precondition(((int)(columnCount(ZV)) == (int)numVoxels) && ((int)(rowCount(ZV)) == (int)numComponents),
      "PLSA::decompose(): The output matrix ZV has to be of dimension NUMCOMPONENTSxNUMVOXELS.");

	// random initialization of result matrices, subsequent normalization
    UniformRandomFunctor<Random> randf(random);
	for (unsigned int i = 0; i < (unsigned int)(numFeatures*numComponents); ++i)
		fz.data () [i] = (U)randf();
	for (unsigned int i = 0; i < (unsigned int)(numComponents*numVoxels); ++i)
		zv.data () [i] = 1;
	fz = normalizeColumns(fz);
	zv = normalizeColumns(zv);

	// init vars
	double eps = 1/NumericTraits<U>::max(); // epsilon > 0
    double lastChange = NumericTraits<U>::max(); // infinity
    double err = 0;
	double err_old;
    int iteration = 0;

    // expectation maximization (EM) algorithm
    Matrix<U> voxelSums = feats.sum(0);
	Matrix<U> fzv = fz*zv;
	Matrix<U> factor;
	Matrix<U> model;
	Matrix<U> ones(numFeatures, 1, 1);
	//clock_t start, finish;
	while(iteration < options_.max_iterations && (lastChange > options_.min_rel_gain))
	{
		if(iteration%25 == 0)
		{
			std::cout << "iteration: " << iteration << std::endl;
			//std::cout << "last relative change: " << lastChange << std::endl;
		}
		factor = pdiv(feats, fzv + (U)eps);
		//start = clock();
		zv *= (fz.transpose() * factor); //zv = (pmul(zv, (fz.transpose() * factor)));
		fz *= (factor * zv.transpose()); //fz = (pmul(fz, (factor * zv.transpose())));
		//finish = clock();
		//std::cout << "mult: " << finish-start << std::endl;
		//start = clock();
		zv = normalizeColumns(zv);
		fz = normalizeColumns(fz);
		//finish = clock();
		//std::cout << "norm: " << finish-start << std::endl;
        fzv = fz*zv; // pre-calculate for next iteration

		// check relative change in least squares model fit
        model = (pmul((ones * voxelSums), fzv));
        err_old = err;
        err = (feats - model).squaredNorm();
		//std::cout << "error: " << err << std::endl;
        lastChange = abs((err-err_old) / (U)(err + eps));
		//std::cout << "lastChange: " << lastChange << std::endl;
		 
        iteration += 1;
	}
	FZ = fz;
	ZV = zv;
	std::cout << "Terminated after " << iteration << " iterations." << std::endl;
	std::cout << "Last relative change was " << lastChange << "." << std::endl;
}

} // namespace vigra


#endif // VIGRA_UNSUPERVISED_DECOMPOSITION_HXX

