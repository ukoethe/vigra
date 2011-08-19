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

#include <numeric>
#include "mathutil.hxx"
#include "matrix.hxx"
#include "singular_value_decomposition.hxx"
#include "random.hxx"

namespace vigra
{

/** \addtogroup Unsupervised_Decomposition Unsupervised Decomposition

    Unsupervised matrix decomposition methods.
**/
//@{

/*****************************************************************/
/*                                                               */
/*              principle component analysis (PCA)               */
/*                                                               */
/*****************************************************************/

   /** \brief Decompose a matrix according to the PCA algorithm. 

        This function implements the PCA algorithm (principle component analysis).

        \arg features must be a matrix with shape <tt>(numFeatures * numSamples)</tt>, which is
        decomposed into the matrices 
        \arg fz with shape <tt>(numFeatures * numComponents)</tt> and
        \arg zv with shape <tt>(numComponents * numSamples)</tt>
        
        such that
        \f[
            \mathrm{features} \approx \mathrm{fz} * \mathrm{zv}
        \f]
        (this formula requires that the features have been center around the mean by  
        \ref <tt>prepareRows(features, features, ZeroMean)</tt>).
       
        The shape parameter <tt>numComponents</tt> determines the complexity of 
        the decomposition model and therefore the approximation quality (if
        <tt>numComponents == numFeatures</tt>, the representation becomes exact). 
        Intuitively, <tt>fz</tt> is a projection matrix from the reduced space
        into the original space, and <tt>zv</tt> is  the reduced representation
        of the data, using just <tt>numComponents</tt> features.

        <b>Declaration:</b>
        
        <b>\#include</b> \<vigra/unsupervised_decomposition.hxx\>

        \code
        namespace vigra {
        
            template <class U, class C1, class C2, class C3>
            void
            principleComponents(MultiArrayView<2, U, C1> const & features,
                                MultiArrayView<2, U, C2> fz, 
                                MultiArrayView<2, U, C3> zv);
        }
        \endcode
        
        <b>Usage:</b>
        \code
        Matrix<double> data(numFeatures, numSamples);
        ... // fill the imput matrix
        
        int numComponents = 3;
        Matrix<double> fz(numFeatures, numComponents),
                       zv(numComponents, numSamples);
                       
        // center the data
        prepareRows(data, data, ZeroMean);
        
        // compute the reduced representation
        principleComponents(data, fz, zv);
        
        Matrix<double> model = fz*zv;
        double meanSquaredError = squaredNorm(data - model) / numSamples;
        \endcode
   */
template <class U, class C1, class C2, class C3>
void
principleComponents(MultiArrayView<2, U, C1> const & features,
                    MultiArrayView<2, U, C2> fz, 
                    MultiArrayView<2, U, C3> zv)
{
    using namespace linalg; // activate matrix multiplication and arithmetic functions

    int numFeatures = rowCount(features);
    int numSamples = columnCount(features);
    int numComponents = columnCount(fz);
    vigra_precondition(numSamples >= numFeatures,
      "principleComponents(): The number of samples has to be larger than the number of features.");
    vigra_precondition(numFeatures >= numComponents && numComponents >= 1,
      "principleComponents(): The number of features has to be larger or equal to the number of components in which the feature matrix is decomposed.");
    vigra_precondition(rowCount(fz) == numFeatures,
      "principleComponents(): The output matrix fz has to be of dimension numFeatures*numComponents.");
    vigra_precondition(columnCount(zv) == numSamples && rowCount(zv) == numComponents,
      "principleComponents(): The output matrix zv has to be of dimension numComponents*numSamples.");

    Matrix<U> U(numSamples, numFeatures), S(numFeatures, 1), V(numFeatures, numFeatures);
    singularValueDecomposition(features.transpose(), U, S, V);
    
    for(int k=0; k<numComponents; ++k)
    {
        rowVector(zv, k) = columnVector(U, k).transpose() * S(k, 0);
        columnVector(fz, k) = columnVector(V, k);
    }
}

/*****************************************************************/
/*                                                               */
/*         probabilistic latent semantic analysis (pLSA)         */
/*         see T Hofmann, Probabilistic Latent Semantic          */
/*         Indexing for details                                  */
/*                                                               */
/*****************************************************************/

   /** \brief Option object for the \ref pLSA algorithm. 
   */
class PLSAOptions
{
  public:
        /** Initialize all options with default values.
        */
    PLSAOptions()
    : min_rel_gain(1e-4),
      max_iterations(50),
      normalized_component_weights(true)
    {}

        /** Maximum number of iterations which is performed by the pLSA algorithm.

            default: 50
        */
    PLSAOptions & maximumNumberOfIterations(unsigned int n)
    {
        vigra_precondition(n >= 1,
            "PLSAOptions::maximumNumberOfIterations(): number must be a positive integer.");
        max_iterations = n;
        return *this;
    }

        /** Minimum relative gain which is required for the algorithm to continue the iterations.

            default: 1e-4
        */
    PLSAOptions & minimumRelativeGain(double g)
    {
        vigra_precondition(g >= 0.0,
            "PLSAOptions::minimumRelativeGain(): number must be positive or zero.");
        min_rel_gain = g;
        return *this;
    }
    
        /** Normalize the entries of the zv result array.
        
            If true, the columns of zv sum to one. Otherwise, they have the same
            column sum as the original feature matrix.
            
            default: true
        */
    PLSAOptions & normalizedComponentWeights(bool v = true)
    {
        normalized_component_weights = v;
        return *this;
    }

    double min_rel_gain;
    int max_iterations;
    bool normalized_component_weights;
};

   /** \brief Decompose a matrix according to the pLSA algorithm. 

        This function implements the pLSA algorithm (probabilistic latent semantic analysis) 
        proposed in

        T. Hofmann: <a href="http://www.cs.brown.edu/people/th/papers/Hofmann-UAI99.pdf">
        <i>"Probabilistic Latent Semantic Analysis"</i></a>,
        in: UAI'99, Proc. 15th Conf. on Uncertainty in Artificial Intelligence,
        pp. 289-296, Morgan Kaufmann, 1999

        \arg features must be a matrix with shape <tt>(numFeatures * numSamples)</tt>, which is
        decomposed into the matrices 
        \arg fz with shape <tt>(numFeatures * numComponents)</tt> and
        \arg zv with shape <tt>(numComponents * numSamples)</tt>
        
        such that
        \f[
            \mathrm{features} \approx \mathrm{fz} * \mathrm{zv}
        \f]
        (this formula applies when pLSA is called with 
        <tt>PLSAOptions.normalizedComponentWeights(false)</tt>. Otherwise, you must 
        normalize the features by calling \ref <tt>prepareColumns(features, features, UnitSum)</tt> to make the formula hold).
       
        The shape parameter <tt>numComponents</tt> determines the complexity of 
        the decomposition model and therefore the approximation quality. 
        Intuitively, features are a set of words, and the samples a set of 
        documents. The entries of the <tt>features</tt> matrix denote the relative 
        frequency of the words in each document. The components represents a 
        (presumably small) set of topics. The matrix <tt>fz</tt> encodes the 
        relative frequency of words in the different topics, and the matrix  
        <tt>zv</tt> encodes to what extend each topic explains the content of each 
        document.

        The option object determines the iteration termination conditions and the ouput
        normalization. In addition, you may pass a random number generator to pLSA()
        which is used to create the initial solution.

        <b>Declarations:</b>
        
        <b>\#include</b> \<vigra/unsupervised_decomposition.hxx\>

        \code
        namespace vigra {
        
            template <class U, class C1, class C2, class C3, class Random>
            void
            pLSA(MultiArrayView<2, U, C1> const & features,
                 MultiArrayView<2, U, C2> & fz, 
                 MultiArrayView<2, U, C3> & zv,
                 Random const& random,
                 PLSAOptions const & options = PLSAOptions());
                 
            template <class U, class C1, class C2, class C3>
            void
            pLSA(MultiArrayView<2, U, C1> const & features, 
                 MultiArrayView<2, U, C2> & fz, 
                 MultiArrayView<2, U, C3> & zv,
                 PLSAOptions const & options = PLSAOptions());
        }
        \endcode
        
        <b>Usage:</b>
        \code
        Matrix<double> words(numWords, numDocuments);
        ... // fill the imput matrix
        
        int numTopics = 3;
        Matrix<double> fz(numWords, numTopics),
                       zv(numTopics, numDocuments);
                       
        pLSA(words, fz, zv, PLSAOptions().normalizedComponentWeights(false));
        
        Matrix<double> model = fz*zv;
        double meanSquaredError = (words - model).squaredNorm() / numDocuments;
        \endcode
   */
doxygen_overloaded_function(template <...> void pLSA)


template <class U, class C1, class C2, class C3, class Random>
void
pLSA(MultiArrayView<2, U, C1> const & features,
     MultiArrayView<2, U, C2> fz, 
     MultiArrayView<2, U, C3> zv,
     Random const& random,
     PLSAOptions const & options = PLSAOptions())
{
    using namespace linalg; // activate matrix multiplication and arithmetic functions

    int numFeatures = rowCount(features);
    int numSamples = columnCount(features);
    int numComponents = columnCount(fz);
    vigra_precondition(numFeatures >= numComponents && numComponents >= 1,
      "pLSA(): The number of features has to be larger or equal to the number of components in which the feature matrix is decomposed.");
    vigra_precondition(rowCount(fz) == numFeatures,
      "pLSA(): The output matrix fz has to be of dimension numFeatures*numComponents.");
    vigra_precondition(columnCount(zv) == numSamples && rowCount(zv) == numComponents,
      "pLSA(): The output matrix zv has to be of dimension numComponents*numSamples.");

    // random initialization of result matrices, subsequent normalization
    UniformRandomFunctor<Random> randf(random);
    initMultiArray(destMultiArrayRange(fz), randf);
    initMultiArray(destMultiArrayRange(zv), randf);
    prepareColumns(fz, fz, UnitSum);
    prepareColumns(zv, zv, UnitSum);

    // init vars
    double eps = 1.0/NumericTraits<U>::max(); // epsilon > 0
    double lastChange = NumericTraits<U>::max(); // infinity
    double err = 0;
    double err_old;
    int iteration = 0;

    // expectation maximization (EM) algorithm
    Matrix<U> columnSums(1, numSamples);
    features.sum(columnSums);
    Matrix<U> expandedSums = ones<U>(numFeatures, 1) * columnSums;
    
    while(iteration < options.max_iterations && (lastChange > options.min_rel_gain))
    {
        Matrix<U> fzv = fz*zv;
        
        //if(iteration%25 == 0)
        //{
            //std::cout << "iteration: " << iteration << std::endl;
            //std::cout << "last relative change: " << lastChange << std::endl;
        //}

        Matrix<U> factor = features / pointWise(fzv + (U)eps);
        zv *= (fz.transpose() * factor);
        fz *= (factor * zv.transpose());
        prepareColumns(fz, fz, UnitSum);
        prepareColumns(zv, zv, UnitSum);

        // check relative change in least squares model fit
        Matrix<U> model = expandedSums * pointWise(fzv);
        err_old = err;
        err = (features - model).squaredNorm();
        //std::cout << "error: " << err << std::endl;
        lastChange = abs((err-err_old) / (U)(err + eps));
        //std::cout << "lastChange: " << lastChange << std::endl;
         
        iteration += 1;
    }
    //std::cout << "Terminated after " << iteration << " iterations." << std::endl;
    //std::cout << "Last relative change was " << lastChange << "." << std::endl;
    
    if(!options.normalized_component_weights)
    {
        // undo the normalization
        for(int k=0; k<numSamples; ++k)
            columnVector(zv, k) *= columnSums(0, k);
    }
}

template <class U, class C1, class C2, class C3>
inline void
pLSA(MultiArrayView<2, U, C1> const & features, 
     MultiArrayView<2, U, C2> & fz, 
     MultiArrayView<2, U, C3> & zv,
     PLSAOptions const & options = PLSAOptions())
{
    RandomNumberGenerator<> generator(RandomSeed);
    pLSA(features, fz, zv, generator, options);
}

//@}

} // namespace vigra


#endif // VIGRA_UNSUPERVISED_DECOMPOSITION_HXX

