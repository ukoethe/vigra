/************************************************************************/
/*                                                                      */
/*                 Copyright 2011 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyoptimization_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/regression.hxx>

namespace python = boost::python;

namespace vigra
{

template <class T >
NumpyAnyArray
pythonLeastSquares(NumpyArray<2, T> A, NumpyArray<2, T> b)
{
    NumpyArray<2, T, UnstridedArrayTag> res(Shape2(A.shape(1), 1));
    
    {
        PyAllowThreads _pythread;
        
        leastSquares(A, b, res);
    }
    return res;
}

template <class T >
NumpyAnyArray
pythonNonnegativeLeastSquares(NumpyArray<2, T> A, NumpyArray<2, T> b)
{
    NumpyArray<2, T, UnstridedArrayTag> res(Shape2(A.shape(1), 1));
    
    {
        PyAllowThreads _pythread;
        
        nonnegativeLeastSquares(A, b, res);
    }
    return res;
}

template <class T >
NumpyAnyArray
pythonRidgeRegression(NumpyArray<2, T> A, NumpyArray<2, T> b, double lambda)
{
    NumpyArray<2, T, UnstridedArrayTag> res(Shape2(A.shape(1), 1));
    
    {
        PyAllowThreads _pythread;
        
        ridgeRegression(A, b, res, lambda);
    }
    return res;
}

template <class T >
python::tuple
pythonlassoRegression(NumpyArray<2, T> A, NumpyArray<2, T> b,
                      bool nonNegative,
                      bool lsqSolutions,
                      bool lassoSolutions,
                      unsigned int maxSolutionCount)
{
    vigra_precondition(lsqSolutions || lassoSolutions,
        "lassoRegression(): At least one of 'lsqSolutions' and 'lassoSolutions' must be 'True'.");
    
    ArrayVector<Matrix<double> > lasso_solutions, lsq_solutions;
    ArrayVector<ArrayVector<MultiArrayIndex> > activeSets;
    unsigned int numSolutions = 0;
    
    {
        PyAllowThreads _pythread;
        
        ArrayVector<Matrix<double> > * plasso = lassoSolutions
                                                     ? &lasso_solutions
                                                     : 0,
                                     * plsq = lsqSolutions
                                                     ? &lsq_solutions
                                                     : 0;
                                                            
        LeastAngleRegressionOptions options;
        if(nonNegative)
            options.nnlasso();
        else
            options.lasso();
        options.maxSolutionCount(maxSolutionCount);
    
        numSolutions = 
            linalg::detail::leastAngleRegressionImpl(A, b, activeSets, plasso, plsq, options);
    }
    
    python::list pyActiveSets;
    for(unsigned int k=0; k<numSolutions; ++k)
        pyActiveSets.append(python::object(activeSets[k]));
    python::list pyLassoSolutions;
    if(lassoSolutions)
    {
        for(unsigned int k=0; k<numSolutions; ++k)
        {
            NumpyArray<2, double, UnstridedArrayTag> sol(Shape2(A.shape(1), 1));
            for(unsigned int m=0; m<activeSets[k].size(); ++m)
            {
                sol(activeSets[k][m], 0) = lasso_solutions[k](m,0);
            }
            pyLassoSolutions.append(python::object(sol));
        }
    }
    python::list pyLsqSolutions;
    if(lsqSolutions)
    {
        for(unsigned int k=0; k<numSolutions; ++k)
        {
            NumpyArray<2, double, UnstridedArrayTag> sol(Shape2(A.shape(1), 1));
            for(unsigned int m=0; m<activeSets[k].size(); ++m)
            {
                sol(activeSets[k][m], 0) = lsq_solutions[k](m,0);
            }
            pyLsqSolutions.append(python::object(sol));
        }
    }
        
    if(lsqSolutions)
    {
        if(lassoSolutions)
            return python::make_tuple(numSolutions, pyActiveSets, pyLsqSolutions, pyLassoSolutions);
        else
            return python::make_tuple(numSolutions, pyActiveSets, pyLsqSolutions, python::object());
    }
    else
    {
        return python::make_tuple(numSolutions, pyActiveSets, python::object(), pyLassoSolutions);
    }
}


void defineOptimization()
{
    using namespace python;

    docstring_options doc_options(true, true, false);
    
    NumpyArrayConverter<NumpyArray<2, double, UnstridedArrayTag> >();

    def("leastSquares", registerConverters(&pythonLeastSquares<double>),
        (arg("A"), arg("b")),
        "Perform plain linear regression.\n"
        "\n"
        "For details see leastSquares_ in the vigra C++ documentation.\n\n");

    def("nonnegativeLeastSquares", registerConverters(&pythonNonnegativeLeastSquares<double>),
        (arg("A"), arg("b")),
        "Perform linear regression where the solution is constrained to be non-negative.\n"
        "\n"
        "For details see nonnegativeLeastSquares_ in the vigra C++ documentation.\n\n");

    def("ridgeRegression", registerConverters(&pythonRidgeRegression<double>),
        (arg("A"), arg("b"), arg("lambda")),
        "Perform linear regression with L2 regularization.\n"
        "\n"
        "'lambda' is the regularization parameter - the larger it is, the more\n"
        "biased towards zero the solution will become.\n"
        "\n"
        "For details see ridgeRegression_ in the vigra C++ documentation.\n\n");

    def("lassoRegression", registerConverters(&pythonlassoRegression<double>),
        (arg("A"), arg("b"), 
         arg("nonNegative")=false, arg("lsq")=true, arg("lasso")=false, arg("maxSolutionCount")=0),
        "Perform linear regression with L1 regularization.\n"
        "\n"
        "If 'nonNegative' is 'True', the solution will be constrained to non-negative\n"
        "values, otherwise values may have arbitrary sign (the default).\n"
        "If 'lsq' is 'True', the algorithm will return the least squares solution\n"
        "for each active set. If 'lasso' is 'True', the LASSO solution will be returned\n"
        "for each active set. Both may be 'True' simultaneously.\n"
        "If 'maxSolutionCount' is non-zero, atr most so many active sets will\n"
        "be computed.\n"
        "\n"
        "The algorithm returns a tuple::\n\n"
        "   (numActiveSets, activeSets, lsqSolutions, lassoSolutions)\n\n"
        "where 'numActiveSets' specifies how many active sets have been computed,\n"
        "'activeSets' is the list of all active sets (ordered by decreasing regularization),\n"
        "and 'lsqSolutions' and 'lassoSolutions' are lists of the corresponding solutions\n"
        "for each active set ('lsqSolutions' and 'lassoSolutions' will be 'None' when\n"
        "the corresponding function argument was 'False'). An active set is a list of\n"
        "indices of all variables whose values are non-zero in the corresponding\n"
        "solution.\n"
        "\n"
        "For details see leastAngleRegression_ in the vigra C++ documentation.\n\n");
}

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(optimization)
{
    import_vigranumpy();
    defineOptimization();
}
