/************************************************************************/
/*                                                                      */
/*       Copyright 2011 by Ullrich Koethe and Michael Hanselmann        */
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
/*    The above copyrigfht notice and this permission notice shall be    */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyunsupervised_PyArray_API
// #define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/unsupervised_decomposition.hxx>
#include <set>
#include <cmath>
#include <memory>
#include <boost/python.hpp>

namespace python = boost::python;
namespace vigra
{

PLSA*
pythonConstructPLSA(int numComponents,
                            int max_iterations,
                            float min_relative_gain)


{
    PLSAOptions options;
    options .numberOfComponents(numComponents)
            .maximumNumberOfIterations(max_iterations)
            .minimumRelativeGain(min_relative_gain);


    PLSA* plsa = new PLSA(options);

    return plsa;
}

template<class U>
python::tuple
pythonDecompose(PLSA & plsa, 
                        NumpyArray<2,U> features)
{
    NumpyArray<2, U> FZ(MultiArrayShape<2>::type(features.shape(0),
                                                           plsa.options_.numComponents)); 
    NumpyArray<2, U> ZV(MultiArrayShape<2>::type(plsa.options_.numComponents, 
		                                                   features.shape(1))); 
    PyAllowThreads _pythread;
	plsa.decompose(features, FZ, ZV);

    return python::make_tuple(FZ, ZV);
}


void definePLSA()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);

    class_<PLSA> plsaclass("PLSA",python::no_init);

    plsaclass
        .def("__init__",python::make_constructor(registerConverters(&pythonConstructPLSA),
                                                 boost::python::default_call_policies(),
                                                 ( arg("numComponents")=1,
                                                   arg("max_iterations")= 100,
                                                   arg("min_relative_gain")=0.001f)),
             "Constructor::\n\n"
             "  PLSA(numComponents = 1, max_iterations=100, min_rel_gain=1e-4,\n")
			 // usage from python:  import vigra as v
			 //						plsa = v.unsupervised.PLSA(3)
			 //						fz, zv = plsa.decompose(data)
        .def("decompose",
             registerConverters(&pythonDecompose<float>),
             (arg("features")),
			 "Decomposes the feature matrix 'features' of size NUMFEATURESxNUMVOXELS. Returns two matrices: The former is of size NUMFEATURESxNUMCOMPONENTS and gives the feature to component distribution, the latter is of size NUMCOMPONENTSxNUMVOXELS and gives the component to voxel distribution. \n")
        ;
}

} // namespace vigra


using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(unsupervised)
{
    import_vigranumpy();
    definePLSA();
}


