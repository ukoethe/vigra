
/************************************************************************/
/*                                                                      */
/*                 Copyright 2010 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpygeometry_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/polygon.hxx>
#include <vigra/array_vector.hxx>
#include <vigra/tinyvector.hxx>

namespace python = boost::python;

namespace vigra
{

template < class Coordinate>
NumpyAnyArray
pyconvexHull(NumpyArray<1, TinyVector<Coordinate, 2>, UnstridedArrayTag > points)
{
    PyAllowThreads _pythread;
    
    ArrayVector<TinyVector<Coordinate, 2> > hull;
	convexHull(ArrayVectorView<TinyVector<Coordinate, 2> >(points.shape(0), points.data()), hull);

	NumpyArray<1, TinyVector<Coordinate, 2> > result(MultiArrayShape<1>::type(hull.size()));
	for(MultiArrayIndex i = 0; i < result.shape(0); ++i)
		result(i) = hull[i];

	return result;
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pythonConvexHull, pyconvexHull)

void defineGeometry()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    multidef("convexHull",
         pythonConvexHull<double, float, Int32>(),
         args("points"),
        "Compute the convex hull of a point set.\n"
        "\n"
        "For details see convexHull_ in the vigra C++ documentation.\n\n");
}

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(geometry)
{
    import_vigranumpy();
    defineGeometry();
}
