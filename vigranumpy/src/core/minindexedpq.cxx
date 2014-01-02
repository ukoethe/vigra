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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyclustering_PyArray_API
#define NO_IMPORT_ARRAY

#include <boost/python.hpp>
#include <vigra/merge_graph/min_indexed_pq.hxx>
namespace python = boost::python;

namespace vigra
{



template<class T>
void defineMinIndexedPqT(const std::string name){
    using namespace python;
    docstring_options doc_options(true, true, false);

    typedef MinIndexedPQ<T>  MinIndexedPQType;

    class_<MinIndexedPQType>(name.c_str(),init<const size_t >())
	    .def("insert",&MinIndexedPQType::insert,
		    (
		    	arg("index"),
		    	arg("value")
		    ),
		    "insert a index with an value"
    	)
    	.def("changeValue",&MinIndexedPQType::changeValue,
		    (
		    	arg("index"),
		    	arg("value")
		    ),
		    "insert a index with an value"
    	)
    	.def("decreaseValue",&MinIndexedPQType::decreaseValue,
		    (
		    	arg("index"),
		    	arg("value")
		    ),
		    "insert a index with an value"
    	)
    	.def("increaseValue",&MinIndexedPQType::increaseValue,
		    (
		    	arg("index"),
		    	arg("value")
		    ),
		    "insert a index with an value"
    	)
    	.def("__len__",&MinIndexedPQType::size,"number of elements in PQ")
    	.def("isEmpty",&MinIndexedPQType::isEmpty,"check if pq is empty")

    ;
}


void defineMinIndexedPq(){
	defineMinIndexedPqT<float>("MinIndexedPQFloat");
}


} // namespace vigra

