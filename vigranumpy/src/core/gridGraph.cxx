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

//#define PY_ARRAY_UNIQUE_SYMBOL vigranumpygraphs_PyArray_API
#define NO_IMPORT_ARRAY

#include "export_graph_visitor.hxx"

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_gridgraph.hxx>
namespace python = boost::python;

namespace vigra{

    template<unsigned int DIM>
    void defineGridGraphT(const std::string & clsName){
        
        typedef GridGraph<DIM,boost::undirected_tag> GridGraph;

        python::class_<GridGraph>(clsName.c_str(),python::init< TinyVector<Int64,DIM> >())
        .def(LemonDirectedGraphCoreVisitor<
                GridGraph,
                python::return_value_policy<python::return_by_value>,
                python::return_value_policy<python::return_by_value>
            >(clsName)
        )
        ;

        
    }

    void defineGridGraph(){
        defineGridGraphT<2>("GridGraphUndirected2d");
    }

    void defineInvalid(){
        python::class_<lemon::Invalid>("Invalid",python::init<>())
        ;
    }

} 


