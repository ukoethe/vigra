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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyutilities_PyArray_API
//#define NO_IMPORT_ARRAY

#include <string>

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/priority_queue.hxx>
namespace python = boost::python;

namespace vigra{

    template<class PQ>
    void pyPush(
        PQ & pq, 
        const NumpyArray<1,UInt32> indices,
        const NumpyArray<1,float>  priorities
    ){
        for(size_t i=0;i<indices.shape(0);++i){
            pq.push(indices(i),priorities(i));
        }
    }





    template<class T,class COMP>
	void defineChangeablePriorityQueue(const std::string & clsName){

        typedef ChangeablePriorityQueue<T,COMP> PQ;

        python::class_<PQ>(clsName.c_str(),python::init<const size_t>())
        .def("push",            registerConverters(&pyPush<PQ>))
        .def("push",            &PQ::push)
        .def("pop",             &PQ::pop)
        .def("top",             &PQ::top)
        .def("topPriority",     &PQ::topPriority)
        .def("deleteItem",      &PQ::deleteItem)
        .def("__len__",         &PQ::size)
        .def("contains",    &PQ::contains)
        .def("__empty__",       &PQ::empty)
        ;

    }

} // namespace vigra

using namespace vigra;
using namespace boost::python;




BOOST_PYTHON_MODULE_INIT(utilities)
{
    import_vigranumpy();

    // all exporters needed for graph exporters (like lemon::INVALID)
    defineChangeablePriorityQueue<float,std::less<float> >("ChangeablePriorityQueueFloat32Min");
}
