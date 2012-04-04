/************************************************************************/
/*                                                                      */
/*            Copyright 2011-2012 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyanalysis_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/accumulator.hxx>
#include <vigra/timing.hxx>

namespace python = boost::python;

namespace vigra
{

namespace acc1 
{

struct GetTagVisitor
{
    mutable python::object result;
    
    template <class TAG, class Accu>
    void exec(Accu & a) const
    {
        result = python::object(get<TAG>(a));
    }
};

template <class T>
struct InspectResultType
{
    typedef typename PromoteTraits<T, double>::Promote type;
};

template <class T>
struct InspectResultType<Singleband<T> >
{
    typedef typename PromoteTraits<T, double>::Promote type;
};

template <class PixelType>
struct PythonAccumulator
: public DynamicAccumulatorChain<PixelType, Select<Count, Mean, Variance, Skewness, Kurtosis> >
{
    typedef typename PythonAccumulator::AccumulatorTags AccumulatorTags;
    
    python::object get(std::string tag)
    {
        GetTagVisitor v;
        
        bool found = detail::ApplyVisitorToTag<AccumulatorTags>::exec(*this, 
                             detail::resolveAlias(normalizeString(tag)), v);
        vigra_precondition(found, std::string("PythonAccumulator::get(): Tag '") + tag + "' not found.");
        return v.result;
    }
    
    void mergeImpl(PythonAccumulator const & o)
    {
        this->merge(o);
    }
};

template <unsigned int ndim, class T>
PythonAccumulator<typename InspectResultType<T>::type> *
pythonInspect(NumpyArray<ndim, T> in, python::object tags)
{
    typedef PythonAccumulator<typename InspectResultType<T>::type> Accu;
    
    std::auto_ptr<Accu> res(new Accu);
    
    if(python::len(tags) == 0)
    {
        res->activateAll();
    }
    else if(PyString_Check(tags.ptr()))
    {
        res->activate(python::extract<std::string>(tags)());
    }
    else
    {
        for(int k=0; k<python::len(tags); ++k)
        {
            res->activate(python::extract<std::string>(tags[k])());
        }
    }
    
    {
        PyAllowThreads _pythread;
        
        // USETICTOC;
        // TIC;
        collectStatistics(in.begin(), in.end(), *res);
        // TOC;
    }
    
    return res.release();
}

} // namespace acc1

template <class T>
void definePythonAccumulator()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename acc1::InspectResultType<T>::type ResultType;
    
    typedef acc1::PythonAccumulator<ResultType> Accu;
    class_<Accu>("Accumulator", python::no_init)
        .def("__getitem__", &Accu::get)
        .def("merge", &Accu::mergeImpl)
        ;
    
    def("extractFeatures", &acc1::pythonInspect<2, T>,
          return_value_policy<manage_new_object>());
}

void defineAccumulators()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    definePythonAccumulator<Singleband<float> >();
    definePythonAccumulator<TinyVector<float, 3> >();
}

} // namespace vigra
