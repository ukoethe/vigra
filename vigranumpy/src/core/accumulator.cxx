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
#include <map>

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

template <class T, int N>
struct InspectResultType<TinyVector<T, N> >
{
    typedef TinyVector<typename PromoteTraits<T, double>::Promote, N> type;
};

template <class T>
struct InspectResultType<Singleband<T> >
{
    typedef typename PromoteTraits<T, double>::Promote type;
};

typedef std::map<std::string, std::string> AliasMap;

AliasMap createAliasMap()
{
    AliasMap res;
    res["DivideByCount<Central<PowerSum<2> > >"] = "Variance";
    res["DivideByCount<Principal<PowerSum<2> > >"] = "Principal<Variance>";
    res["DivideByCount<FlatScatterMatrix>"] = "Covariance";
    res["DivideByCount<PowerSum<1> >"] = "Mean";
    res["PowerSum<1>"] = "Sum";
    res["PowerSum<0>"] = "Count";
    res["StandardQuantiles<AutoRangeHistogram<100 > >"] = "Quantiles";
    return res;
}

static AliasMap const & tagToAlias()
{
    static const AliasMap m = createAliasMap();
    return m;
}

AliasMap createInverseAliasMap()
{
    AliasMap res;
    for(AliasMap::const_iterator k = tagToAlias().begin(); k != tagToAlias().end(); ++k)
        res[normalizeString(k->second)] = normalizeString(k->first);
    return res;
}

static AliasMap const & aliasToTag()
{
    static const AliasMap m = createInverseAliasMap();
    return m;
}

static std::string const & resolveAlias(std::string const & n)
{
    AliasMap::const_iterator k = aliasToTag().find(n);
    if(k == aliasToTag().end())
        return n;
    else
        return k->second;
}

static std::string const & createAlias(std::string const & n)
{
    AliasMap::const_iterator k = tagToAlias().find(n);
    if(k == tagToAlias().end())
        return n;
    else
        return k->second;
}

template <class PixelType, class Accumulators>
struct PythonAccumulator
: public DynamicAccumulatorChain<PixelType, Accumulators>
{
    typedef DynamicAccumulatorChain<PixelType, Accumulators> BaseType;
    typedef typename BaseType::AccumulatorTags AccumulatorTags;
    
    void activate(std::string tag)
    {
        bool found = detail::ApplyVisitorToTag<AccumulatorTags>::exec((BaseType &)*this, 
                                     resolveAlias(normalizeString(tag)), detail::ActivateTagVisitor());
        vigra_precondition(found, std::string("PythonAccumulator::activate(): Tag '") + tag + "' not found.");
    }
    
    python::list namesImpl(bool activeOnly) const
    {
        ArrayVector<std::string> a = BaseType::namesImpl(activeOnly);
        for(unsigned int k=0; k<a.size(); ++k)
            a[k] = createAlias(a[k]);

        std::sort(a.begin(), a.end());
        
        python::list result;
        for(unsigned int k=0; k<a.size(); ++k)
            result.append(python::object(a[k]));
        return result;
    }
    
    python::list activeNames() const
    {
        return namesImpl(true);
    }
    
    python::list names() const
    {
        return namesImpl(false);
    }
    
    python::object get(std::string tag)
    {
        GetTagVisitor v;
        
        bool found = detail::ApplyVisitorToTag<AccumulatorTags>::exec(*this, 
                                           resolveAlias(normalizeString(tag)), v);
        vigra_precondition(found, std::string("PythonAccumulator::get(): Tag '") + tag + "' not found.");
        return v.result;
    }
    
    void mergeImpl(PythonAccumulator const & o)
    {
        this->merge(o);
    }
};

template <class Accumulators, unsigned int ndim, class T>
PythonAccumulator<typename InspectResultType<T>::type, Accumulators> *
pythonInspect(NumpyArray<ndim, T> in, python::object tags)
{
    typedef PythonAccumulator<typename InspectResultType<T>::type, Accumulators> Accu;
    
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
        
        collectStatistics(in.begin(), in.end(), *res);
    }
    
    return res.release();
}

} // namespace acc1

template <class T, class Accumulators>
void definePythonAccumulator()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename acc1::InspectResultType<T>::type ResultType;
    
    typedef acc1::PythonAccumulator<ResultType, Accumulators> Accu;
    class_<Accu>("Accumulator", python::no_init)
        .def("__getitem__", &Accu::get)
        .def("activeNames", &Accu::activeNames)
        .def("names", &Accu::names)
        .def("merge", &Accu::mergeImpl)
        ;
    
    def("extractFeatures", &acc1::pythonInspect<Accumulators, 2, T>,
          (arg("image"), arg("tags") = ""),
          return_value_policy<manage_new_object>());
    
    def("extractFeatures", &acc1::pythonInspect<Accumulators, 3, T>,
          (arg("volume"), arg("tags") = ""),
          return_value_policy<manage_new_object>());
}

void defineAccumulators()
{
    using namespace python;
    using namespace vigra::acc1;

    docstring_options doc_options(true, true, false);

    definePythonAccumulator<Singleband<float>, 
                            Select<Count, Mean, Variance, Skewness, Kurtosis, Minimum, Maximum,
                                   StandardQuantiles<AutoRangeHistogram<100> > > >();
    
    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, Minimum, Maximum, 
                   Covariance, Principal<Variance>, Principal<Skewness>, Principal<Kurtosis> 
                   > VectorAccumulators;
    definePythonAccumulator<TinyVector<float, 2>, VectorAccumulators>();
    definePythonAccumulator<TinyVector<float, 3>, VectorAccumulators>();
    definePythonAccumulator<TinyVector<float, 4>, VectorAccumulators>();
}

} // namespace vigra
