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

#ifdef _MSC_VER
#pragma warning (disable: 4503)
#endif

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

struct GetTag_Visitor
{
    mutable python::object result;
    
    template <class T>
    void to_python(T const & t) const
    {
        result = python::object(t);
    }
    
    template <class T, int N>
    void to_python(TinyVector<T, N> const & t) const
    {
        NumpyArray<1, T> a = NumpyArray<1, T>(Shape1(N));
        for(int k=0; k<N; ++k)
            a(k) = t[k];
        result = python::object(a);
    }
    
    template <class TAG, class Accu>
    void exec(Accu & a) const
    {
        to_python(get<TAG>(a));
    }
};

struct GetArrayTag_Visitor
: public GetTag_Visitor
{
    template <class TAG, class T, class Accu>
    struct ToPythonArray
    {
        static python::object exec(Accu & a)
        {
            unsigned int n = a.regionCount();
            Shape1 s(n);
            NumpyArray<1, T> res(s);
            
            for(unsigned int k=0; k<n; ++k)
                res(k) = get<TAG>(a, k);
            return python::object(res);
        }
    };
    
    template <class TAG, class T, int N, class Accu>
    struct ToPythonArray<TAG, TinyVector<T, N>, Accu>
    {
        static python::object exec(Accu & a)
        {
            unsigned int n = a.regionCount();
            Shape2 s(n, N);
            NumpyArray<2, T> res(s);
            
            for(unsigned int k=0; k<n; ++k)
                for(int j=0; j<N; ++j)
                    res(k, j) = get<TAG>(a, k)[j];
            return python::object(res);
        }
    };
    
    template <class TAG, class T, class Accu>
    struct ToPythonArray<TAG, Matrix<T>, Accu>
    {
        static python::object exec(Accu & a)
        {
            unsigned int n = a.regionCount();
            Shape2 m = get<TAG>(a, 0).shape();
            Shape3 s(n, m[0], m[1]);
            NumpyArray<3, T> res(s);
            
            for(unsigned int k=0; k<n; ++k)
                for(int j=0; j<m[0]; ++j)
                    for(int i=0; i<m[1]; ++i)
                        res(k, j, i) = get<TAG>(a, k)(j,i);
            return python::object(res);
        }
    };
    
    template <class TAG, class T, class Accu>
    struct ToPythonArray<TAG, Error__Attempt_to_access_inactive_statistic<T>, Accu>
    {
        static python::object exec(Accu & a)
        {
            return python::object();
        }
    };
    
    template <class TAG, class T, int N, class Accu>
    struct ToPythonArray<TAG, std::pair<TinyVector<T,N>, Matrix<T> >, Accu>
    {
        static python::object exec(Accu & a)
        {
            return python::object();
        }
    };
    
    template <class TAG, class Accu>
    void exec(Accu & a) const
    {
        exec(a, (TAG *)0);
    }
    
    template <class Accu, class TAG>
    void exec(Accu & a, TAG *) const
    {
        this->result = ToPythonArray<TAG, typename LookupTag<TAG, Accu>::value_type, Accu>::exec(a);
    }
    
    template <class Accu, class TAG>
    void exec(Accu & a, Global<TAG> *) const
    {
        to_python(get<Global<TAG> >(a));
    }
    
};

template <class T>
struct AccumulatorValueTypeTraits
{
    // typedef typename PromoteTraits<T, double>::Promote type;
    typedef T type;
};

// template <class T, int N>
// struct AccumulatorValueTypeTraits<TinyVector<T, N> >
// {
    // typedef TinyVector<typename PromoteTraits<T, double>::Promote, N> type;
// };

template <class T>
struct AccumulatorValueTypeTraits<Singleband<T> >
{
    // typedef typename PromoteTraits<T, double>::Promote type;
    typedef T type;
};

typedef std::map<std::string, std::string> AliasMap;

AliasMap defineAliasMap()
{
    AliasMap res;
    res["DivideByCount<Central<PowerSum<2> > >"] = "Variance";
    res["DivideUnbiased<Central<PowerSum<2> > >"] = "UnbiasedVariance";
    res["DivideByCount<Principal<PowerSum<2> > >"] = "Principal<Variance>";
    res["DivideByCount<FlatScatterMatrix>"] = "Covariance";
    res["DivideByCount<PowerSum<1> >"] = "Mean";
    res["PowerSum<1>"] = "Sum";
    res["PowerSum<0>"] = "Count";
    res["StandardQuantiles<AutoRangeHistogram<100> >"] = "Quantiles";
    res["StandardQuantiles<GlobalRangeHistogram<100> >"] = "Quantiles";
    res["Coord<DivideByCount<PowerSum<1> > >"] = "GeometricCenter";
    res["Coord<RootDivideByCount<Principal<PowerSum<2> > > >"] = "PrincipalRadii";
    res["Coord<Principal<CoordinateSystem> >"] = "PrincipalCoordSystem";
    res["Weighted<Coord<DivideByCount<PowerSum<1> > > >"] = "CenterOfMass";
    res["Weighted<Coord<DivideByCount<Principal<PowerSum<2> > > > >"] = "MomentsOfInertia";
    res["Weighted<Coord<Principal<CoordinateSystem> > >"] = "CoordSystemOfInertia";
    return res;
}

AliasMap createTagToAlias(ArrayVector<std::string> const & names)
{
    static const AliasMap aliases = defineAliasMap();
    AliasMap res;
    for(unsigned int k=0; k<names.size(); ++k)
    {
        AliasMap::const_iterator a = aliases.find(names[k]);
        if(a == aliases.end())
            res[names[k]] = names[k];
        else
            res[names[k]] = a->second;
    }
    return res;   
}

AliasMap createAliasToTag(AliasMap const & tagToAlias)
{
    AliasMap res;
    for(AliasMap::const_iterator k = tagToAlias.begin(); k != tagToAlias.end(); ++k)
        res[normalizeString(k->second)] = normalizeString(k->first);
    return res;
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
                                     resolveAlias(tag), detail::ActivateTag_Visitor());
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
        GetTag_Visitor v;
        
        bool found = detail::ApplyVisitorToTag<AccumulatorTags>::exec(*this, resolveAlias(tag), v);
        vigra_precondition(found, std::string("PythonAccumulator::get(): Tag '") + tag + "' not found.");
        return v.result;
    }
    
    void merge(PythonAccumulator const & o)
    {
        BaseType::merge(o);
    }
    
    static std::string createAlias(std::string const & n)
    {
        AliasMap::const_iterator k = tagToAlias().find(n);
        if(k == tagToAlias().end())
            return n;
        else
            return k->second;
    }
    
    static std::string resolveAlias(std::string const & n)
    {
        AliasMap::const_iterator k = aliasToTag().find(normalizeString(n));
        if(k == aliasToTag().end())
            return n;
        else
            return k->second;
    }
    
    static AliasMap const & tagToAlias()
    {
        static const AliasMap a = createTagToAlias(tagNames());
        return a;   
    }
    
    static AliasMap const & aliasToTag()
    {
        static const AliasMap a = createAliasToTag(tagToAlias());
        return a;   
    }
};

template <class Accumulators, unsigned int ndim, class T>
PythonAccumulator<typename AccumulatorValueTypeTraits<T>::type, Accumulators> *
pythonInspect(NumpyArray<ndim, T> in, python::object tags)
{
    typedef PythonAccumulator<typename AccumulatorValueTypeTraits<T>::type, Accumulators> Accu;
    
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

template <class PixelType, class Accumulators>
struct PythonAccumulatorArray
: public DynamicAccumulatorChainArray<PixelType, Accumulators>
{
    typedef DynamicAccumulatorChainArray<PixelType, Accumulators> BaseType;
    typedef typename BaseType::AccumulatorTags AccumulatorTags;
    
    void activate(std::string tag)
    {
        bool found = detail::ApplyVisitorToTag<AccumulatorTags>::exec((BaseType &)*this, 
                                             resolveAlias(tag), detail::ActivateTag_Visitor());
        vigra_precondition(found, std::string("PythonAccumulatorArray::activate(): Tag '") + tag + "' not found.");
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
        ArrayVector<std::string> a;
        for(AliasMap::const_iterator k = tagToAlias().begin(); k != tagToAlias().end(); ++k)
            a.push_back(k->second);
        std::sort(a.begin(), a.end());
        
        python::list result;
        for(unsigned int k=0; k<a.size(); ++k)
            result.append(python::object(a[k]));
        return result;
    }
    
    static std::string createAlias(std::string const & n)
    {
        AliasMap::const_iterator k = tagToAlias().find(n);
        if(k == tagToAlias().end())
            return n;
        else
            return k->second;
    }
    
    static std::string resolveAlias(std::string const & n)
    {
        AliasMap::const_iterator k = aliasToTag().find(normalizeString(n));
        if(k == aliasToTag().end())
            return n;
        else
            return k->second;
    }
    
    static AliasMap const & tagToAlias()
    {
        static const AliasMap a = createTagToAlias(tagNames());
        return a;   
    }
    
    static AliasMap const & aliasToTag()
    {
        static const AliasMap a = createAliasToTag(tagToAlias());
        return a;   
    }
    
    python::object get(std::string tag)
    {
        GetArrayTag_Visitor v;
        
        bool found = detail::ApplyVisitorToTag<AccumulatorTags>::exec((BaseType &)*this, resolveAlias(tag), v);
        vigra_precondition(found, std::string("PythonAccumulatorArray::get(): Tag '") + tag + "' not found.");
        return v.result;
    }
    
    void merge(PythonAccumulatorArray const & o)
    {
        BaseType::merge(o);
    }
};

template <class Accumulators, unsigned int ndim, class T>
PythonAccumulatorArray<typename CoupledIteratorType<ndim, 
                                      typename AccumulatorValueTypeTraits<T>::type, npy_uint32>::type::value_type, Accumulators> *
pythonRegionInspect(NumpyArray<ndim, T> in, 
                    NumpyArray<ndim, Singleband<npy_uint32> > labels,
                    python::object tags)
{
    typedef typename CoupledIteratorType<ndim, typename AccumulatorValueTypeTraits<T>::type, npy_uint32>::type Iterator;
    typedef Iterator::value_type Handle;
    typedef PythonAccumulatorArray<Handle, Accumulators> Accu;
    
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
        
        Iterator i     = createCoupledIterator(in, labels),
                 end   = i.getEndIterator();
        collectStatistics(i, end, *res);
    }
    
    return res.release();
}

} // namespace acc1

template <class T, class Accumulators>
void definePythonAccumulator()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename acc1::AccumulatorValueTypeTraits<T>::type ResultType;
    
    typedef acc1::PythonAccumulator<ResultType, Accumulators> Accu;
    class_<Accu>("Accumulator", python::no_init)
        .def("__getitem__", &Accu::get)
        .def("activeNames", &Accu::activeNames)
        .def("names", &Accu::names)
        .def("merge", &Accu::merge)
        ;
    
    def("extractFeatures", &acc1::pythonInspect<Accumulators, 2, T>,
          (arg("image"), arg("tags") = ""),
          return_value_policy<manage_new_object>());
    
    def("extractFeatures", &acc1::pythonInspect<Accumulators, 3, T>,
          (arg("volume"), arg("tags") = ""),
          return_value_policy<manage_new_object>());
}

template <class T, class Accumulators>
void definePythonAccumulatorArray()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename acc1::AccumulatorValueTypeTraits<T>::type ResultType;
    typedef typename CoupledIteratorType<2, ResultType, npy_uint32>::type Iterator;
    typedef Iterator::value_type Handle;
    
    typedef acc1::PythonAccumulatorArray<Handle, Accumulators> Accu;
    class_<Accu>("AccumulatorArray2", python::no_init)
        .def("__getitem__", &Accu::get)
        .def("activeNames", &Accu::activeNames)
        .def("names", &Accu::names)
        .def("merge", &Accu::merge)
        ;
    
    def("extractRegionFeatures", &acc1::pythonRegionInspect<Accumulators, 2, T>,
          (arg("image"), arg("labels"), arg("tags") = ""),
          return_value_policy<manage_new_object>());
    
    // def("extractFeatures", &acc1::pythonInspect<Accumulators, 3, T>,
          // (arg("volume"), arg("tags") = ""),
          // return_value_policy<manage_new_object>());
}

void defineAccumulators()
{
    using namespace python;
    using namespace vigra::acc1;

    docstring_options doc_options(true, true, false);
    
    NumpyArrayConverter<NumpyArray<1, float> >();
    NumpyArrayConverter<NumpyArray<2, MultiArrayIndex> >();
    NumpyArrayConverter<NumpyArray<3, double> >();

    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, 
                   UnbiasedVariance, UnbiasedSkewness, UnbiasedKurtosis,
                   Minimum, Maximum, StandardQuantiles<AutoRangeHistogram<100> > 
                   > ScalarAccumulators;
    definePythonAccumulator<Singleband<float>, ScalarAccumulators>();
    
    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, 
                   Covariance, Principal<Variance>, Principal<Skewness>, Principal<Kurtosis>,
                   Minimum, Maximum, Principal<Minimum>, Principal<Maximum>
                   > VectorAccumulators;
    definePythonAccumulator<TinyVector<float, 2>, VectorAccumulators>();
    definePythonAccumulator<TinyVector<float, 3>, VectorAccumulators>();
    definePythonAccumulator<TinyVector<float, 4>, VectorAccumulators>();

    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, 
                   Minimum, Maximum, StandardQuantiles<GlobalRangeHistogram<100> >,
                   GeometricCenter, PrincipalRadii, PrincipalCoordSystem,
                   CenterOfMass, MomentsOfInertia, CoordSystemOfInertia,
                   Select<Coord<Minimum>, Coord<Maximum>, Coord<ArgMinWeight>, Coord<ArgMaxWeight> >,
                   DataArg<1>, WeightArg<1>, LabelArg<2>
                   > ScalarRegionAccumulators;
    definePythonAccumulatorArray<Singleband<float>, ScalarRegionAccumulators>();
}

// TODO:
//  * nested Select
//  * Multiband support
//  * implement PythonAccumulatorArray::merge()
//  * implement label remapping in merge()
//  * is there a good implementation of merge for histogramms with different mapping?
//  * multiband histograms
//  * ensure that accumulators promote float arguments to double
//  * general refactoring
//  * better names for PrincipalRadii, PrincipalCoordSystem, MomentsOfInertia, CoordSystemOfInertia
//  * tests and docu

} // namespace vigra
