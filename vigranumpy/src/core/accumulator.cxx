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
    
    python::object to_python(signed char t) const { return python::object(t); }
    python::object to_python(signed short t) const { return python::object(t); }
    python::object to_python(signed int t) const { return python::object(t); }
    python::object to_python(signed long t) const { return python::object(t); }
    python::object to_python(signed long long t) const { return python::object(t); }
    python::object to_python(unsigned char t) const { return python::object(t); }
    python::object to_python(unsigned short t) const { return python::object(t); }
    python::object to_python(unsigned int t) const { return python::object(t); }
    python::object to_python(unsigned long t) const { return python::object(t); }
    python::object to_python(unsigned long long t) const { return python::object(t); }
    python::object to_python(float t) const { return python::object(t); }
    python::object to_python(double t) const { return python::object(t); }
    python::object to_python(long double t) const { return python::object(t); }
    
    template <class T, int N>
    python::object to_python(TinyVector<T, N> const & t) const
    {
        NumpyArray<1, T> a = NumpyArray<1, T>(Shape1(N));
        for(int k=0; k<N; ++k)
            a(k) = t[k];
        return python::object(a);
    }
    
    template <class T, class Stride>
    python::object to_python(MultiArrayView<1, T, Stride> const & t) const
    {
        NumpyArray<1, T> a(t);
        return python::object(a);
    }
    
    template <class T>
    python::object to_python(Matrix<T> const & t) const
    {
        return python::object(t);
    }
    
    template <class T1, class T2>
    python::object to_python(std::pair<T1, T2> const & t) const
    {
        return python::make_tuple(to_python(t.first), to_python(t.second));
    }
    
    template <class TAG, class Accu>
    void exec(Accu & a) const
    {
        result = to_python(get<TAG>(a));
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
    
    template <class TAG, class T, class Alloc, class Accu>
    struct ToPythonArray<TAG, MultiArray<1, T, Alloc>, Accu>
    {
        static python::object exec(Accu & a)
        {
            unsigned int n = a.regionCount();
            MultiArrayIndex N = get<TAG>(a, 0).shape(0);
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
    
    template <class TAG, class T1, class T2, class Accu>
    struct ToPythonArray<TAG, std::pair<T1, T2>, Accu>
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
        this->result = to_python(get<Global<TAG> >(a));
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

template <class Accu>
void pythonActivateTags(Accu & a, python::object tags)
{
    if(python::len(tags) == 0)
    {
        a.activateAll();
    }
    else if(PyString_Check(tags.ptr()))
    {
        a.activate(python::extract<std::string>(tags)());
    }
    else
    {
        for(int k=0; k<python::len(tags); ++k)
        {
            a.activate(python::extract<std::string>(tags[k])());
        }
    }
}

template <class Accumulators, unsigned int ndim, class T>
PythonAccumulator<typename AccumulatorValueTypeTraits<T>::type, Accumulators> *
pythonInspect(NumpyArray<ndim, T> in, python::object tags)
{
    typedef PythonAccumulator<typename AccumulatorValueTypeTraits<T>::type, Accumulators> Accu;
    
    std::auto_ptr<Accu> res(new Accu);
    pythonActivateTags(*res, tags);
    
    {
        PyAllowThreads _pythread;
        
        collectStatistics(in.begin(), in.end(), *res);
    }
    
    return res.release();
}

template <class Accumulators, unsigned int ndim, class T>
PythonAccumulator<typename CoupledIteratorType<ndim, Multiband<T> >::type::value_type, Accumulators> *
pythonInspectMultiband(NumpyArray<ndim, Multiband<T> > in, python::object tags)
{
    typedef typename CoupledIteratorType<ndim, Multiband<T> >::type Iterator;
    typedef Iterator::value_type Handle;
    typedef PythonAccumulator<Handle, Accumulators> Accu;
    
    std::auto_ptr<Accu> res(new Accu);
    pythonActivateTags(*res, tags);
    
    {
        PyAllowThreads _pythread;
        
        Iterator i   = createCoupledIterator(MultiArrayView<ndim, Multiband<T>, StridedArrayTag>(in)),
                 end = i.getEndIterator();
        collectStatistics(i, end, *res);
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
    
    void remappingMerge(PythonAccumulatorArray const & o, NumpyArray<1, npy_uint32> labelMapping)
    {
        BaseType::merge(o, labelMapping);
    }
    
    void mergeRegions(npy_uint32 i, npy_uint32 j)
    {
        BaseType::merge(i, j);
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
    pythonActivateTags(*res, tags);
    
    {
        PyAllowThreads _pythread;
        
        Iterator i     = createCoupledIterator(in, labels),
                 end   = i.getEndIterator();
        collectStatistics(i, end, *res);
    }
    
    return res.release();
}

template <class Accumulators, unsigned int ndim, class T>
PythonAccumulatorArray<typename CoupledIteratorType<ndim, Multiband<T>, npy_uint32>::type::value_type, Accumulators> *
pythonRegionInspectMultiband(NumpyArray<ndim, Multiband<T> > in, 
                             NumpyArray<ndim-1, Singleband<npy_uint32> > labels,
                             python::object tags)
{
    typedef typename CoupledIteratorType<ndim, Multiband<T>, npy_uint32>::type Iterator;
    typedef Iterator::value_type Handle;
    typedef PythonAccumulatorArray<Handle, Accumulators> Accu;
    
    std::auto_ptr<Accu> res(new Accu);
    pythonActivateTags(*res, tags);
    
    {
        PyAllowThreads _pythread;
        
        Iterator i     = createCoupledIterator(MultiArrayView<ndim, Multiband<T>, StridedArrayTag>(in), labels),
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

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorMultiband()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef acc1::PythonAccumulator<typename CoupledIteratorType<N, Multiband<T> >::type::value_type, 
                                    Accumulators> Accu;
    class_<Accu>("Accumulator", python::no_init)
        .def("__getitem__", &Accu::get)
        .def("activeNames", &Accu::activeNames)
        .def("names", &Accu::names)
        .def("merge", &Accu::merge)
        ;
        
    std::string argname = N == 3 
                             ? "image"
                             : "volume";
    
    def("extractFeatures", &acc1::pythonInspectMultiband<Accumulators, N, T>,
          (arg(argname.c_str()), arg("tags") = ""),
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
        .def("merge", &Accu::remappingMerge)
        .def("merge", &Accu::mergeRegions)
        ;
    
    def("extractRegionFeatures", &acc1::pythonRegionInspect<Accumulators, 2, T>,
          (arg("image"), arg("labels"), arg("tags") = ""),
          return_value_policy<manage_new_object>());
    
    def("extractRegionFeatures", &acc1::pythonRegionInspect<Accumulators, 3, T>,
          (arg("volume"), arg("labels"), arg("tags") = ""),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorArrayMultiband()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef acc1::PythonAccumulatorArray<typename CoupledIteratorType<N, Multiband<T>, npy_uint32>::type::value_type, 
                                         Accumulators> Accu;
    class_<Accu>("Accumulator", python::no_init)
        .def("__getitem__", &Accu::get)
        .def("activeNames", &Accu::activeNames)
        .def("names", &Accu::names)
        .def("merge", &Accu::merge)
        .def("merge", &Accu::remappingMerge)
        .def("merge", &Accu::mergeRegions)
        ;
        
    std::string argname = N == 3 
                             ? "image"
                             : "volume";
    
    def("extractRegionFeatures", &acc1::pythonRegionInspectMultiband<Accumulators, N, T>,
          (arg(argname.c_str()), arg("labels"), arg("tags") = ""),
          return_value_policy<manage_new_object>());
}

void defineGlobalAccumulators()
{
    using namespace python;
    using namespace vigra::acc1;

    docstring_options doc_options(true, true, false);
    
    NumpyArrayConverter<NumpyArray<1, npy_uint32> >();
    NumpyArrayConverter<NumpyArray<1, float> >();
    NumpyArrayConverter<NumpyArray<2, MultiArrayIndex> >();
    NumpyArrayConverter<NumpyArray<3, float> >();
    NumpyArrayConverter<NumpyArray<3, double> >();
    
    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, Covariance, 
                   Principal<Variance>, Principal<Skewness>, Principal<Kurtosis>,
                   Principal<CoordinateSystem>,
                   Minimum, Maximum, Principal<Minimum>, Principal<Maximum>
                   > VectorAccumulators;

    definePythonAccumulatorMultiband<3, float, VectorAccumulators>();
    definePythonAccumulatorMultiband<4, float, VectorAccumulators>();
    
    definePythonAccumulator<TinyVector<float, 2>, VectorAccumulators>();
    definePythonAccumulator<TinyVector<float, 3>, VectorAccumulators>();
    definePythonAccumulator<TinyVector<float, 4>, VectorAccumulators>();

    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, 
                   UnbiasedVariance, UnbiasedSkewness, UnbiasedKurtosis,
                   Minimum, Maximum, StandardQuantiles<AutoRangeHistogram<100> > 
                   > ScalarAccumulators;
    definePythonAccumulator<Singleband<float>, ScalarAccumulators>();
}

void defineRegionAccumulators()
{
    using namespace python;
    using namespace vigra::acc1;

    docstring_options doc_options(true, true, false);
    
    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, Covariance, 
                   Principal<Variance>, Principal<Skewness>, Principal<Kurtosis>,
                   Principal<CoordinateSystem>,
                   Minimum, Maximum, Principal<Minimum>, Principal<Maximum>,
                   Select<GeometricCenter, PrincipalRadii, PrincipalCoordSystem,
                          Coord<Minimum>, Coord<Maximum>, Principal<Coord<Skewness> >, Principal<Coord<Kurtosis> > >,
                   DataArg<1>, LabelArg<2>
                   > VectorRegionAccumulators;

    definePythonAccumulatorArrayMultiband<3, float, VectorRegionAccumulators>();
    // definePythonAccumulatorArrayMultiband<4, float, VectorRegionAccumulators>();
    
    // definePythonAccumulatorArray<TinyVector<float, 2>, VectorRegionAccumulators>();
    // definePythonAccumulatorArray<TinyVector<float, 3>, VectorRegionAccumulators>();
    // definePythonAccumulatorArray<TinyVector<float, 4>, VectorRegionAccumulators>();

    // typedef Select<Count, Mean, Variance, Skewness, Kurtosis, 
                   // Minimum, Maximum, StandardQuantiles<GlobalRangeHistogram<100> >,
                   // GeometricCenter, PrincipalRadii, PrincipalCoordSystem,
                   // CenterOfMass, MomentsOfInertia, CoordSystemOfInertia,
                   // Select<Coord<Minimum>, Coord<Maximum>, Coord<ArgMinWeight>, Coord<ArgMaxWeight>, 
                          // Principal<Coord<Skewness> >, Principal<Coord<Kurtosis> >, 
                          // Principal<Weighted<Coord<Skewness> > >, Principal<Weighted<Coord<Kurtosis> > > >,
                   // DataArg<1>, WeightArg<1>, LabelArg<2>
                   // > ScalarRegionAccumulators;
    // definePythonAccumulatorArray<Singleband<float>, ScalarRegionAccumulators>();
}

void defineAccumulators()
{
    defineGlobalAccumulators();
    defineRegionAccumulators();
}

// TODO:
//  * nested Select
//  * Multiband support
//  * implement PythonAccumulatorArray::merge()
//  * check that merge skips inactive accumulators
//  * implement label remapping in merge()
//  * is there a good implementation of merge for histogramms with different mapping?
//  * multiband histograms
//  * ensure that accumulators promote float arguments to double
//  * general refactoring
//  * better names for PrincipalRadii, PrincipalCoordSystem, MomentsOfInertia, CoordSystemOfInertia
//  * tests and docu
//  * speed-up compilation of Python bindings

} // namespace vigra
