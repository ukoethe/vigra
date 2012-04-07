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
    res["Coord<DivideByCount<PowerSum<1> > >"] = "RegionCenter";
    res["Coord<RootDivideByCount<Principal<PowerSum<2> > > >"] = "RegionRadii";
    res["Coord<Principal<CoordinateSystem> >"] = "RegionAxes";
    res["Weighted<Coord<DivideByCount<PowerSum<1> > > >"] = "Weighted<RegionCenter>";
    res["Weighted<Coord<RootDivideByCount<Principal<PowerSum<2> > > > >"] = "Weighted<RegionRadii>";
    res["Weighted<Coord<Principal<CoordinateSystem> > >"] = "Weighted<RegionAxes>";
    return res;
}

AliasMap createTagToAlias(ArrayVector<std::string> const & names)
{
    static const AliasMap aliases = defineAliasMap();
    AliasMap res;
    for(unsigned int k=0; k<names.size(); ++k)
    {
            // treat ScatterMatrixEigensystem as internal 
        if(names[k].find("ScatterMatrixEigensystem") != std::string::npos)
            continue;
        
            // lookup alias names
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

ArrayVector<std::string> createSortedNames(AliasMap const & tagToAlias)
{
    ArrayVector<std::string> res;
    for(AliasMap::const_iterator k = tagToAlias.begin(); k != tagToAlias.end(); ++k)
        res.push_back(k->second);
    std::sort(res.begin(), res.end());
    return res;
}

template <class BaseType, class GetVisitor>
struct PythonAccumulator
: public BaseType
{
    typedef typename BaseType::AccumulatorTags AccumulatorTags;
    
    void activate(std::string const & tag)
    {
        vigra_precondition(this->activateImpl(resolveAlias(tag)), 
                            "PythonAccumulator::activate(): Tag '" + tag + "' not found.");
    }
    
    bool isActive(std::string const & tag) const
    {
        detail::TagIsActive_Visitor v;
        vigra_precondition(isActiveImpl(resolveAlias(tag), v), 
                           "PythonAccumulator::isActive(): Tag '" + tag + "' not found.");
        return v.result;
    }
        
    python::list activeNames() const
    {
        python::list result;
        for(unsigned int k=0; k<nameList().size(); ++k)
            if(isActive(nameList()[k]))
                result.append(python::object(nameList()[k]));
        return result;
    }
    
    python::list names() const
    {
        python::list result;
        for(unsigned int k=0; k<nameList().size(); ++k)
            result.append(python::object(nameList()[k]));
        return result;
    }
    
    python::object get(std::string const & tag)
    {
        GetVisitor v;
        
        vigra_precondition(isActive(tag), "PythonAccumulator::get(): Tag '" + tag + "' is not active.");
        detail::ApplyVisitorToTag<AccumulatorTags>::exec((BaseType &)*this, resolveAlias(tag), v);
        return v.result;
    }
    
    void merge(PythonAccumulator const & o)
    {
        BaseType::merge(o);
    }
    
    void remappingMerge(PythonAccumulator const & o, NumpyArray<1, npy_uint32> labelMapping)
    {
        BaseType::merge(o, labelMapping);
    }
    
    void mergeRegions(npy_uint32 i, npy_uint32 j)
    {
        BaseType::merge(i, j);
    }
    
    static void definePythonClass()
    {
        python::class_<PythonAccumulator>("Accumulator", python::no_init)
            .def("__getitem__", &get)
            .def("isActive", &isActive)
            .def("activeNames", &activeNames)
            .def("names", &names)
            .def("merge", &merge)
            ;
    }
    
    static void definePythonArrayClass()
    {
        python::class_<PythonAccumulator>("Accumulator", python::no_init)
            .def("__getitem__", &get)
            .def("isActive", &isActive)
            .def("activeNames", &activeNames)
            .def("names", &names)
            .def("merge", &merge)
            .def("merge", &remappingMerge)
            .def("merge", &mergeRegions)
            ;
    }
    
  private:
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
    
    static ArrayVector<std::string> const & nameList()
    {
        static const ArrayVector<std::string> n = createSortedNames(tagToAlias());
        return n;
    }
};

template <class T>
struct StripSinglebandTag
{
    typedef T type;
};

template <class T>
struct StripSinglebandTag<Singleband<T> >
{
    typedef T type;
};

template <class Accu>
bool pythonActivateTags(Accu & a, python::object tags)
{
    if(tags == python::object() || python::len(tags) == 0)
        return false;

    if(PyString_Check(tags.ptr()))
    {
        std::string tag = python::extract<std::string>(tags)();
        if(normalizeString(tag) == "all")
            a.activateAll();
        else
            a.activate(tag);
    }
    else
    {
        for(int k=0; k<python::len(tags); ++k)
        {
            a.activate(python::extract<std::string>(tags[k])());
        }
    }
    return true;
}

template <class Accumulator, unsigned int ndim, class T>
Accumulator *
pythonInspect(NumpyArray<ndim, T> in, python::object tags)
{
    std::auto_ptr<Accumulator> res(new Accumulator);
    if(pythonActivateTags(*res, tags))
    {
        PyAllowThreads _pythread;
        
        collectStatistics(in.begin(), in.end(), *res);
    }
    
    return res.release();
}

template <class Accumulator, unsigned int ndim, class T>
Accumulator *
pythonInspectMultiband(NumpyArray<ndim, Multiband<T> > in, python::object tags)
{
    typedef typename CoupledIteratorType<ndim, Multiband<T> >::type Iterator;
    
    std::auto_ptr<Accumulator> res(new Accumulator);
    if(pythonActivateTags(*res, tags))
    {
        PyAllowThreads _pythread;
        
        Iterator i   = createCoupledIterator(MultiArrayView<ndim, Multiband<T>, StridedArrayTag>(in)),
                 end = i.getEndIterator();
        collectStatistics(i, end, *res);
    }
    
    return res.release();
}

template <class Accumulator, unsigned int ndim, class T>
Accumulator *
pythonRegionInspect(NumpyArray<ndim, T> in, 
                    NumpyArray<ndim, Singleband<npy_uint32> > labels,
                    python::object tags)
{
    typedef typename CoupledIteratorType<ndim, typename StripSinglebandTag<T>::type, npy_uint32>::type Iterator;
    
    std::auto_ptr<Accumulator> res(new Accumulator);
    if(pythonActivateTags(*res, tags))
    {
        PyAllowThreads _pythread;
        
        Iterator i     = createCoupledIterator(in, labels),
                 end   = i.getEndIterator();
        collectStatistics(i, end, *res);
    }
    
    return res.release();
}

template <class Accumulator, unsigned int ndim, class T>
Accumulator *
pythonRegionInspectMultiband(NumpyArray<ndim, Multiband<T> > in, 
                             NumpyArray<ndim-1, Singleband<npy_uint32> > labels,
                             python::object tags)
{
    typedef typename CoupledIteratorType<ndim, Multiband<T>, npy_uint32>::type Iterator;
    
    std::auto_ptr<Accumulator> res(new Accumulator);
    if(pythonActivateTags(*res, tags))
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

    typedef typename acc1::StripSinglebandTag<T>::type ResultType;    
    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChain<ResultType, Accumulators>, acc1::GetTag_Visitor> Accu;
    
    Accu::definePythonClass();
        
    def("extractFeatures", &acc1::pythonInspect<Accu, 2, T>,
          (arg("image"), arg("tags") = "all"),
          return_value_policy<manage_new_object>());
    
    def("extractFeatures", &acc1::pythonInspect<Accu, 3, T>,
          (arg("volume"), arg("tags") = "all"),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorMultiband()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, Multiband<T> >::type::value_type ResultType;    
    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChain<ResultType, Accumulators>, acc1::GetTag_Visitor> Accu;
    
    Accu::definePythonClass();
            
    std::string argname = N == 3 
                             ? "image"
                             : "volume";
    
    def("extractFeatures", &acc1::pythonInspectMultiband<Accu, N, T>,
          (arg(argname.c_str()), arg("tags") = "all"),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorArray()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename acc1::StripSinglebandTag<T>::type ResultType;
    typedef typename CoupledIteratorType<N, ResultType, npy_uint32>::type Iterator;
    typedef Iterator::value_type Handle;
    
    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChainArray<Handle, Accumulators>, acc1::GetArrayTag_Visitor> Accu;
    Accu::definePythonArrayClass();
    
    std::string argname = N == 2 
                             ? "image"
                             : "volume";
    
    def("extractRegionFeatures", &acc1::pythonRegionInspect<Accu, N, T>,
          (arg(argname.c_str()), arg("labels"), arg("tags") = "all"),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorArrayMultiband()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, Multiband<T>, npy_uint32>::type::value_type Handle;
    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChainArray<Handle, Accumulators>, acc1::GetArrayTag_Visitor> Accu;
    Accu::definePythonArrayClass();
        
    std::string argname = N == 3 
                             ? "image"
                             : "volume";
    
    def("extractRegionFeatures", &acc1::pythonRegionInspectMultiband<Accu, N, T>,
          (arg(argname.c_str()), arg("labels"), arg("tags") = "all"),
          return_value_policy<manage_new_object>());
}

void defineGlobalAccumulators()
{
    using namespace python;
    using namespace vigra::acc1;

    docstring_options doc_options(true, true, false);
    
    NumpyArrayConverter<NumpyArray<1, npy_uint32> >();
    NumpyArrayConverter<NumpyArray<1, double> >();
    NumpyArrayConverter<NumpyArray<2, MultiArrayIndex> >();
    NumpyArrayConverter<NumpyArray<2, double> >();
    NumpyArrayConverter<NumpyArray<3, double> >();
    
    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, Covariance, 
                   Principal<Variance>, Principal<Skewness>, Principal<Kurtosis>,
                   Principal<CoordinateSystem>,
                   Minimum, Maximum, Principal<Minimum>, Principal<Maximum>
                   > VectorAccumulators;

    definePythonAccumulatorMultiband<3, float, VectorAccumulators>();
    definePythonAccumulatorMultiband<4, float, VectorAccumulators>();
    
    definePythonAccumulator<TinyVector<float, 3>, VectorAccumulators>();

    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, 
                   UnbiasedVariance, UnbiasedSkewness, UnbiasedKurtosis,
                   Minimum, Maximum, StandardQuantiles<AutoRangeHistogram<64> > 
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
                   Select<RegionCenter, RegionRadii, RegionAxes,
                          Coord<Minimum>, Coord<Maximum>, Principal<Coord<Skewness> >, Principal<Coord<Kurtosis> > >,
                   DataArg<1>, LabelArg<2>
                   > VectorRegionAccumulators;

    definePythonAccumulatorArrayMultiband<3, float, VectorRegionAccumulators>();
    definePythonAccumulatorArrayMultiband<4, float, VectorRegionAccumulators>();
    
    definePythonAccumulatorArray<2, TinyVector<float, 3>, VectorRegionAccumulators>();
    definePythonAccumulatorArray<3, TinyVector<float, 3>, VectorRegionAccumulators>();

    typedef Select<Count, Mean, Variance, Skewness, Kurtosis, 
                   Minimum, Maximum, StandardQuantiles<GlobalRangeHistogram<64> >,
                   RegionCenter, RegionRadii, RegionAxes,
                   Weighted<RegionCenter>, Weighted<RegionRadii>, Weighted<RegionAxes>,
                   Select<Coord<Minimum>, Coord<Maximum>, Coord<ArgMinWeight>, Coord<ArgMaxWeight>, 
                          Principal<Coord<Skewness> >, Principal<Coord<Kurtosis> >, 
                          Principal<Weighted<Coord<Skewness> > >, Principal<Weighted<Coord<Kurtosis> > > >,
                   DataArg<1>, WeightArg<1>, LabelArg<2>
                   > ScalarRegionAccumulators;
    definePythonAccumulatorArray<2, Singleband<float>, ScalarRegionAccumulators>();
    definePythonAccumulatorArray<3, Singleband<float>, ScalarRegionAccumulators>();
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
