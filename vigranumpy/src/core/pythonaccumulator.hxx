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

#ifndef VIGRA_PYTHONACCUMULATOR_HXX
#define VIGRA_PYTHONACCUMULATOR_HXX

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

namespace acc 
{

struct GetTag_Visitor
{
    mutable python::object result;
    
    GetTag_Visitor()
    {}
    
    template <class Permutation>
    GetTag_Visitor(Permutation const & p)
    {}

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
        template <class Permutation>
        static python::object exec(Accu & a, Permutation const &)
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
        template <class Permutation>
        static python::object exec(Accu & a, Permutation const & p)
        {
            unsigned int n = a.regionCount();
            Shape2 s(n, N);
            NumpyArray<2, T> res(s);
            
            for(unsigned int k=0; k<n; ++k)
                for(int j=0; j<N; ++j)
                    res(k, p(j)) = get<TAG>(a, k)[j];
            return python::object(res);
        }
    };
    
    template <class TAG, class T, class Alloc, class Accu>
    struct ToPythonArray<TAG, MultiArray<1, T, Alloc>, Accu>
    {
        template <class Permutation>
        static python::object exec(Accu & a, Permutation const & p)
        {
            unsigned int n = a.regionCount();
            MultiArrayIndex N = get<TAG>(a, 0).shape(0);
            Shape2 s(n, N);
            NumpyArray<2, T> res(s);
            
            for(unsigned int k=0; k<n; ++k)
                for(int j=0; j<N; ++j)
                    res(k, p(j)) = get<TAG>(a, k)[j];
            return python::object(res);
        }
    };
    
    template <class TAG, class T, class Accu>
    struct ToPythonArray<TAG, Matrix<T>, Accu>
    {
        template <class Permutation>
        static python::object exec(Accu & a, Permutation const & p)
        {
            unsigned int n = a.regionCount();
            Shape2 m = get<TAG>(a, 0).shape();
            Shape3 s(n, m[0], m[1]);
            NumpyArray<3, T> res(s);
            
            for(unsigned int k=0; k<n; ++k)
                for(int i=0; i<m[0]; ++i)
                    for(int j=0; j<m[1]; ++j)
                        res(k, p(i), p(j)) = get<TAG>(a, k)(i, j);
            return python::object(res);
        }
    };
    
    template <class TAG, class T, class Accu>
    struct ToPythonArray<TAG, Error__Attempt_to_access_inactive_statistic<T>, Accu>
    {
        template <class Permutation>
        static python::object exec(Accu & a, Permutation const & p)
        {
            vigra_precondition(false, "PythonAccumulator::get(): Attempt to access inactive statistic.");
            return python::object();
        }
    };
    
    template <class TAG, class T1, class T2, class Accu>
    struct ToPythonArray<TAG, std::pair<T1, T2>, Accu>
    {
        template <class Permutation>
        static python::object exec(Accu & a, Permutation const & p)
        {
            vigra_precondition(false, "PythonAccumulator::get(): Export for this statistic is not implemented, sorry.");
            return python::object();
        }
    };
    
    struct CoordPermutation
    {
        ArrayVector<npy_intp> permutation_;
        
        CoordPermutation()
        {}
        
        template <class Permute>
        CoordPermutation(Permute const & p)
        : permutation_(p.begin(), p.end())
        {}
        
        template <class T>
        T operator()(T const & t) const
        {
            return permutation_[t];
        }
    };
    
    struct IdentityPermutation
    {
        template <class T>
        T operator()(T const & t) const
        {
            return t;
        }
    };
    
    CoordPermutation coord_permutation_;
    
    GetArrayTag_Visitor()
    {}
    
    template <class Permute>
    GetArrayTag_Visitor(Permute const & p)
    : coord_permutation_(p)
    {}
    
    template <class TAG, class Accu>
    void exec(Accu & a) const
    {
        exec(a, (TAG *)0);
    }
    
    template <class Accu, class TAG>
    void exec(Accu & a, TAG *) const
    {
        if(IsCoordinateFeature<TAG>::value && !IsPrincipalFeature<TAG>::value)
            this->result = ToPythonArray<TAG, typename LookupTag<TAG, Accu>::value_type, Accu>::exec(a, coord_permutation_);
        else
            this->result = ToPythonArray<TAG, typename LookupTag<TAG, Accu>::value_type, Accu>::exec(a, IdentityPermutation());
    }
    
    template <class Accu, class TAG>
    void exec(Accu & a, Global<TAG> *) const
    {
        vigra_precondition(IsPrincipalFeature<TAG>::value || !IsCoordinateFeature<TAG>::value,
            "PythonAccumulator::get(): Export of global coordinate features unsupported, sorry.");
        this->result = to_python(get<Global<TAG> >(a));
    }
};

typedef std::map<std::string, std::string> AliasMap;

AliasMap defineAliasMap();

AliasMap * createTagToAlias(ArrayVector<std::string> const & names);

AliasMap * createAliasToTag(AliasMap const & tagToAlias);

ArrayVector<std::string> * createSortedNames(AliasMap const & tagToAlias);

struct PythonFeatureAccumulator
{
    virtual void activate(std::string const & tag) { throw std::runtime_error("abstract function called."); }   
    virtual bool isActive(std::string const & tag) const { throw std::runtime_error("abstract function called."); return false; }
    virtual python::list activeNames() const { throw std::runtime_error("abstract function called."); return python::list(); }
    virtual python::list names() const { throw std::runtime_error("abstract function called."); return python::list(); }
    virtual python::object get(std::string const & tag) { throw std::runtime_error("abstract function called."); return python::object(); }
    virtual void merge(PythonFeatureAccumulator const & o) { throw std::runtime_error("abstract function called."); }
    virtual PythonFeatureAccumulator * create() const { throw std::runtime_error("abstract function called."); return 0; }
    
    static void definePythonClass()
    {
        python::class_<PythonFeatureAccumulator>(
               "FeatureAccumulator", 
               "An instance of this accumulator class is returned by"
               " :func:`extractFeatures`. "
               "The object contains the computed features "
               "(i.e. the selected features and their dependencies).\n"
               ,python::no_init)
            .def("__getitem__", &PythonFeatureAccumulator::get, 
               "accumulator[feature] returns the value of the 'feature'. The return type is a"
               " float or a numpy array of appropriate shape.\n",
               python::arg("feature") )
            .def("isActive", &PythonFeatureAccumulator::isActive,
               "Returns True if 'feature' has been computed and False otherwise.\n",
               python::arg("feature"))
            .def("activeFeatures", &PythonFeatureAccumulator::activeNames,
               "Returns a list with the names of all computed features.\n")
            .def("keys", &PythonFeatureAccumulator::activeNames,
               "Returns a list with the names of all computed features.\n")
            .def("supportedFeatures", &PythonFeatureAccumulator::names,
               "Returns a list of all supported features for the given input data array.\n")
            .def("merge", &PythonFeatureAccumulator::merge,
               "Merge features with the features from accumulator 'other'. Raises a "
               "TypeError when 'other' is incompatible with 'self'.\n",
               python::arg("other"))
            .def("createAccumulator", &PythonFeatureAccumulator::create,
               "Create an empty accumulator with the same active features as 'self'. "
               "This is useful for merging.\n", 
               python::return_value_policy<python::manage_new_object>())
            ;
    }
};

struct PythonRegionFeatureAccumulator
: public PythonFeatureAccumulator
{
    virtual MultiArrayIndex maxRegionLabel() { throw std::runtime_error("abstract function called."); }
    virtual void mergeAll(PythonRegionFeatureAccumulator const & o) { throw std::runtime_error("abstract function called."); }
    virtual void remappingMerge(PythonFeatureAccumulator const & o, NumpyArray<1, npy_uint32> labelMapping) { throw std::runtime_error("abstract function called."); }
    virtual void mergeRegions(npy_uint32 i, npy_uint32 j) { throw std::runtime_error("abstract function called."); }
    virtual PythonRegionFeatureAccumulator * create() const { throw std::runtime_error("abstract function called."); return 0; }
    
    static void definePythonClass()
    {
        python::class_<PythonRegionFeatureAccumulator>(
                "RegionFeatureAccumulator", 
                "An instance of this accumulator class is returned by "
                ":func:`extractRegionFeatures()` and contains the computed"
                " global and per-region features. \n",
                python::no_init)
            .def("__getitem__", &PythonRegionFeatureAccumulator::get,
               "accumulator[feature] returns the value of the 'feature'. "
               "The return type is a numpy array of appropriate shape. "
               "The first index of the returned arrays is the region label.\n",
               python::arg("feature"))
            .def("maxRegionLabel", &PythonRegionFeatureAccumulator::maxRegionLabel,
                "Return the highest region label in this accumulator.\n")
            .def("isActive", &PythonRegionFeatureAccumulator::isActive, 
               "Returns True if 'feature' has been computed"
               " and False otherwise.\n",
               python::arg("feature"))
            .def("activeFeatures", &PythonRegionFeatureAccumulator::activeNames,
               "Returns a list with the names of all selected features.\n")
            .def("keys", &PythonRegionFeatureAccumulator::activeNames,
               "Returns a list with the names of all selected features.\n")
            .def("supportedFeatures", &PythonRegionFeatureAccumulator::names,
               "Returns a list with the names of all supported features for the given input arrays.\n")
            .def("merge", &PythonRegionFeatureAccumulator::mergeAll,
               "Merge features with the features from accumulator 'other'. "
               "'self' and 'other' must have the same `maxRegionLabel`(), or "
               "'self' must be an empty accumulator (as returned by `create`).\n",
               python::arg("other"))
            .def("merge", &PythonRegionFeatureAccumulator::remappingMerge,
               "Merge features with the features from accumulator 'other'. "
               "The 'labelMap' determines the correspondence of regions between "
               "'self' and 'other' (i.e. region k of accumulator 'other' is "
               "merged into region labelMap[k] of 'self').\n",
               (python::arg("other"), python::arg("labelMap")))
            .def("merge", &PythonRegionFeatureAccumulator::mergeRegions,
               "Merge features from region 'j' into region 'i' of this accumulator.\n",
               (python::arg("i"), python::arg("j")))
            .def("createAccumulator", &PythonRegionFeatureAccumulator::create,
               "Create an empty accumulator with the same active features as 'self'. "
               "This is useful for merging.\n", 
               python::return_value_policy<python::manage_new_object>())
            ;
    }  
};

template <class BaseType, class PythonBaseType, class GetVisitor>
struct PythonAccumulator
: public BaseType, public PythonBaseType
{
    typedef typename BaseType::AccumulatorTags AccumulatorTags;
    typedef PythonBaseType PythonBase;
    
    ArrayVector<npy_intp> permutation_;
    
    PythonAccumulator()
    {}
    
    template <class Permutation>
    PythonAccumulator(Permutation const & p)
    : permutation_(p.begin(), p.end())
    {}
    
    void activate(std::string const & tag)
    {
        vigra_precondition(this->activateImpl(resolveAlias(tag)), 
                            "FeatureAccumulator::activate(): Tag '" + tag + "' not found.");
    }
    
    bool isActive(std::string const & tag) const
    {
        acc_detail::TagIsActive_Visitor v;
        vigra_precondition(this->isActiveImpl(resolveAlias(tag), v), 
                           "FeatureAccumulator::isActive(): Tag '" + tag + "' not found.");
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
        GetVisitor v(permutation_);
        
        vigra_precondition(isActive(tag), "FeatureAccumulator::get(): Tag '" + tag + "' is not active.");
        acc_detail::ApplyVisitorToTag<AccumulatorTags>::exec((BaseType &)*this, resolveAlias(tag), v);
        return v.result;
    }
    
    void merge(PythonFeatureAccumulator const & o)
    {
        PythonAccumulator const * p = dynamic_cast<PythonAccumulator const *>(&o);
        if(p == 0)
        {
            PyErr_SetString(PyExc_TypeError, "FeatureAccumulator::merge(): accumulators are incompatible.");
            python::throw_error_already_set();
        }
        BaseType::merge(*p);
    }
    
    void mergeAll(PythonRegionFeatureAccumulator const & o)
    {
        merge(o);
    }
    
    void remappingMerge(PythonFeatureAccumulator const & o, NumpyArray<1, npy_uint32> labelMapping)
    {
        PythonAccumulator const * p = dynamic_cast<PythonAccumulator const *>(&o);
        if(p == 0)
        {
            PyErr_SetString(PyExc_TypeError, "FeatureAccumulator::merge(): accumulators are incompatible.");
            python::throw_error_already_set();
        }
        BaseType::merge(*p, labelMapping);
    }
    
    void mergeRegions(npy_uint32 i, npy_uint32 j)
    {
        BaseType::merge(i, j);
    }
    
    PythonAccumulator * create() const
    {
        VIGRA_UNIQUE_PTR<PythonAccumulator> a(new PythonAccumulator(permutation_));
        pythonActivateTags(*a, activeNames());
        return a.release();
    }
    
    MultiArrayIndex maxRegionLabel() 
    {
        return BaseType::maxRegionLabel();
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
        static AliasMap * a = VIGRA_SAFE_STATIC(a, createTagToAlias(PythonAccumulator::tagNames()));
        return *a;   
    }
    
    static AliasMap const & aliasToTag()
    {
        static AliasMap * a = VIGRA_SAFE_STATIC(a, createAliasToTag(tagToAlias()));
        return *a;   
    }
    
    static ArrayVector<std::string> const & nameList()
    {
        static ArrayVector<std::string> * n = VIGRA_SAFE_STATIC(n, createSortedNames(tagToAlias()));
        return *n;
    }
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

template <class Accu>
void pythonHistogramOptions(Accu & a, python::object minmax, int binCount)
{
    HistogramOptions options;
    options.setBinCount(binCount);
    
    if(PyString_Check(minmax.ptr()))
    {
        std::string spec = normalizeString(python::extract<std::string>(minmax)());
        if(spec == "globalminmax")
            options.globalAutoInit();
        else if(spec == "regionminmax")
            options.regionAutoInit();
        else
            vigra_precondition(false, 
                "extractFeatures(): invalid histogramRange.");
    }
    else if(python::len(minmax) == 2)
    {
        options.setMinMax(python::extract<double>(minmax[0])(), python::extract<double>(minmax[1])());
    }
    else
        vigra_precondition(false, "extractFeatures(): invalid histogramRange.");
    a.setHistogramOptions(options);
}

template <class Accumulator, unsigned int ndim, class T>
typename Accumulator::PythonBase *
pythonInspect(NumpyArray<ndim, T> in, python::object tags)
{
    VIGRA_UNIQUE_PTR<Accumulator> res(new Accumulator);
    if(pythonActivateTags(*res, tags))
    {
        PyAllowThreads _pythread;
        
        extractFeatures(in.begin(), in.end(), *res);
    }
    
    return res.release();
}

template <class Accumulator, unsigned int ndim, class T>
typename Accumulator::PythonBase *
pythonInspectWithHistogram(NumpyArray<ndim, Singleband<T> > in, python::object tags,
                           python::object histogramRange, int binCount)
{
    VIGRA_UNIQUE_PTR<Accumulator> res(new Accumulator);
    if(pythonActivateTags(*res, tags))
    {
        pythonHistogramOptions(*res, histogramRange, binCount);
        
        PyAllowThreads _pythread;
        
        extractFeatures(in.begin(), in.end(), *res);
    }
    
    return res.release();
}

template <class Accumulator, unsigned int ndim, class T>
typename Accumulator::PythonBase *
pythonInspectMultiband(NumpyArray<ndim, Multiband<T> > in, python::object tags)
{
    typedef typename CoupledIteratorType<ndim, Multiband<T> >::type Iterator;
    
    VIGRA_UNIQUE_PTR<Accumulator> res(new Accumulator);
    if(pythonActivateTags(*res, tags))
    {
        PyAllowThreads _pythread;
        
        Iterator i   = createCoupledIterator(MultiArrayView<ndim, Multiband<T>, StridedArrayTag>(in)),
                 end = i.getEndIterator();
        extractFeatures(i, end, *res);
    }
    
    return res.release();
}

template <class Accumulator, unsigned int ndim, class T>
typename Accumulator::PythonBase *
pythonRegionInspect(NumpyArray<ndim, T> in, 
                    NumpyArray<ndim, Singleband<npy_uint32> > labels,
                    python::object tags,
                    python::object ignore_label)
{
    typedef typename CoupledIteratorType<ndim, T, npy_uint32>::type Iterator;
    
    TinyVector<npy_intp, ndim> permutation = in.template permuteLikewise<ndim>();
    
    VIGRA_UNIQUE_PTR<Accumulator> res(new Accumulator(permutation));
    if(pythonActivateTags(*res, tags))
    {
        if(ignore_label != python::object())
            res->ignoreLabel(python::extract<MultiArrayIndex>(ignore_label)());
            
        PyAllowThreads _pythread;
        
        Iterator i     = createCoupledIterator(in, labels),
                 end   = i.getEndIterator();
        extractFeatures(i, end, *res);
    }
    
    return res.release();
}

template <class Accumulator, unsigned int ndim, class T>
typename Accumulator::PythonBase *
pythonRegionInspectWithHistogram(NumpyArray<ndim, Singleband<T> > in, 
                    NumpyArray<ndim, Singleband<npy_uint32> > labels,
                    python::object tags, python::object histogramRange, int binCount,
                    python::object ignore_label)
{
    typedef typename CoupledIteratorType<ndim, T, npy_uint32>::type Iterator;
    
    TinyVector<npy_intp, ndim> permutation = in.template permuteLikewise<ndim>();
    
    VIGRA_UNIQUE_PTR<Accumulator> res(new Accumulator(permutation));
    if(pythonActivateTags(*res, tags))
    {
        pythonHistogramOptions(*res, histogramRange, binCount);
        if(ignore_label != python::object())
            res->ignoreLabel(python::extract<MultiArrayIndex>(ignore_label)());
                    
        PyAllowThreads _pythread;
        
        Iterator i     = createCoupledIterator(in, labels),
                 end   = i.getEndIterator();
        extractFeatures(i, end, *res);
    }
    
    return res.release();
}

template <class Accumulator, unsigned int ndim, class T>
typename Accumulator::PythonBase *
pythonRegionInspectMultiband(NumpyArray<ndim, Multiband<T> > in, 
                             NumpyArray<ndim-1, Singleband<npy_uint32> > labels,
                             python::object tags,
                             python::object ignore_label)
{
    typedef typename CoupledIteratorType<ndim, Multiband<T>, npy_uint32>::type Iterator;
    
    TinyVector<npy_intp, ndim-1> permutation = in.template permuteLikewise<ndim-1>();
    
    VIGRA_UNIQUE_PTR<Accumulator> res(new Accumulator(permutation));
    if(pythonActivateTags(*res, tags))
    {
        if(ignore_label != python::object())
            res->ignoreLabel(python::extract<MultiArrayIndex>(ignore_label)());
            
        PyAllowThreads _pythread;
        
        Iterator i     = createCoupledIterator(MultiArrayView<ndim, Multiband<T>, StridedArrayTag>(in), labels),
                 end   = i.getEndIterator();
        extractFeatures(i, end, *res);
    }
    
    return res.release();
}

} // namespace acc

template <class T, class Accumulators>
void definePythonAccumulatorSingleband()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef acc::PythonAccumulator<acc::DynamicAccumulatorChain<T, Accumulators>, 
                                    acc::PythonFeatureAccumulator, acc::GetTag_Visitor> Accu;
    
    def("extractFeatures", &acc::pythonInspectWithHistogram<Accu, 2, T>,
          (arg("image"), arg("features") = "all", 
           arg("histogramRange") = "globalminmax", arg("binCount") = 64),
          "\nThis overload of extractFeatures() computes global statistics for a\n"
          "2D scalar input array, e.g. :class:`vigra.ScalarImage`\n\n"
          "Features 'Histogram' and 'Quantiles' are supported for this input.\nOptions are:\n\n"
          "    - histogramRange: lower and upper bound of the histogram\n\n"
          "        + 'globalminmax':  compute and use global minimum/maximum (default)\n"
          "        + [lower, upper]:  provide explicit bounds (float numbers),\n"
          "                           useful to ensure that merge will be allowed.\n\n"
          "    - binCount: number of bins (default: 64).\n\n"
          "Histogram options are ignored when the histogram feature is not selected.\n"
          "Quantiles (0%, 10%, 25%, 50%, 75%, 90%, 100%) are computed from\n"
          "the specified histogram.\n\n",
          return_value_policy<manage_new_object>());
    
    def("extractFeatures", &acc::pythonInspectWithHistogram<Accu, 3, T>,
          (arg("volume"), arg("features") = "all", 
           arg("histogramRange") = "globalminmax", arg("binCount") = 64),
          "Likewise for a scalar 3D input array, e.g. :class:`vigra.ScalarVolume`.\n\n",
          return_value_policy<manage_new_object>());
}

template <class T, class Accumulators>
void definePythonAccumulator()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef acc::PythonAccumulator<acc::DynamicAccumulatorChain<T, Accumulators>, 
                                    acc::PythonFeatureAccumulator, acc::GetTag_Visitor> Accu;
    
    def("extractFeatures", &acc::pythonInspect<Accu, 2, T>,
          (arg("image"), arg("features") = "all"),
          "Likewise for 2D arrays with 3 channels.\n"
          "Histograms and quantiles are not supported for this input.\n\n",
          return_value_policy<manage_new_object>());
    
    def("extractFeatures", &acc::pythonInspect<Accu, 3, T>,
          (arg("volume"), arg("features") = "all"),
          "Likewise for 3D arrays with 3 channels.\n"
          "Histograms and quantiles are not supported for this input.\n\n",
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorMultiband()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, Multiband<T> >::type::value_type ResultType;    
    typedef acc::PythonAccumulator<acc::DynamicAccumulatorChain<ResultType, Accumulators>, 
                                    acc::PythonFeatureAccumulator, acc::GetTag_Visitor> Accu;
    
    std::string argname = N == 3 
                             ? "image"
                             : "volume";
    
    std::string doc_string;
    if (N==3) {
      doc_string +=
        "Extract global features (e.g. Mean, Variance, Minimum, etc.)\n"
        "from the input array ('image' or 'volume'). An accumulator object\n"
        "of type :class:`FeatureAccumulator` is returned that holds the computed\nfeatures.\n\n" 
        "The overloaded function extractFeatures() supports 2D or 3D\n"
        "arrays with arbitrary many channels. The element type of the\n"
        "input array must be **dtype=numpy.float32**. The set of available features\n"
        "depends on the input array. The 'Histogram' feature, for example,\n"
        "is only supported for singleband arrays. Call :func:`supportedFeatures`\n"
        "with the same input array to get a list of all available features\n"
        "for this input.\n\n"
        "The argument 'features' can take the following values:\n\n"
        "   - 'all': compute all supported features (default)\n\n"
        "   - name:  compute a single feature (and its dependencies)\n\n"
        "   - [name1, name2,...]:  compute the given features plus dependencies\n\n"
        "   - None or '':  return an empty accumulator, whose method \n"
        "                  :meth:`~.FeatureAccumulator.supportedFeatures`\n"
        "                  tells you the list of supported features for the\n"
        "                  given input array.\n\n"
        "To compute per-region features, use :func:`extractRegionFeatures`.\n\n"
        "This overload is called for 2D input arrays two or more than\n"
        "four channels. Histograms and quantiles are not supported for\n"
        "this input.\n\n"
        "For further details about the meaning of the features, see\n"
        "`Feature Accumulators <../vigra/group__FeatureAccumulators.html>`_ in the vigra C++ documentation.\n\n";
    } else {
      doc_string +=
        "Overload for 3D arrays with arbitrary many channels.\n"
        "Histograms and quantiles are not supported for this input.\n\n";
    }
    def("extractFeatures", &acc::pythonInspectMultiband<Accu, N, T>,
          (arg(argname.c_str()), arg("features") = "all"),
          doc_string.c_str(),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorArraySingleband()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, T, npy_uint32>::type Iterator;
    typedef typename Iterator::value_type Handle;
    
    typedef acc::PythonAccumulator<acc::DynamicAccumulatorChainArray<Handle, Accumulators>, 
                                    acc::PythonRegionFeatureAccumulator, acc::GetArrayTag_Visitor> Accu;
    
    std::string argname = N == 2 
                             ? "image"
                             : "volume";
    
    std::string doc_string;
    if (N==2) {
      doc_string +=
         "\nThis overload of extractRegionFeatures() computes region statistics\n"
         "for a scalar 2D input array, e.g. :class:`vigra.ScalarImage`.\n\n"
         "Features 'Histogram' and 'Quantiles' are supported for this input. Options are:\n\n"
         "    - histogramRange: lower and upper bound of the histogram\n\n"
         "        + 'globalminmax':  compute and use global minimum/maximum (default)\n"
         "        + 'regionminmax':   use minimum/maximum within each region\n"
         "        + [lower, upper]:  provide explicit bounds (float numbers),\n"
         "                           useful to ensure that merge will be allowed.\n\n"
         "    - binCount: number of bins (default: 64).\n\n"
         "Histogram options are ignored when Histogram feature is not selected.\n"
         "Quantiles (0%, 10%, 25%, 50%, 75%, 90%, 100%) are computed from\n"
         "the specified histogram.\n\n";
    } else {
      doc_string += 
         "Likewise for 3D scalar arrays, e.g. :class:`vigra.ScalarVolume`.\n\n";
    }
    
    def("extractRegionFeatures", &acc::pythonRegionInspectWithHistogram<Accu, N, T>,
          (arg(argname.c_str()), arg("labels"), arg("features") = "all", 
           arg("histogramRange") = "globalminmax", arg("binCount") = 64, arg("ignoreLabel")=python::object()),
          doc_string.c_str(),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorArray()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, T, npy_uint32>::type Iterator;
    typedef typename Iterator::value_type Handle;
    
    typedef acc::PythonAccumulator<acc::DynamicAccumulatorChainArray<Handle, Accumulators>, 
                                    acc::PythonRegionFeatureAccumulator, acc::GetArrayTag_Visitor> Accu;
    
    std::string argname = N == 2 
                             ? "image"
                             : "volume";
    
    std::string doc_string;
    if (N==2) {
      doc_string +=
         "This overload of extractRegionFeatures() is called for\n"
         "2D input arrays with 3 channels.\n\n";
    } else {
      doc_string +=
         "This overload of extractRegionFeatures() is called for\n"
         "3D input arrays with 3 channels.\n\n";
    }
    
    def("extractRegionFeatures", &acc::pythonRegionInspect<Accu, N, T>,
          (arg(argname.c_str()), arg("labels"), arg("features") = "all", arg("ignoreLabel")=python::object()),
          doc_string.c_str(),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorArrayMultiband()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, Multiband<T>, npy_uint32>::type::value_type Handle;
    typedef acc::PythonAccumulator<acc::DynamicAccumulatorChainArray<Handle, Accumulators>, 
                                    acc::PythonRegionFeatureAccumulator, acc::GetArrayTag_Visitor> Accu;
        
    std::string argname = N == 3 
                             ? "image"
                             : "volume";
    
    std::string doc_string;
    if (N==3) {
      doc_string +=
        "\nExtract region features from an input array with **dtype=numpy.float32**\n"
        "and return a :class:`RegionFeatureAccumulator` object.\n\n"
        "Membership of the array elements (pixels) to regions is specified\n"
        "by a 'labels' array with element type **dtype=uint32**.\n\n"
        "The set of available features depends on the input array.\n"
        "Call :func:`supportedRegionFeatures` with the same input and label\n"
        "arrays to get a list of all available features for these inputs.\n\n"
        "The argument 'features' can take the following values:\n\n"
        "   - 'all': compute all supported features (default)\n\n"
        "   - name:  compute a single feature (and its dependencies)\n\n"
        "   - [name1, name2,...]:  compute the given features plus dependencies\n\n"
        "   - None or '':  return an empty accumulator, whose method \n"
        "                  :meth:`~.RegionFeatureAccumulator.supportedFeatures`\n"
        "                  tells you the list of supported features for the\n"
        "                  given input array.\n\n"
        "When the feature name starts with 'Global', the feature is computed\n"
        "globally, i.e. without considering region membership.\n\n"
        "The argument 'ignoreLabel' is useful when the label array contains\n"
        "a background region (usually label 0) that should be ignored during\n"
        "feature computation. If 'ignoreLabel' is None (the default), all\n"
        "region labels are used.\n\n"
        "This overload is called for 2D input arrays with two or more than\n"
        "four channels. Histograms and quantiles are not supported for this\n"
        "input.\n\n"
        "For further details about the meaning of the features, see\n"
        "`Feature Accumulators <../vigra/group__FeatureAccumulators.html>`_ in the vigra C++ documentation.\n\n";

    } else {
      doc_string +=
         "Likewise for a 3D input array  with two or more than four channels.\n"
         "Histograms and quantiles are not supported for this input.\n\n";
    }
    
    def("extractRegionFeatures", &acc::pythonRegionInspectMultiband<Accu, N, T>,
          (arg(argname.c_str()), arg("labels"), arg("features") = "all", arg("ignoreLabel")=python::object()),
          doc_string.c_str(),
          return_value_policy<manage_new_object>());
}

} // namespace vigra

#endif // VIGRA_PYTHONACCUMULATOR_HXX
