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

AliasMap createTagToAlias(ArrayVector<std::string> const & names);

AliasMap createAliasToTag(AliasMap const & tagToAlias);

ArrayVector<std::string> createSortedNames(AliasMap const & tagToAlias);

template <class BaseType, class GetVisitor>
struct PythonAccumulator
: public BaseType
{
    typedef typename BaseType::AccumulatorTags AccumulatorTags;
    
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
        GetVisitor v(permutation_);
        
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
    
    static void definePythonClass(char const * classname)
    {

      //docstring_options doc_options(true, true, false); //me (?)

     // std::string classname_str(classname);
     // if (strcmp("MultibandFeatures2D", classname)==0) {
      
        python::class_<PythonAccumulator>
	  (classname, 
	   "An instance of this accumulator class is returned by\n"
	   " :func:`extractFeatures` with input data type ... . \n"
           "The object contains the computed features \n"
           "(i.e. the selected features and their dependencies).\n"
	   "For details see `Feature Accumulators <../vigra/group__FeatureAccumulators.html>`_ \n"
           "in the vigra C++ documentation. \n\n",
	   python::no_init)
	  .def("__getitem__", &PythonAccumulator::get, 
	       "Returns the value of the 'feature'. The return type is a \n"
	       " float or a numpy array of appropriate shape.\n",
	       python::arg("feature") )
	  .def("isActive", &PythonAccumulator::isActive,
	       "Returns True if 'feature' has been computed and False otherwise.\n",
	       python::arg("feature") )
	  .def("activeNames", &PythonAccumulator::activeNames,
	       "Returns a list of all computed features.\n")
	  .def("names", &PythonAccumulator::names,
	       "Returns a list of all supported features for the given input data array.\n"
	       )
	  .def("merge", &PythonAccumulator::merge,
	       "Merge features with the features from accumulator 'b'.\n"
	       "Merging is not supported for all statistics  (?) , see the\n"
	       "C++ documentation for more information (?) .\n",
	       python::arg("b"))
            ;
     // } else {
//	python::class_<PythonAccumulator>
//	  (classname, 
//	   " (?) An instance of this accumulator class is returned by \n"
//	   "extractFeatures() and contains the computed features. \n"
//	   "See MultibandFeatures2D for documentation of member functions.\n",
//	   python::no_init)
//	  .def("__getitem__", &PythonAccumulator::get)
//	  .def("isActive", &PythonAccumulator::isActive)
//	  .def("activeNames", &PythonAccumulator::activeNames)
//	  .def("names", &PythonAccumulator::names)
//	  .def("merge", &PythonAccumulator::merge)
//            ;
  //    }
    }
    
    static void definePythonArrayClass(char const * classname)
    {
        python::class_<PythonAccumulator>
	  (classname, 
	   "An instance of this accumulator class is returned by \n"
	   "extractRegionFeatures() and contains the computed features. \n",
	   python::no_init)
	  .def("__getitem__", &PythonAccumulator::get,
	       "Returns the value of the 'feature'. The return type is a \n"
	       " a numpy array of appropriate shape. The first array index \n"
	       "is the region label. \n",
	       python::arg("feature"))
	  .def("isActive", &PythonAccumulator::isActive, 
	       "Returns True if 'feature' is selected (i.e. \n"
	       "'feature' has been computed) and False otherwise.\n",
	       python::arg("feature"))
	  .def("activeNames", &PythonAccumulator::activeNames,
	       "Returns a list of all selected features.\n")
	  .def("names", &PythonAccumulator::names,
	       "Returns a list of all supported features for the given input.\n")
	  .def("merge", &PythonAccumulator::merge,
	       "Merge statistics with the statistics from accumulator 'b'.\n"
	       "Merging is not supported for all statistics, see the\n"
	       "C++ documentation for more information (?) .\n",
	       python::arg("b"))
	  .def("merge", &PythonAccumulator::remappingMerge,
	       "Merge statistics with the statistics from accumulator 'b'\n"
	       "using label remapping. The label map 'labelmap' maps the \n"
	       "labels of accumulator 'b' (i.e. the index is the old label \n"
	       "and the value is the new label). 'labelmap' must be a \n"
	       "numpy array with dtype=numpy.uint32.\n",
	       (python::arg("b"), python::arg("labelmap")))
	  .def("merge", &PythonAccumulator::mergeRegions,
	       "Merge statistics from region 'k' into region 'j'.\n",
	       (python::arg("j"), python::arg("k")))
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
        static const AliasMap a = createTagToAlias(PythonAccumulator::tagNames());
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
pythonInspectWithHistogram(NumpyArray<ndim, Singleband<T> > in, python::object tags,
                           python::object histogramRange, int binCount)
{
    std::auto_ptr<Accumulator> res(new Accumulator);
    if(pythonActivateTags(*res, tags))
    {
        pythonHistogramOptions(*res, histogramRange, binCount);
        
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
                    python::object tags,
                    python::object ignore_label)
{
    typedef typename CoupledIteratorType<ndim, T, npy_uint32>::type Iterator;
    
    TinyVector<npy_intp, ndim> permutation = in.template permuteLikewise<ndim>();
    
    std::auto_ptr<Accumulator> res(new Accumulator(permutation));
    if(pythonActivateTags(*res, tags))
    {
        if(ignore_label != python::object())
            res->ignoreLabel(python::extract<MultiArrayIndex>(ignore_label)());
            
        PyAllowThreads _pythread;
        
        Iterator i     = createCoupledIterator(in, labels),
                 end   = i.getEndIterator();
        collectStatistics(i, end, *res);
    }
    
    return res.release();
}

template <class Accumulator, unsigned int ndim, class T>
Accumulator *
pythonRegionInspectWithHistogram(NumpyArray<ndim, Singleband<T> > in, 
                    NumpyArray<ndim, Singleband<npy_uint32> > labels,
                    python::object tags, python::object histogramRange, int binCount,
                    python::object ignore_label)
{
    typedef typename CoupledIteratorType<ndim, T, npy_uint32>::type Iterator;
    
    TinyVector<npy_intp, ndim> permutation = in.template permuteLikewise<ndim>();
    
    std::auto_ptr<Accumulator> res(new Accumulator(permutation));
    if(pythonActivateTags(*res, tags))
    {
        pythonHistogramOptions(*res, histogramRange, binCount);
        if(ignore_label != python::object())
            res->ignoreLabel(python::extract<MultiArrayIndex>(ignore_label)());
                    
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
                             python::object tags,
                             python::object ignore_label)
{
    typedef typename CoupledIteratorType<ndim, Multiband<T>, npy_uint32>::type Iterator;
    
    TinyVector<npy_intp, ndim-1> permutation = in.template permuteLikewise<ndim-1>();
    
    std::auto_ptr<Accumulator> res(new Accumulator(permutation));
    if(pythonActivateTags(*res, tags))
    {
        if(ignore_label != python::object())
            res->ignoreLabel(python::extract<MultiArrayIndex>(ignore_label)());
            
        PyAllowThreads _pythread;
        
        Iterator i     = createCoupledIterator(MultiArrayView<ndim, Multiband<T>, StridedArrayTag>(in), labels),
                 end   = i.getEndIterator();
        collectStatistics(i, end, *res);
    }
    
    return res.release();
}

} // namespace acc1

template <class T, class Accumulators>
void definePythonAccumulatorSingleband(char const * classname)
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChain<T, Accumulators>, acc1::GetTag_Visitor> Accu;
    
    Accu::definePythonClass(classname);

    std::string test_string = "Docu_test\n\n";
    
    def("extractFeatures", &acc1::pythonInspectWithHistogram<Accu, 2, T>,
          (arg("image"), arg("features") = "all", 
           arg("histogramRange") = "globalminmax", arg("binCount") = 64),
	return_value_policy<manage_new_object>(),
	"");
    
    def("extractFeatures", &acc1::pythonInspectWithHistogram<Accu, 3, T>,
          (arg("volume"), arg("features") = "all", 
           arg("histogramRange") = "globalminmax", arg("binCount") = 64),
          return_value_policy<manage_new_object>());
}

template <class T, class Accumulators>
void definePythonAccumulator(char const * classname)
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChain<T, Accumulators>, acc1::GetTag_Visitor> Accu;

    Accu::definePythonClass(classname);
    
    def("extractFeatures", &acc1::pythonInspect<Accu, 2, T>,
          (arg("image"), arg("features") = "all"),
          return_value_policy<manage_new_object>());
    
    def("extractFeatures", &acc1::pythonInspect<Accu, 3, T>,
          (arg("volume"), arg("features") = "all"),
          return_value_policy<manage_new_object>());

}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorMultiband(char const * classname)
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, Multiband<T> >::type::value_type ResultType;    
    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChain<ResultType, Accumulators>, acc1::GetTag_Visitor> Accu;
    
    Accu::definePythonClass(classname);
            
    std::string argname = N == 3 
                             ? "image"
                             : "volume";

    //for documentation only
    std::string doc_string;
    if (N==3) {
      doc_string.append
	("Extract global features (e.g. Mean, Variance, Minimum etc.) from the input array. As input array the overloaded extractFeatures() function supports 2D or 3D numpy arrays with arbitrary many channels. For different input arrays a different accumulator class is returned. The element type of the input array must be **dtype=float32**. The set of available features depends on the input array, e.g. the 'Histogram' feature is only supported for singleband arrays. Call :func:`MultibandFeatures2D.names` (and others...) to get a list of all available features for the respective input array. Features are can be:\n\n    + 'feature'    compute single feature (and its dependencies)\n\n    + ['feat1', 'feat2',...]    compute given feature set\n\n    + 'all'    compute all features (default)\n\n    + None or ''    compute nothing (usefull to get list of available features)\n\nTo compute per-region features, use :func:`extractRegionFeatures`.\n\nThis instance is called for input array :class:`vigra.VigraArray` with two spatial axes and two or more than four channels. The function returns an instance of :class:`MultibandFeatures2D`.\n\nFor details see `Feature Accumulators <../vigra/group__FeatureAccumulators.html>`_ in the vigra C++ documentation.\n\n");
    } else {
      doc_string.append("This instance of :func:`extractFeatures` is called for input array :class:`vigra.VigraArray` with three spatial axes and two or more than four channels. An instance of :class:`MultibandFeatures3D` is returned.\n\n");
    }
    
    def("extractFeatures", &acc1::pythonInspectMultiband<Accu, N, T>,
          (arg(argname.c_str()), arg("features") = "all"),
	doc_string.c_str(),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorArraySingleband(char const * classname)
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, T, npy_uint32>::type Iterator;
    typedef typename Iterator::value_type Handle;
    
    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChainArray<Handle, Accumulators>, acc1::GetArrayTag_Visitor> Accu;
    Accu::definePythonArrayClass(classname);
    
    std::string argname = N == 2 
                             ? "image"
                             : "volume";
    
    def("extractRegionFeatures", &acc1::pythonRegionInspectWithHistogram<Accu, N, T>,
          (arg(argname.c_str()), arg("labels"), arg("features") = "all", 
           arg("histogramRange") = "globalminmax", arg("binCount") = 64, arg("ignoreLabel")=python::object()),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorArray(char const * classname)
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, T, npy_uint32>::type Iterator;
    typedef typename Iterator::value_type Handle;
    
    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChainArray<Handle, Accumulators>, acc1::GetArrayTag_Visitor> Accu;
    Accu::definePythonArrayClass(classname);
    
    std::string argname = N == 2 
                             ? "image"
                             : "volume";
    
    def("extractRegionFeatures", &acc1::pythonRegionInspect<Accu, N, T>,
          (arg(argname.c_str()), arg("labels"), arg("features") = "all", arg("ignoreLabel")=python::object()),
          return_value_policy<manage_new_object>());
}

template <unsigned int N, class T, class Accumulators>
void definePythonAccumulatorArrayMultiband(char const * classname)
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    typedef typename CoupledIteratorType<N, Multiband<T>, npy_uint32>::type::value_type Handle;
    typedef acc1::PythonAccumulator<acc1::DynamicAccumulatorChainArray<Handle, Accumulators>, acc1::GetArrayTag_Visitor> Accu;
    Accu::definePythonArrayClass(classname);
        
    std::string argname = N == 3 
                             ? "image"
                             : "volume";
    
    def("extractRegionFeatures", &acc1::pythonRegionInspectMultiband<Accu, N, T>,
          (arg(argname.c_str()), arg("labels"), arg("features") = "all", arg("ignoreLabel")=python::object()),
          return_value_policy<manage_new_object>());
}

} // namespace vigra

#endif // VIGRA_PYTHONACCUMULATOR_HXX
