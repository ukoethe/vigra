#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyanalysis_PyArray_API
#define NO_IMPORT_ARRAY

// all provoke 'warning: "_POSIX_C_SOURCE" redefined' ...
//#include <functional> // unary_function
//#include <cstddef> // size_t 
//#include <algorithm>

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include <boost/unordered_map.hpp>

#include <vigra/multi_pointoperators.hxx>
#include <vigra/coordinate_iterator.hxx>
#include <vigra/object_features.hxx>
#include <vigra/inspectimage.hxx> // A_O_R
#include <vigra/tinyvector.hxx>
#include <vigra/algorithm.hxx>

#include <sstream>

namespace python = boost::python;

namespace vigra
{

python::str py_str(python::object x)
{
    python::str ret(boost::python::handle<>(PyObject_Str(x.ptr())));
    return ret;
}

template <class AccuType, unsigned N>
struct PythonObjectFeatures_extra
    : public type_lists::use_template_list<NestedAccumulators,
                                           StridePairPointer<N, AccuType>,
                                           acc::OctetHistogram
                                          >
{};

using namespace vigra::acc;
using vigra::type_lists::cons;
using vigra::type_lists::nil;

template <class T>
struct PythonObjectFeaturesList
{
    typedef
        cons<CenterOfMass<T>,
        cons<Variance<T>,
        cons<UnbiasedVariance<T>,
        cons<Skewness<T>,
        cons<Kurtosis<T>,
        cons<Sum<T>,
        cons<Min<T>,
        cons<Max<T>,
        cons<CoordVariance<T>,
        cons<CoordUnbiasedVariance<T>,
        cons<CoordSkewness<T>,
        cons<CoordKurtosis<T>,
        cons<CoordSum<T>,
        cons<CoordMin<T>,
        cons<CoordMax<T>,
        cons<WeightedVariance<T>,
        cons<WeightedUnbiasedVariance<T>,
        cons<WeightedSkewness<T>,
        cons<WeightedKurtosis<T>,
        cons<WeightedSum<T>,
        cons<WeightedMin<T>,
        cons<WeightedMax<T>,
             nil> > > > > > > > > > > > > > > > > > > > > >
    type;
};

template <template<class, class> class BASE, class T, template<class> class L>
struct PythonUseTemplateList
    : public BASE<T, typename L<T>::type>
{};

template <class AccuType, unsigned N>
struct PythonObjectFeatures
    : public PythonUseTemplateList<NestedAccumulators,
                                   StridePairPointer<N, AccuType>,
                                   PythonObjectFeaturesList>
{};

template <unsigned N>
struct PyCoordinatePermutation
{
    typedef typename detail::ExtractorInfo ExtractorInfo;
    typedef TinyVector<unsigned, N>        perm_type;
    bool      is_valid;
    perm_type perm;

    PyCoordinatePermutation() : is_valid(false) {}

    template <class Array>
    PyCoordinatePermutation(const Array & array)
        : is_valid(true)
    {
        perm_type p;
        linearSequence(p.begin(), p.end(), 0);
        p = array.permuteLikewise(p);
        inversePermutation(p.begin(), p.end(), perm.begin());
    }
    PyCoordinatePermutation & operator=(const PyCoordinatePermutation & x)
    {
        vigra_precondition(x.is_valid,
            "PythonAccumulators: assigning undefined permutation.");
        if (is_valid)
        {
            vigra_precondition(perm == x.perm,
                "PythonAccumulators: using incompatible axes orders.");
        }
        else
        {
            is_valid = x.is_valid;
            perm = x.perm;
        }
        return *this;
    }
    void check() const
    {
        vigra_precondition(is_valid,
            "PythonAccumulators::extract(): using undefined axes order.");
    }

    template <class V>
    V permuteVector(const V & v, ExtractorInfo::Flags flags)
    {
        V ret;
        if (is_valid && (flags & ExtractorInfo::IsCoordinate))
            applyPermutation(perm.begin(), perm.end(), v.begin(), ret.begin());
        else
            ret = v;
        return ret;
    }
};

template <class NUMPYARRAY>
void registerConvertersDummy(NUMPYARRAY &) {}
template <class NUMPYARRAY>
bool registerConvertersOnce(bool)
{
    registerConverters(&registerConvertersDummy<NUMPYARRAY>);
    return true; // employ "register only on first use" semantics
}

template <class EXBASE, unsigned N>
struct PythonExtractor : public EXBASE
{
    typedef typename detail::ExtractorInfo ExtractorInfo;
    typedef PyCoordinatePermutation<N>     perm_type;

    python::object val;
    perm_type      perm;

    operator python::object() const { return val; }
    
    PythonExtractor() {}
    PythonExtractor(const perm_type & p) : perm(p) {}
    
    // for count or AccuType
    void operator()(const double & v, ExtractorInfo::Flags)
    {
        val = python::object(v);
    }
    // for AccuType, most likely
    void operator()(const float & v, ExtractorInfo::Flags)
    {
        val = python::object(v);
    }
    // (weighted) coordinates:
    void operator()(const TinyVector<double, N> & v,
                     ExtractorInfo::Flags flags)
    {
        static bool reg = registerConvertersOnce<TinyVector<double, N> >(reg);
        val = python::object(perm.permuteVector(v, flags));
    }
    // e.g., octet histogram
    void operator()(const TinyVector<double, 256> & v,
                     ExtractorInfo::Flags flags)
    {
        typedef NumpyArray<1, double> val_type;
        static bool reg = registerConvertersOnce<val_type>(reg);
        typedef typename val_type::difference_type shape_type;
        val_type x(shape_type(256));
        std::copy(v.begin(), v.end(), x.begin());
        val = python::object(x);
    }
};

template <class T>
struct PyToFloat
{
    typedef float type;
};
template <class X, int N>
struct PyToFloat<TinyVector<X, N> >
{
    typedef TinyVector<float, N> type;
};

template <class EXBASE, unsigned N>
struct PythonArrayExtractor : public EXBASE
{
    typedef typename detail::ExtractorInfo ExtractorInfo;
    typedef PyCoordinatePermutation<N>     perm_type;

    unsigned       i;
    unsigned       size;
    python::object val; // The use of python::object for val is slow inasmuch
                        // it provokes reference count adjustments for each
                        // call of operator(). A proper replacement would be
                        // a type not unlike boost::any, with the added
                        // functionality of conversion to Python objects.
    perm_type      perm;

    PythonArrayExtractor(unsigned s, const perm_type & p)
        : i(0), size(s), perm(p) {}

    operator python::object() const { return val; }
    
    template <class T>
    void put(const T & v)
    {
        typedef typename PyToFloat<T>::type F;
        typedef typename
                If<typename NumericTraits<T>::isScalar, Singleband<F>, F>::type
            single_type;
        typedef NumpyArray<1, single_type> val_type;
        static bool reg = registerConvertersOnce<val_type>(reg);
        typedef typename val_type::difference_type shape_type;
        if (i == 0)
            val = python::object(val_type(shape_type(size)));
        val_type z = val_type(val.ptr());
        z(i++) = v;
    }

    // for count or AccuType
    void operator()(const double & v, ExtractorInfo::Flags)
    {
        put(float(v));
    }
    // for AccuType, most likely
    void operator()(const float & v, ExtractorInfo::Flags)
    {
        put(v);
    }
    // (weighted) coordinates:
    void  operator()(const TinyVector<double, N> & v,
                     ExtractorInfo::Flags flags)
    {
        put(perm.permuteVector(v, flags));
    }
    // e.g., octet histogram
    void operator()(const TinyVector<double, 256> & v, ExtractorInfo::Flags)
    {
        put(v);
    }
};

// convenience functions for sequences
template <class X>
void select_any(X & x, python::str t)
{
    x.select(t);
}
template <class X, class SEQ>
void select_any(X & x, SEQ seq)
{
    for (unsigned i = 0; i != len(seq); ++i)
        x.select(py_str(seq[i]));
}

template <class X, class SEQ>
bool isSelected_seq(const X & x, SEQ seq)
{
    for (unsigned i = 0; i != len(seq); ++i)
        if (!x.isSelected(py_str(seq[i])))
            return false;
    return true;
}

template <class X, class SEQ>
python::dict extract_seq(const X & x, SEQ seq)
{
    python::dict dict;
    for (unsigned i = 0; i != len(seq); ++i)
    {
        python::str t = py_str(seq[i]);
        if (x.isSelected(t))
            dict[t] = x.extract(t);
    }
    return dict;
}

template <class A>
struct PythonFeatureHash;

template <>
struct PythonFeatureHash<python::str>
    : public std::unary_function<std::string, std::size_t>
{
    std::size_t operator()(const python::str & x) const
    {
        return PyObject_Hash(x.ptr());
    }
};

template <class A, class B>
struct PythonFeatureMap
    : public boost::unordered_map<A, B, PythonFeatureHash<A> >
{
    std::string name_convert(const A & x) const
    {
        // produce nice error messages within Python
        return python::extract<std::string>(x);
    }
    void reserve(unsigned size) { this->rehash(3 * size); }
};

template <class ACX>
struct PythonDispatch
   : public detail::selector_name_dispatch<ACX, python::str, PythonFeatureMap>
{};


// employ construct on first use to avert the static initialisation order fiasco
// for shared libraries:
template <class ACX>
const PythonDispatch<ACX> & PythonSelector()
{
    static const PythonDispatch<ACX>* p = new PythonDispatch<ACX>();
    return *p;
}

template <class VoxelType, unsigned N, template<class, unsigned> class Features>
struct PythonAccumulatorsBase
    : public Features<VoxelType, N>
{
    typedef NumpyArray<N, Singleband<VoxelType> >  volume_type;
    typedef PythonDispatch<PythonAccumulatorsBase> dispatch_base;
    typedef const dispatch_base &                  dispatch_type;
    typedef typename dispatch_base::selector_type  selector_type;
    typedef PyCoordinatePermutation<N>             permutation_type;

    static const unsigned size = PythonAccumulatorsBase::feature_type::size;

    void select_base(python::str t, dispatch_type dispatch)
    {
        dispatch.select(t, *this);
    }
    bool isSelected_base(python::str t, dispatch_type dispatch) const
    {
        return dispatch.isSelected(t, *this);
    }
    unsigned numberOfPasses_str_base(python::str t, dispatch_type dispatch)
        const
    {
        return dispatch.numberOfPasses(t, *this);
    }

    static python::list names()
    {
        python::str input[size];
        PythonSelector<PythonAccumulatorsBase>().names(input);
        python::list ret;
        for (unsigned i = 0; i != size; ++i)
            ret.append(input[i]);
        ret.sort();
        return ret;
    }
    python::list selected_base(dispatch_type dispatch) const
    {
        python::list seq = names();
        python::list ret;
        for (unsigned i = 0; i != len(seq); ++i)
        {
            python::str z = py_str(seq[i]);
            if (isSelected_base(z, dispatch))
                ret.append(z);
        }
        return ret;
    }
};

template <class VoxelType, unsigned N, template<class, unsigned> class Features>
struct PythonAccumulators
    : public PythonAccumulatorsBase<VoxelType, N, Features>
{
    typedef PythonAccumulatorsBase<VoxelType, N, Features> base_type;
    typedef typename PythonAccumulators::dispatch_type     dispatch_type;
    typedef typename PythonAccumulators::volume_type       volume_type;
    typedef typename PythonAccumulators::extractor_type    extractor_base;
    typedef typename PythonAccumulators::permutation_type  permutation_type;
    typedef PythonExtractor<extractor_base, N>             extractor_type;

    static const unsigned dim = N;
    static const unsigned size = PythonAccumulators::feature_type::size;

    dispatch_type     dispatch;
    permutation_type  my_permutation;
    
    PythonAccumulators() : dispatch(PythonSelector<base_type>()) {}

    void select(python::str t)
    {
        this->select_base(t, dispatch);
    }
    bool isSelected(python::str t) const
    {
        return this->isSelected_base(t, dispatch);
    }
    python::list selected() const
    {
        return this->selected_base(dispatch);
    }
    unsigned numberOfPasses_str(python::str t) const
    {
        return this->numberOfPasses_str_base(t, dispatch);
    }

    python::object extract(python::str t) const
    {
        extractor_type g(my_permutation);
        if (isSelected(t))
            dispatch.extract(t, *this, g);
        return g;
    }

    void inspect(volume_type volume)
    {
        my_permutation = permutation_type(volume);
        inspectMultiArray(srcCoordinateMultiArrayRange(volume), *this);
    }

    template <class ANY>
    void inspect_any(volume_type volume, ANY any)
    {
        select_any(*this, any);
        inspect(volume);
    }

    python::object inspect_extract_str(volume_type volume, python::str t)
    {
        inspect_any(volume, t);
        return extract(t);
    }
    template <class SEQ>
    python::dict inspect_extract_seq(volume_type volume, SEQ seq)
    {
        inspect_any(volume, seq);
        return extract_seq(*this, seq);
    }
    python::dict inspect_extract(volume_type volume)
    {
        return inspect_extract_seq(volume, selected());
    }
    python::dict extract_selected()
    {
        return extract_seq(*this, selected());
    }
        /** merge single accumulator
        */
    void merge(const PythonAccumulators & x)
    {
        // only merge the states of being selected from empty accumulators
        if (x.my_permutation.is_valid)
        {
            my_permutation = x.my_permutation;
            (*this)(x); // includes selectSelected()
        }
        else
        {
            selectSelected(x);
        }
    }
};

struct PythonAccumulatorsArrayExtender
{
    template <class AORS>
    static void call(AORS & x, typename AORS::iterator k)
    {
        for (; k != x.end(); ++k)
            k->selectSelected(*x.begin());
    }
};

template <class VoxelType, unsigned N, template<class, unsigned> class Features>
struct PythonAccumulatorsArray
{
    typedef PythonAccumulatorsBase<VoxelType, N, Features> base_type;
    typedef typename base_type::dispatch_type              dispatch_type;
    typedef typename base_type::volume_type                volume_type;
    typedef typename base_type::extractor_type             extractor_base;
    typedef typename base_type::permutation_type           permutation_type;
    typedef PythonArrayExtractor<extractor_base, N>        extractor_type;

    typedef NumpyArray<N, Singleband<npy_uint32> >         labels_type;

    typedef ArrayOfRegionStatistics<base_type, npy_uint32, true,
                                    PythonAccumulatorsArrayExtender>
        stat_type;
    typedef typename stat_type::iterator                   iterator;
    typedef typename stat_type::const_iterator             const_iterator;

    static const unsigned dim = N;
    static const unsigned size = base_type::feature_type::size;

    dispatch_type    dispatch;
    stat_type        objects;
    permutation_type my_permutation;

    PythonAccumulatorsArray(unsigned max_region_label = 0) // at least 1 accu..
        : dispatch(PythonSelector<base_type>()), objects(max_region_label) {}

    static python::list names() { return base_type::names(); }
    unsigned labelCount() const { return objects.size(); }

    const_iterator begin() const { return objects.begin(); }
          iterator begin()       { return objects.begin(); }
    const_iterator end()   const { return objects.end(); }
          iterator end()         { return objects.end(); }

    void sweep_selected(iterator k)
    {
        stat_type::extender_type::call(objects, k);
    }
    void sweep_selected()
    {
        sweep_selected(begin());
    }
    void select(python::str t)
    {
        begin()->select_base(t, dispatch);
        sweep_selected();
    }
    bool isSelected(python::str t) const
    {
        return begin()->isSelected_base(t, dispatch);
    }
    python::list selected() const
    {
        return begin()->selected_base(dispatch);
    }
    unsigned numberOfPasses_str(python::str t) const
    {
        return begin()->numberOfPasses_str_base(t, dispatch);
    }

    python::object extract(python::str t) const
    {
        extractor_type g(labelCount(), my_permutation);
        typename base_type::selector_type selector = dispatch.ref(t);
        if (isSelected(t))
            for (const_iterator k = begin(); k != end(); ++k)
                selector.extract(*k, g);
        return g;
    }

    void inspect(volume_type volume, labels_type labels)
    {
       my_permutation = permutation_type(volume);
       inspectTwoMultiArrays(srcCoordinateMultiArrayRange(volume),
                             srcMultiArray(labels),
                             objects);
    }

    template <class ANY>
    void inspect_any(volume_type volume, labels_type labels, ANY any)
    {
        select_any(*this, any);
        inspect(volume, labels);
    }

    python::object inspect_extract_str(volume_type volume, labels_type labels,
                                       python::str t)
    {
        inspect_any(volume, labels, t);
        return extract(t);
    }
    template <class SEQ>
    python::dict inspect_extract_seq(volume_type volume, labels_type labels,
                                     SEQ seq)
    {
        inspect_any(volume, labels, seq);
        return extract_seq(*this, seq);
    }
    python::dict inspect_extract(volume_type volume, labels_type labels)
    {
        return inspect_extract_seq(volume, labels, selected());
    }
    python::dict extract_selected()
    {
        return extract_seq(*this, selected());
    }

    void size_test(unsigned label) const
    {
        vigra_precondition(label < objects.size(),
            "PythonAccumulatorsArray: label value too large.");
    }
        /** merge second region into first
        */
    void merge_labels(unsigned label1, unsigned label2)
    {
        if (label1 == label2)
            return;
        size_test(label1);
        size_test(label2);
        objects.merge(label1, label2);
        // clear merged object
        objects[label2] = base_type();
        // restore selected features correctly
        objects[label2].selectSelected(objects[label1]);
    }
        /** merge all accumulators
        */
    void merge(const PythonAccumulatorsArray & x)
    {
        // only merge the states of being selected from empty accumulators
        if (x.my_permutation.is_valid)
        {
            my_permutation = x.my_permutation;
            unsigned old_size = labelCount();

            objects.extend(x.labelCount() - 1); // iterator invalidation
            iterator k = begin();
            for (const_iterator q = x.begin(); q != x.end(); ++q, ++k)
                (*k)(*q); // includes selectSelected()

            sweep_selected(begin() + old_size);
        }
        else
        {
            begin()->selectSelected(*x.begin());
            sweep_selected();
        }
    }
    unsigned numberOfPasses()
    {
        return objects.numberOfPasses();
    }
};

// The overloaded set of free-standing "extractFeatures" Python functions
// must somehow create instances of wrapped C++ classes. However, returning
// naked pointers from the C++ functions corresponding to "extractFeatures" is
// not an option, since that would violate exception safety. It is thus easiest
// to use the Python interpreter's mechanisms to manage object lifetime.

// Construct on first use again: for global functions returning wrapped classes.
template <class CLASS>
python::class_<CLASS> & PythonClassSingleton(const std::string class_name = "",
                                             const std::string class_info = "",
                                             const std::string init_info  = "")
{
    using namespace python;
    typedef class_<CLASS> class_wrapper;
    static class_wrapper* p = new class_wrapper(class_name.c_str(),
                                                class_info.c_str(),
                                                init<>(init_info.c_str()));
    return *p;
}
template <class CLASS>
python::object PythonClassInstance()
{
    return python::object(PythonClassSingleton<CLASS>())();
}

template <class CLASS>
struct PyAccumulatorExtraDefs
{
    typedef typename CLASS::volume_type volume_type;

    static std::string class_name() { return "Accumulators"; }
    static std::string class_type() { return "accumulator"; }

    template <class DEF>
    static void exec(DEF & defined_class, const std::string & class_name) {}

    template <class ANY>
    static python::object inspect_any(volume_type volume, ANY any)
    {
        using namespace python;
        object ret = PythonClassInstance<CLASS>();
        extract<CLASS &>(ret)().inspect_any(volume, any);
        return ret;
    }
};

template <class C, class KEYWORDS>
void defineObjectFeaturesClass(const KEYWORDS &    images,
                               const std::string & what,
                               const std::string & name_ext)
{
    using namespace python;
    docstring_options doc_options(true, true, false);
    typedef PyAccumulatorExtraDefs<C> extras;

    std::ostringstream dim_s;
    dim_s << C::dim;
    const std::string dim = dim_s.str();
    const std::string class_name
        = extras::class_name() + dim + "D" + name_ext;
    const std::string class_type = extras::class_type();
    const arg & key = arg(what.c_str());
    const std::string whats = what + "s";
    const arg & keys = arg(whats.c_str());
    
    class_<C> & defined_class
        = PythonClassSingleton<C>(class_name,
            dim + "-dimensional " + class_type + "s.\n\n"
            "For more details, see the C++ documentation.\n\n",
            "Standard constructor::\n\n   " + class_name + "()\n\n"
                    "Creates an empty " + class_type + " object.\n");
    defined_class
        .def(init<C>(args(class_type.c_str()),
             ("Copy constructor::\n\n"
             "   " + class_name + "(other_accumulator)\n\n").c_str()))
        .def("inspect", registerConverters(&C::inspect), images)
        .def("inspect", &C::template inspect_any<str>,   (images, key))
        .def("inspect", &C::template inspect_any<list>,  (images, keys))
        .def("inspect", &C::template inspect_any<tuple>, (images, keys))
        .def("extract",     &C::extract,           key)
        .def("extract",     extract_seq<C, list>,  keys)
        .def("extract",     extract_seq<C, tuple>, keys)
        .def("__getitem__", &C::extract,           key)
        .def("__getitem__", extract_seq<C, list>,  keys)
        .def("__getitem__", extract_seq<C, tuple>, keys)
        .def("extract", &C::inspect_extract_str,   (images, key))
        .def("extract", &C::template inspect_extract_seq<list>,  (images, keys))
        .def("extract", &C::template inspect_extract_seq<tuple>, (images, keys))
        .def("extract", &C::inspect_extract, images)
        .def("extract", &C::extract_selected)
        .def("select", &C::select, key)
        .def("select", select_any<C, list>,  keys)
        .def("select", select_any<C, tuple>, keys)
        .def("isSelected", &C::isSelected, key)
        .def("isSelected", isSelected_seq<C, list>,  keys)
        .def("isSelected", isSelected_seq<C, tuple>, keys)
        .def("selected", &C::selected)
        .def("merge", &C::merge, arg("other_accumulator"))
        .def("numberOfPasses", &C::numberOfPasses_str, key)
        .def("numberOfPasses", &C::numberOfPasses)
        .def("names", &C::names).staticmethod("names")
        ;
    extras::exec(defined_class, class_name);

    // global functions
    if (!name_ext.empty())
        return;
    def("extractFeatures", extras::template inspect_any<str>,   (images, key));
    def("extractFeatures", extras::template inspect_any<list>,  (images, keys));
    def("extractFeatures", extras::template inspect_any<tuple>, (images, keys));
}

template <class T, unsigned N, template<class, unsigned> class Features>
struct PyAccumulatorExtraDefs<PythonAccumulatorsArray<T, N, Features> >
{
    typedef PythonAccumulatorsArray<T, N, Features> C;
    typedef typename C::volume_type volume_type;
    typedef typename C::labels_type labels_type;

    static std::string class_name() { return "AccumulatorsArray"; }
    static std::string class_type() { return "accumulator array"; }

    template <class DEF>
    static void exec(DEF & defined_class, const std::string & class_name)
    {
        using namespace python;
        docstring_options doc_options(true, true, false);
        defined_class
            .def(init<unsigned>(args("max_region_label"),
                ("Constructor with reservation of maximum region label::\n\n"
                "   " + class_name + "(max_region_label)\n\n").c_str()))
            .def("labelCount", &C::labelCount)
            .def("merge", &C::merge_labels, (arg("label1"), arg("label1")))
            ;
    }

    template <class ANY>
    static python::object inspect_any(volume_type volume, labels_type labels,
                                      ANY any)
    {
        using namespace python;
        object ret = PythonClassInstance<C>();
        extract<C &>(ret)().inspect_any(volume, labels, any);
        return ret;
    }
};

const char *const what = "feature_name";
const char *const data[] = { "point", "line", "image", "volume" };

template <class DATA, unsigned dim, template<class, unsigned> class Features>
void defineObjectFeaturesAccumulators(const std::string & name_ext)
{
    using namespace python;
    defineObjectFeaturesClass<PythonAccumulators<DATA, dim, Features> >(
        arg(data[dim]), what, name_ext);
}
template <class DATA, unsigned dim, template<class, unsigned> class Features>
void defineObjectFeaturesAccumulatorsArray(const std::string & name_ext)
{
    using namespace python;
    defineObjectFeaturesClass<PythonAccumulatorsArray<DATA, dim, Features> >(
        (arg(data[dim]), arg("labels")), what, name_ext);
}

template <class DATA, template<class, unsigned> class Features>
void defineObjectFeaturesAny(const std::string & name_ext = "")
{
    defineObjectFeaturesAccumulators     <DATA, 2, Features>(name_ext);
    defineObjectFeaturesAccumulators     <DATA, 3, Features>(name_ext);
    defineObjectFeaturesAccumulatorsArray<DATA, 2, Features>(name_ext);
    defineObjectFeaturesAccumulatorsArray<DATA, 3, Features>(name_ext);
}

void defineObjectFeatures()
{
    defineObjectFeaturesAny<float, PythonObjectFeatures>();
    // defineObjectFeaturesAny<float, PythonObjectFeatures_extra>("_extra");
}

} // namespace vigra
