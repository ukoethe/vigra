#ifndef VIGRA_OBJECT_FEATURES_HXX
#define VIGRA_OBJECT_FEATURES_HXX

#include <cmath>
#include <string>
#include <utility>

#include "error.hxx"
#include "numerictraits.hxx"
#include "functortraits.hxx"
#include "type_lists.hxx"
#include "algorithm.hxx"
#include "mathutil.hxx"
#include "tinyvector.hxx"
#include "sized_int.hxx"

namespace vigra {

// mask cl.exe shortcomings [begin]
#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4503 )
#endif

namespace detail {

template <unsigned n> // test the number of passes a feature requires
struct passes
{
    template <class X>
    struct at_least
    {
        static const bool value = X::number_of_passes >= n;
    };
};

template <class ACX>
struct passes_number
{
    static const unsigned value = ACX::number_of_passes;
};

template <class, unsigned n>
struct passes_test
{
    static bool call(unsigned pass)
    {
        return pass == n;
    }
};
// compile-time optimisation for completely 2-pass accumulators:
template <class ACX>
struct passes_test<ACX, 2>
{
    static bool call(unsigned pass)
    {
        return ACX::max_passes == 2 || pass == 2;
    }
};
template <class ACX, unsigned n>
inline bool test_pass(unsigned pass)
{
    return passes_test<ACX, n>::call(pass);
}

template <class X>
struct get_result_type
{
    typedef typename X::result_type type;
};

template <class RTL> // build up a base class over a result-type list
struct extractor_base : public extractor_base<typename RTL::rest>
{
    using extractor_base<typename RTL::rest>::operator();
    virtual void operator()(const typename RTL::first &, ExtractorInfo::Flags)
        = 0;
};
template <>
struct extractor_base<type_lists::nil>
{
    virtual void operator()() {}
};

template <class T>
struct features_double_promote
{
    typedef typename NumericTraits<T>::DoublePromote TraitsDoublePromote;
    typedef typename
        IfBool<
            IsSameType<
                Error_NumericTraits_not_specialized_for_this_case,
                TraitsDoublePromote>::value,
            T,
            TraitsDoublePromote>::type
        type;
};

template <class V, ExtractorInfo::Flags flags = ExtractorInfo::None>
struct feature_base
{
        /** the functor's result type
        */
    typedef V result_type;
    static const unsigned number_of_passes = 0;
protected:
    mutable result_type value;
public:
        /** retrieve calculated or accumulated value
        */
    const result_type & soft_value() const
    {
        return this->value;
    }
        /** dummy calculation for accumulators, called before extract()
        */
    template <class ACX>
    void calculate(const ACX &) const {}
        /** copy calculated or accumulated value to extractor object
        */
    template <class VG>
    void extract(VG & g) const
    {
        g(this->value, flags);
    }
        /** (re-)init
        */
    void reset()
    {
        this->value = result_type(); // guaranteed to be zero for built-in types
    }
    template <class X>
    result_type cast(const X & x) const
    {
        return static_cast<result_type>(x);
    }
    template <class X>
    void set_val(const X & x)
    {
        this->value = static_cast<const result_type &>(x);
    }
    template <class X>
    void add_val(const X & x)
    {
        this->value += static_cast<const result_type &>(x);
    }
    // variant for calculators
    template <class X>
    void set_val(const X & x) const
    {
        this->value = static_cast<const result_type &>(x);
    }
    feature_base()
    {
        reset();
    }
    feature_base(const V & v) : value(v) {}
    // provide a default static name just in case.
    static std::string name() { return "[invalid name]"; }
};

// all collectors / accumulators must implement void reset()

template <class T, class V, ExtractorInfo::Flags flags = ExtractorInfo::None>
struct collector_base : public feature_base<V, flags>
{
    collector_base() {}
    collector_base(const V & v) : feature_base<V>(v) {}

    static const unsigned number_of_passes = 1;

        /** the functor's argument type
        */
    typedef T argument_type;

        /** the functor's first argument type (for calls with a weight)
        */
    typedef T first_argument_type;

        /** the functor's second argument type (for calls with a weight)
        */
    typedef double second_argument_type;

    typedef typename collector_base::result_type result_type;
    typedef          collector_base              collector_base_type;

        /** retrieve accumulated value
        */
    const result_type & get_value() const
    {
        return this->value;
    }
        /** automated re-init
        */
    void operator()()
    {
        this->reset();
    }
    void fail(const std::string & what)
    {
        vigra_fail("collector_base::operator(): " +
                    what + " is not implemented for this accumulator.");
    }
    template <class ACX>
    void updatePass1(ACX &, const T &) { fail("regular update"); }
    template <class ACX>
    void updatePass1(ACX &, const T &, double) { fail("weighted update"); }
    template <class ACX>
    void merge(const ACX &, const ACX &, const ACX &) { fail("merge"); }
    // Multi-pass features template classes must
    // a) declare their number of passes (@n) as
    //    "static const unsigned number_of_passes = @n;"
    // b) implement the routine for passes >= 2, namely,
    //    template <class ACX> void call(ACX & x, const T & v, unsigned m),
    //    which must dispatch different pass numbers m by itself, such as
    //    using "if (test_pass<ACX, 2>(m))" for pass 2 -- even if there is
    //    only pass 2 to handle, since calls with higher pass number 'm' may
    //    occur, for example due to user-written accumulators.
    //
    // The multi-pass 'inspect' algorithms are supposed to call
    // calc_sync() on the accumulator objects between passes,
    // thus it is possible to effectively use pre-calculated values of
    // calculators in call() via
    // "use<feature_template_name>::c_val(x.calculators)", whereas the use
    // of 'val' would imply a (likely expensive) re-calculation for
    // each update.
    template <class ACX>
    void updatePassN(ACX &, const T &, unsigned) { fail("pass n update"); }
};

template <class T, class A = T,
          ExtractorInfo::Flags flags = ExtractorInfo::None>
struct numeric_accumulator_base
    : public collector_base<T, typename features_double_promote<A>::type, flags>
{
    typedef typename numeric_accumulator_base::collector_base_type
        collector_base_type;
    typedef typename numeric_accumulator_base::argument_type
        argument_type;
    typedef typename numeric_accumulator_base::first_argument_type
        first_argument_type;
    typedef typename numeric_accumulator_base::second_argument_type
        second_argument_type;
    typedef typename numeric_accumulator_base::result_type
        result_type;

    typedef numeric_accumulator_base accumulator_base_type;
};

template <class T, ExtractorInfo::Flags flags = ExtractorInfo::None>
struct calculator_base
    : public feature_base<typename features_double_promote<T>::type, flags>
{
    // the actual calculator class templates should re-implement
    // template <class ACX> void calculate(const ACX & x) const,
    // wherein the values of other features are to be accessed
    // as 'soft values' via "use<feature_template_name>::val(x)"
};

// The 'use' and 'use_d' helper classes encapsulate the verbose
// 'template' disambiguators that all feature implementations
// would be forced to employ otherwise.

template <template<class> class V>
struct use
{
    template <class ACX>
    static const typename ACX::template get<V>::result_type &
    f_val(const ACX & x)
    {
        return x.template value_template<V>();
    }
    template <class ACX>
    static void
    set(ACX & x)
    {
        x.template set_template<V>();
    }
    template <class ACX>
    static bool
    is_set(const ACX & x)
    {
        return x.template is_set_template<V>();
    }
    // "soft value" for merging and output / calculation purposes:
    // return "empty" for inexistent tuple elements
    template <class ACX>
    static const typename ACX::template get<V>::result_type &
    val(const ACX & x)
    {
        x.template calculate_template<V>(x);
        typedef typename ACX::template get<V>::result_type type;
        static type t;
        if (is_set(x))
            return x.template soft_value_template<V>();
        else
            return t;
    }
    // pre-calculated (such as via calc_sync()) values of calculators.
    template <class CCX>
    static const typename CCX::template get<V>::result_type &
    c_val(const CCX & x)
    {
        return x.template soft_value_template<V>();
    }
    // direct access to accumulators or calculators as a means of last resort.
    template <class ACX>
    static const typename ACX::template get<V>::type &
    member(const ACX & x)
    {
        return x.template member_template<V>();
    }
    template <class ACX>
    static       typename ACX::template get<V>::type &
    member(      ACX & x)
    {
        return x.template member_template<V>();
    }
};
template <class V>
struct use_d
{
    template <class ACX>
    static const typename V::result_type &
    f_val(const ACX & x)
    {
        return x.template value<V>();
    }
    template <class ACX>
    static bool
    is_set(const ACX & x)
    {
        return x.template is_set<V>();
    }
    template <class ACX>
    static const typename V::result_type &
    val(const ACX & x)
    {
        x.template calculate<V>(x);
        typedef typename V::result_type type;
        static type t;
        if (is_set(x))
            return x.template soft_value<V>();
        else
            return t;
    }
};

// accessor classes that are specialised for StridePairPointer<N, T>
template<class T>
struct accumulable_value_access
{
    typedef T accumulable_type;
    typedef T type;
    static const T & get(const T & v) { return v; }
    static const ExtractorInfo::Flags flags = ExtractorInfo::None;
};
template<class> // undefined class to provoke usable error messages
struct Error_use_StridePairPointer_for_Accumulator__coordinates_absent_in_type_;
// arbitrary accumulated types have no embedded coordinates: flag an error.
template<class T>
struct accumulable_coord_access
{
    typedef T accumulable_type;
    typedef T type;
    static const T & get(const T &)
    {
        return 
        Error_use_StridePairPointer_for_Accumulator__coordinates_absent_in_type_
            <T>();
    }
    static const ExtractorInfo::Flags flags = ExtractorInfo::IsCoordinate;
};
template<class T>
struct accumulable_weighted_access : public accumulable_coord_access<T> {};

} // namespace detail

// -- start feature template classes --

namespace acc {

template <class T>
struct Count : public detail::collector_base<T, double>
{
    /** update
     */
    template <class ACX>
    void updatePass1(ACX &, const T &)
    {
        this->add_val(1);
    }
    /** update using weighted input
     */
    template <class ACX>
    void updatePass1(ACX & x, const T &, double weight)
    {
        this->add_val(weight);
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        using namespace detail;
        this->set_val(use_d<Count>::val(x) + use_d<Count>::val(y));
    }
    static std::string name() { return "count"; }
};

} // namespace acc

namespace detail {

template <class T, template<class> class CV>
struct gen_accumulator_base
    : public numeric_accumulator_base<T, typename CV<T>::type, CV<T>::flags>
{
    typedef typename CV<T>::type type;
    type get(const T & v) { return CV<T>::get(v); }
};

template <class T>
struct value_accumulator_base
    : public gen_accumulator_base<T, accumulable_value_access> {};

template <class T>
struct coord_accumulator_base
    : public gen_accumulator_base<T, accumulable_coord_access> {};

template <class T>
struct weighted_accumulator_base
    : public gen_accumulator_base<T, accumulable_weighted_access> {};


template <class T, template<class> class CV,
                   template<class> class SELF,
                   template<class> class SUM_SQUARED_DIFF,
                   template<class> class COVAR = SUM_SQUARED_DIFF>
struct gen_mean_acc : public gen_accumulator_base<T, CV>,
                      public type_lists::uses_template<T, acc::Count>
                                       ::template follows<SUM_SQUARED_DIFF,
                                                          COVAR>
{
    /** update average
    */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v)
    {
        this->add_val((this->get(v) - this->value)
                      / use<acc::Count>::f_val(accumulators));
    }
    /** update average, using weighted input.
     * <tt>stats(value, 1.0)</tt> is equivalent to the unweighted
     * call <tt>stats(value)</tt>, and <tt>stats(value, 2.0)</tt>
     * is equivalent to two unweighted calls.
     */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v, double weight)
    {
        this->add_val((this->get(v) - this->value) * weight
                      / use<acc::Count>::f_val(accumulators));
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        this->set_val((   use<acc::Count>::val(x) * use<SELF>::val(x)
                        + use<acc::Count>::val(y) * use<SELF>::val(y)
                      ) / use<acc::Count>::val(z));
    }
};

template <class T, template<class> class SELF,
                   template<class> class MEAN>
struct gen_sum_squared_diff_merge
{
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        double count_x = use<acc::Count>::val(x);
        double count_y = use<acc::Count>::val(y);
        static_cast<SELF<T>*>(this)->set_val(
            use<SELF>::val(x) + use<SELF>::val(y)
            + count_x * count_y / (count_x + count_y)
                      * power<2>(use<MEAN>::val(y) - use<MEAN>::val(x)));
    }
};

template <class T, template<class> class CV,
                   template<class> class SELF,
                   template<class> class MEAN>
struct gen_sum_squared_diff_acc
    : public gen_accumulator_base<T, CV>,
      public gen_sum_squared_diff_merge<T, SELF, MEAN>,
      public type_lists::uses_template<T, acc::Count>
                       ::template implies<MEAN>
{
    typedef typename gen_sum_squared_diff_acc::result_type result_type;
    using gen_sum_squared_diff_merge<T, SELF, MEAN>::merge;

    /** update sum of squared differences
    */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v)
    {
        result_type t1 = static_cast<result_type>(
                         this->get(v) - use<MEAN>::f_val(accumulators));
        result_type t2 = static_cast<result_type>(
                         t1 / use<acc::Count>::f_val(accumulators));
        this->add_val(t1 * t2 * (use<acc::Count>::f_val(accumulators) - 1.0));
    }
    /** update sum of squared differences, using weighted input.
     * <tt>stats(value, 1.0)</tt> is equivalent to the unweighted
     * call <tt>stats(value)</tt>, and <tt>stats(value, 2.0)</tt>
     * is equivalent to two unweighted calls.
     */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v, double weight)
    {
        result_type t1 = this->get(v) - use<MEAN>::f_val(accumulators);
        result_type t2 = t1 * weight / use<acc::Count>::f_val(accumulators);
        // count > weight is ensured by preceding acc::Count::operator() call
        this->add_val(t1 * t2
                      * (use<acc::Count>::f_val(accumulators) - weight));
    }
};

} // namespace detail

template <class V>
struct covariance_result
{
    typedef V type;
};

namespace detail {

// covariance_result<V>:
// if is_derived V of TinyVector<X, N>: (StridePairCoord, RGBValue, whatever...)
//    typedef TinyVector<X, N * (N - 1)> type;
//
/////   ->  new is_derived meta with X, N template parameters.....

template <class T, template<class> class CV,
                   template<class> class SELF,
                   template<class> class MEAN>
struct gen_covariance_acc
    : public numeric_accumulator_base<T, typename
                                         covariance_result<typename
                                                           CV<T>::type>::type,
                                         CV<T>::flags>,
      public type_lists::uses_template<T, acc::Count>
                       ::template implies<MEAN>
{
    typedef typename gen_covariance_acc::result_type result_type;

    /** update sum of squared differences
    */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v)
    {
//         result_type t1 = this->get(v) - use<MEAN>::f_val(accumulators);
//         result_type t2 = t1 / use<acc::Count>::f_val(accumulators);
//         this->add_val(t1 * t2 * (use<acc::Count>::f_val(accumulators) - 1.));
    }
    /** update sum of squared differences, using weighted input.
     * <tt>stats(value, 1.0)</tt> is equivalent to the unweighted
     * call <tt>stats(value)</tt>, and <tt>stats(value, 2.0)</tt>
     * is equivalent to two unweighted calls.
     */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v, double weight)
    {
//         result_type t1 = this->get(v) - use<MEAN>::f_val(accumulators);
//         result_type t2 = t1 * weight / use<acc::Count>::f_val(accumulators);
//         // count > weight is ensured by preceding acc::Count::operator() call
//         this->add_val(t1 * t2
//                       * (use<acc::Count>::f_val(accumulators) - weight));
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        // ...
    }
};

} // namespace detail

namespace acc {

// mean and second moment for values
template <class T>
struct Moment2;

template <class T>
struct Mean
    : public detail::gen_mean_acc<T, detail::accumulable_value_access, Mean,
                                     Moment2, Moment2>
{
    static std::string name() { return "mean"; }
};
template <class T>
struct Moment2
    : public detail::gen_sum_squared_diff_acc<T,
                                              detail::accumulable_value_access,
                                              Moment2, Mean
                                             >
{
    static std::string name() { return "moment2"; }
};

// mean and second moment for coordinates

template <class T>
struct CoordMoment2;

template <class T>
struct CoordMean
    : public detail::gen_mean_acc<T, detail::accumulable_coord_access,
                                     CoordMean, CoordMoment2, CoordMoment2>
{
    static std::string name() { return "coord mean"; }
};

template <class T>
struct CoordMoment2
    : public detail::gen_sum_squared_diff_acc<T,
                                              detail::accumulable_coord_access,
                                              CoordMoment2, CoordMean
                                             >
{
    static std::string name() { return "coord moment2"; }
};

// mean and second moment for coordinates weighted by values

template <class T>
struct WeightedMoment2;

template <class T>
struct WeightedMean
    : public detail::gen_mean_acc<T, detail::accumulable_weighted_access,
                                     WeightedMean, WeightedMoment2,
                                                   WeightedMoment2>
{
    static std::string name() { return "weighted mean"; }
};
template <class T>
struct WeightedMoment2
    : public detail::gen_sum_squared_diff_acc<T,
                                            detail::accumulable_weighted_access,
                                            WeightedMoment2, WeightedMean>
{
    static std::string name() { return "weighted moment2"; }
};

} // namespace acc

namespace detail {

// mean and second moment for templated things [design study]
// - to use like temp_mean<accumulable_weighted_access>::template acc
// in place of WeightedMean
// - sample use may be accumulable_weighted_grey_access / etc.
//
template <template<class> class ACCU_ACCESS>
struct temp_sum_squared_diff;

template <template<class> class ACCU_ACCESS>
struct temp_mean
{
    template <class T>
    struct acc
        : public gen_mean_acc<T, ACCU_ACCESS, acc,
                              temp_sum_squared_diff<ACCU_ACCESS>::template acc,
                              temp_sum_squared_diff<ACCU_ACCESS>::template acc>
    {};
};

template <template<class> class ACCU_ACCESS>
struct temp_sum_squared_diff
{
    template <class T>
    struct acc
        : public gen_sum_squared_diff_acc<T, ACCU_ACCESS, acc,
                                          temp_mean<ACCU_ACCESS>::template acc>
    {};
};
//
// end of design study for mean and second moment for templated things


template <class T, template<class> class CV>
struct gen_calculator_base
    : public calculator_base<typename CV<T>::type, CV<T>::flags>
{
    typedef typename CV<T>::type type;
};

template <class T>
struct value_calculator_base
    : public gen_calculator_base<T, accumulable_value_access> {};

template <class T>
struct coord_calculator_base
    : public gen_calculator_base<T, accumulable_coord_access> {};

template <class T>
struct weighted_calculator_base
    : public gen_calculator_base<T, accumulable_weighted_access> {};


template <class T, template<class> class CV,
                   template<class> class SUM_SQUARED_DIFF>
struct gen_variance_calc
    : public gen_calculator_base<T, CV>,
      public type_lists::implies_template<T, SUM_SQUARED_DIFF>
{
    template <class ACX>
    void calculate(const ACX & x) const // calculator tuple == x.calculators
    {
        this->set_val(use<SUM_SQUARED_DIFF>::val(x)
                          / use<acc::Count>::val(x));
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Variance
    : public detail::gen_variance_calc<T, detail::accumulable_value_access,
                                          Moment2>
{
    static std::string name() { return "variance"; }
};

template <class T>
struct CoordVariance
    : public detail::gen_variance_calc<T, detail::accumulable_coord_access,
                                          CoordMoment2>
{
    static std::string name() { return "coord variance"; }
};

template <class T>
struct WeightedVariance
    : public detail::gen_variance_calc<T, detail::accumulable_weighted_access,
                                          WeightedMoment2>
{
    static std::string name() { return "weighted variance"; }
};

} // namespace acc

namespace detail {

template <class T, template<class> class CV,
                   template<class> class SUM_SQUARED_DIFF>
struct gen_unbiased_variance_calc
    : public gen_variance_calc<T, CV, SUM_SQUARED_DIFF>
{
    template <class ACX>
    void calculate(const ACX & x) const
    {
        this->set_val(use<SUM_SQUARED_DIFF>::val(x)
                         / (use<acc::Count>::val(x) - 1.0));
    }
};

} // namespace detail

namespace acc {

template <class T>
struct UnbiasedVariance
    : public detail::gen_unbiased_variance_calc<T,
                                               detail::accumulable_value_access,
                                               Moment2>
{
    static std::string name() { return "unbiased variance"; }
};
template <class T>
struct CoordUnbiasedVariance
    : public detail::gen_unbiased_variance_calc<T,
                                               detail::accumulable_coord_access,
                                               CoordMoment2>
{
    static std::string name() { return "coord unbiased variance"; }
};
template <class T>
struct WeightedUnbiasedVariance
    : public detail::gen_unbiased_variance_calc<T,
                                            detail::accumulable_weighted_access,
                                            WeightedMoment2>
{
    static std::string name() { return "weighted unbiased variance"; }
};

} // namespace acc

namespace detail {

template <class T, template<class> class CV, template<class> class SELF>
struct gen_sum_acc : public gen_accumulator_base<T, CV>
{
    /** update sum
    */
    template <class ACX>
    void updatePass1(ACX &, const T & v)
    {
        this->add_val(this->get(v));
    }
    /** update average, using weighted input.
     * <tt>stats(value, 1.0)</tt> is equivalent to the unweighted
     * call <tt>stats(value)</tt>, and <tt>stats(value, 2.0)</tt>
     * is equivalent to two unweighted calls.
     */
    template <class ACX>
    void updatePass1(ACX &, const T & v, double weight)
    {
        this->add_val(this->get(v) * weight);
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        this->set_val(use<SELF>::val(x) + use<SELF>::val(y));
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Sum
    : public detail::gen_sum_acc<T, detail::accumulable_value_access, Sum>
{
    static std::string name() { return "sum"; }
};

template <class T>
struct CoordSum
    : public detail::gen_sum_acc<T, detail::accumulable_coord_access, CoordSum>
{
    static std::string name() { return "coord sum"; }
};

template <class T>
struct WeightedSum
    : public detail::gen_sum_acc<T, detail::accumulable_weighted_access,
                                    WeightedSum>
{
    static std::string name() { return "weighted sum"; }
};

template <class T>
struct CenterOfMass
    : public detail::coord_calculator_base<T>,
      public type_lists::implies_template<T, Sum, WeightedSum>
{
    template <class ACX>
    void calculate(const ACX & x) const
    {
        using namespace detail;
        this->set_val(use<WeightedSum>::val(x) / use<Sum>::val(x));
    }
    static std::string name() { return "center of mass"; }
};

} // namespace acc

namespace detail {

template <unsigned n, class T, template<class> class CV,
                               template<class> class MEAN,
                               template<class> class SELF>
struct gen_central_moment_acc2
    : public gen_accumulator_base<T, CV>
{
    static const unsigned number_of_passes = 2;
    // reset() and unimplemented variants
    using gen_accumulator_base<T, CV>::operator();

    template <class ACX> void updatePass1(ACX &, const T &) {} // pass 1
    template <class ACX> void updatePass1(ACX &, const T &, double) {} // pass 1
    
    template <class ACX>
    void updatePassN(ACX & x, const T & v, unsigned pass) // pass 2
    {
        if (test_pass<ACX, 2>(pass))
            this->add_val(power<n>(this->get(v) - use<MEAN>::f_val(x)));
    }
};

template <class T, template<class> class CV,
                   template<class> class MEAN,
                   template<class> class SELF>
struct gen_central_moment2_acc2
    : public gen_central_moment_acc2<2, T, CV, MEAN, SELF>,
      public gen_sum_squared_diff_merge<T, SELF, MEAN>,
      public type_lists::implies_template<T, MEAN>
{
    using gen_sum_squared_diff_merge<T, SELF, MEAN>::merge;
};

template <class T, template<class> class CV,
                   template<class> class MEAN,
                   template<class> class MOMENT2,
                   template<class> class SELF>
struct gen_central_moment3_acc2
    : public gen_central_moment_acc2<3, T, CV, MEAN, SELF>,
      public type_lists::implies_template<T, MEAN, MOMENT2>
{
    typedef typename gen_central_moment3_acc2::result_type result_type;
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        double count_x = use<acc::Count>::val(x);
        double count_y = use<acc::Count>::val(y);
        double count = count_x + count_y;
        result_type delta = use<MEAN>::val(y) - use<MEAN>::val(x);
        this->set_val(  use<SELF>::val(x)
                      + use<SELF>::val(y)
                      + count_x * count_y * (count_x - count_y)
                        * power<3>(delta) / power<2>(count)
                      + 3 * (  count_x * use<MOMENT2>::val(y)
                             - count_y * use<MOMENT2>::val(x)) * delta / count);
    }
};

template <class T, template<class> class CV,
                   template<class> class MEAN,
                   template<class> class MOMENT2,
                   template<class> class MOMENT3,
                   template<class> class SELF>
struct gen_central_moment4_acc2
    : public gen_central_moment_acc2<4, T, CV, MEAN, SELF>,
      public type_lists::implies_template<T, MEAN, MOMENT2, MOMENT3>
{
    typedef typename gen_central_moment4_acc2::result_type result_type;
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        double count_x = use<acc::Count>::val(x);
        double count_y = use<acc::Count>::val(y);
        double count_x_2 = power<2>(use<acc::Count>::val(x));
        double count_y_2 = power<2>(use<acc::Count>::val(y));
        double count = count_x + count_y;
        result_type delta = use<MEAN>::val(y) - use<MEAN>::val(x);
        this->set_val(  use<SELF>::val(x)
                      + use<SELF>::val(y)
                      + count_x * count_y
                        * (count_x_2 - count_x * count_y + count_y_2)
                        * power<4>(delta) / power<3>(count)
                      + 6 * (  count_x_2 * use<MOMENT2>::val(y)
                             + count_y_2 * use<MOMENT2>::val(x))
                          * power<2>(delta) / power<2>(count)
                      + 4 * (  count_x * use<MOMENT3>::val(y)
                             - count_y * use<MOMENT3>::val(x)) * delta / count);
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Moment2_2Pass
    : public detail::gen_central_moment2_acc2<T,
                                              detail::accumulable_value_access,
                                              Mean, Moment2_2Pass>
{
    static std::string name() { return "moment2 2-pass"; }
};

template <class T>
struct Moment3
    : public detail::gen_central_moment3_acc2<T,
                                              detail::accumulable_value_access,
                                              Mean, Moment2_2Pass, Moment3>
{
    static std::string name() { return "moment3"; }
};

template <class T>
struct Moment4
    : public detail::gen_central_moment4_acc2<T,
                                              detail::accumulable_value_access,
                                              Mean, Moment2_2Pass, Moment3,
                                              Moment4>
{
    static std::string name() { return "moment4"; }
};

template <class T>
struct CoordMoment2_2Pass
    : public detail::gen_central_moment2_acc2<T,
                                              detail::accumulable_coord_access,
                                              CoordMean, CoordMoment2_2Pass>
{
    static std::string name() { return "coord moment2 2-pass"; }
};

template <class T>
struct CoordMoment3
    : public detail::gen_central_moment3_acc2<T,
                                              detail::accumulable_coord_access,
                                              CoordMean, CoordMoment2_2Pass,
                                              CoordMoment3>
{
    static std::string name() { return "coord moment3"; }
};

template <class T>
struct CoordMoment4
    : public detail::gen_central_moment4_acc2<T,
                                              detail::accumulable_coord_access,
                                              CoordMean, CoordMoment2_2Pass,
                                              CoordMoment3, CoordMoment4>
{
    static std::string name() { return "coord moment4"; }
};

template <class T>
struct WeightedMoment2_2Pass
    : public detail::gen_central_moment2_acc2<T,
                                              detail::accumulable_weighted_access,
                                              WeightedMean,
                                              WeightedMoment2_2Pass>
{
    static std::string name() { return "weighted moment2 2-pass"; }
};

template <class T>
struct WeightedMoment3
    : public
        detail::gen_central_moment3_acc2<T, detail::accumulable_weighted_access,
                                         WeightedMean, WeightedMoment2_2Pass,
                                         WeightedMoment3>
{
    static std::string name() { return "weighted moment3"; }
};
                                         
template <class T>
struct WeightedMoment4
    : public
        detail::gen_central_moment4_acc2<T, detail::accumulable_weighted_access,
                                         WeightedMean, WeightedMoment2_2Pass,
                                         WeightedMoment3, WeightedMoment4>
{
    static std::string name() { return "weighted moment4"; }
};

} // namespace acc

namespace detail {

template <class T, template<class> class CV,
                   template<class> class M2_ACC2,
                   template<class> class M3_ACC2>
struct gen_skewness_calc2
    : public gen_calculator_base<T, CV>,
      public type_lists::implies_template<T, M2_ACC2, M3_ACC2>
{
    template <class ACX>
    void calculate(const ACX & x) const
    {
        using std::sqrt;
        this->set_val(use<M3_ACC2>::val(x)
                      / sqrt(power<3>(use<M2_ACC2>::val(x))));
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Skewness
    : public detail::gen_skewness_calc2<T, detail::accumulable_value_access,
                                           Moment2_2Pass, Moment3>
{
    static std::string name() { return "skewness"; }
};

template <class T>
struct CoordSkewness
    : public detail::gen_skewness_calc2<T, detail::accumulable_coord_access,
                                           CoordMoment2_2Pass, CoordMoment3>
{
    static std::string name() { return "coord skewness"; }
};

template <class T>
struct WeightedSkewness
    : public detail::gen_skewness_calc2<T, detail::accumulable_weighted_access,
                                           WeightedMoment2_2Pass,
                                           WeightedMoment3>
{
    static std::string name() { return "weighted skewness"; }
};

} // namespace acc

namespace detail {

template <class T, template<class> class CV,
                   template<class> class M2_ACC2,
                   template<class> class M4_ACC2>
struct gen_kurtosis_calc2
    : public gen_calculator_base<T, CV>,
      public type_lists::implies_template<T, M2_ACC2, M4_ACC2>
{
    template <class ACX>
    void calculate(const ACX & x) const
    {
        this->set_val(use<M4_ACC2>::val(x) / power<2>(use<M2_ACC2>::val(x)));
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Kurtosis
    : public detail::gen_kurtosis_calc2<T, detail::accumulable_value_access,
                                           Moment2_2Pass, Moment4>
{
    static std::string name() { return "kurtosis"; }
};
                                           
template <class T>
struct CoordKurtosis
    : public detail::gen_kurtosis_calc2<T, detail::accumulable_coord_access,
                                           CoordMoment2_2Pass, CoordMoment4>
{
    static std::string name() { return "coord kurtosis"; }
};
                                           
template <class T>
struct WeightedKurtosis
    : public detail::gen_kurtosis_calc2<T, detail::accumulable_weighted_access,
                                           WeightedMoment2_2Pass,
                                           WeightedMoment4>
{
    static std::string name() { return "weighted kurtosis"; }
};


} // namespace acc

namespace detail {

template<class X> // use the dereferenced (update) type for accumulation
struct plain_collector_load
{
    typedef X type;
};

template<class X> // use an explicitly specified type for accumulation
struct collector_load_is
{
    template<class>
    struct temp
    {
        typedef X type;
    };
};

template <class T, template<class> class CV,
                   template<class> class LOAD = plain_collector_load>
struct gen_collector_base
    : public collector_base<T, typename
                               LOAD<typename CV<T>::type>::type, CV<T>::flags>
{
    typedef typename CV<T>::type type;
    type get(const T & v) { return CV<T>::get(v); }
};

template <class T, template<class> class LOAD = plain_collector_load>
struct value_collector_base
    : public gen_collector_base<T, accumulable_value_access, LOAD> {};

template <class T, template<class> class LOAD = plain_collector_load>
struct coord_collector_base
    : public gen_collector_base<T, accumulable_coord_access, LOAD> {};

template <class T, template<class> class LOAD = plain_collector_load>
struct weighted_collector_base
    : public gen_collector_base<T, accumulable_weighted_access, LOAD> {};

template <class T, class V>
struct counting_value_collector_base
    : public value_collector_base<T, collector_load_is<V>::template temp>
{};
template <class T, class V>
struct counting_coord_collector_base
    : public coord_collector_base<T, collector_load_is<V>::template temp>
{};
template <class T, class V>
struct counting_weighted_collector_base
    : public weighted_collector_base<T, collector_load_is<V>::template temp>
{};

// no features_double_promote for min and max.
template <class T, template<class> class CV>
struct gen_plain_collector_base
    : public gen_collector_base<T, CV> {};


template <class T, template<class> class CV, template<class> class SELF>
struct gen_min_acc : public gen_plain_collector_base<T, CV>
{
    typedef typename gen_min_acc::type type;
    void reset()
    {
        this->set_val(NumericTraits<type>::max());
    }
    gen_min_acc()
    {
        reset();
    }
    /** update
    */
    template <class ACX>
    void updatePass1(ACX &, const T & v, double = 0)
    {
        updateMin(this->value, this->get(v));
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        this->set_val(use<SELF>::val(x));
        updateMin(this->value, use<SELF>::val(y));
    }
};

template <class T, template<class> class CV, template<class> class SELF>
struct gen_max_acc : public gen_plain_collector_base<T, CV>
{
    typedef typename gen_max_acc::type type;
    void reset()
    {
        this->set_val(NumericTraits<type>::min());
    }
    gen_max_acc()
    {
        reset();
    }
    /** update
    */
    template <class ACX>
    void updatePass1(ACX &, const T & v, double = 0)
    {
        updateMax(this->value, this->get(v));
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        this->set_val(use<SELF>::val(x));
        updateMax(this->value, use<SELF>::val(y));
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Min
    : public detail::gen_min_acc<T, detail::accumulable_value_access, Min>
{
    static std::string name() { return "min"; }
};
template <class T>
struct CoordMin
    : public detail::gen_min_acc<T, detail::accumulable_coord_access, CoordMin>
{
    static std::string name() { return "coord min"; }
};

template <class T>
struct WeightedMin
    : public detail::gen_min_acc<T, detail::accumulable_weighted_access,
                                    WeightedMin>
{
    static std::string name() { return "weighted min"; }
};

template <class T>
struct Max
    : public detail::gen_max_acc<T, detail::accumulable_value_access, Max>
{
    static std::string name() { return "max"; }
};
template <class T>
struct CoordMax
    : public detail::gen_max_acc<T, detail::accumulable_coord_access, CoordMax>
{
    static std::string name() { return "coord max"; }
};
template <class T>
struct WeightedMax
    : public detail::gen_max_acc<T, detail::accumulable_weighted_access,
                                    WeightedMax>
{
    static std::string name() { return "weighted max"; }
};

template <class T>
struct OctetHistogram
    : public detail::counting_value_collector_base<T, TinyVector<double, 256> >
{
    typedef typename OctetHistogram::type type;
    /** update
    */
    template <class ACX>
    void updatePass1(ACX &, const T & v, double weight = 1)
    {
        type pos = this->get(v);
        if (pos > 255)
            pos = 255;
        if (pos < 0)
            pos = 0;
        this->value[static_cast<unsigned>(pos)] += weight;
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        using namespace detail;
        this->set_val(  use_d<OctetHistogram>::val(x)
                      + use_d<OctetHistogram>::val(y));
    }
    static std::string name() { return "octet histogram"; }
};

} // namespace acc

// -- end feature template classes --

namespace detail {

struct BinderReset
{
    template <class F>
    void operator()(F & first) // via exec_bound
    {
        first.reset();
    }
};
template <class A, class B>
struct BinderPass1
{
    A & v;
    const B & w;
    BinderPass1(A & v_, const B & w_) : v(v_), w(w_) {}
    template <class F>
    void operator()(F & first) // via exec_bound
    {
        first.updatePass1(v, w);
    }
};
template <class A, class B, class C>
struct BinderPassN
{
    A & v;
    const B & w;
    const C x;
    BinderPassN(A & v_, const B & w_, C x_) : v(v_), w(w_), x(x_) {}
    template <class F>
    void operator()(F & first) // via exec_bound
    {
        first.updatePass1(v, w, x); // Pass1, not N.
    }
    template <class F>
    void call(F & first) // via call_bound
    {
        first.updatePassN(v, w, x);
    }
};
template <class A>
struct BinderMerge
{
    const A & v;
    const A & w;
    const A & x;
    BinderMerge(const A & v_, const A & w_, const A & x_) : v(v_), w(w_), x(x_)
    {}
    template <class F>
    void operator()(F & first) // via exec_bound
    {
        first.merge(v, w, x);
    }
};

// glue code for virtual function calls
// from type_lists::cond_virtual_tuple
struct v_feature_dispatch
{
    template <class ACX, class T> 
    struct member_base_type
    {
        virtual void reset() = 0;
        virtual void updatePass1(ACX &, const T &) = 0;
        virtual void updatePass1(ACX &, const T &, double) = 0;
        virtual void updatePassN(ACX &, const T &, unsigned) = 0;
        virtual void merge(const ACX &, const ACX &, const ACX &) = 0;
        virtual ~member_base_type() {}
    };
    template <class ACX, class T, class V>
    struct load_type : public member_base_type<ACX, T>, public V
    {
        void reset()
        {
            V::reset();
        }
        void updatePass1(ACX & x, const T & v)
        {
            V::updatePass1(x, v);
        }
        void updatePass1(ACX & x, const T & v, double w)
        {
            V::updatePass1(x, v, w);
        }
        void updatePassN(ACX & x, const T & v, unsigned n)
        {
            V::updatePassN(x, v, n);
        }
        void merge(const ACX & z, const ACX & x, const ACX & y)
        {
            V::merge(z, x, y);
        }
    };
};

template <class T, class L, class EXTRACTOR = type_lists::nil,
          template<class, template<class> class TES>
          class TUPLE = type_lists::tuple>
struct feature_context_base
    : public TUPLE<L, passes<2>::template at_least>
{
    typedef EXTRACTOR            extractor_type;
    typedef feature_context_base feature_context_base_type;

    template <class V>
    const V & member() const
    {
        return static_cast<const V &>(this->template ref<V>());
    }
    template <class V>
    V & member()
    {
        return static_cast<V &>(this->template ref<V>());
    }
    template <template<class> class V>
    const V<T> & member_template() const
    {
        return member< V<T> >();
    }
    template <template<class> class V>
    V<T> & member_template()
    {
        return member< V<T> >();
    }

    template <template<class> class V>
    struct get
    {
        typedef V<T> type;
        typedef typename type::result_type result_type;
    };

    template <class V>
    const typename V::result_type & value() const
    {
        return member<V>().get_value();
    }
    template <template<class> class V>
    const typename V<T>::result_type & value_template() const
    {
        return value< V<T> >();
    }

    template <class V>
    const typename V::result_type & soft_value()
    {
        return member<V>().soft_value();
    }
    template <template<class> class V>
    const typename V<T>::result_type & soft_value_template()
    {
        return soft_value< V<T> >();
    }
    template <class V>
    const typename V::result_type & soft_value() const
    {
        return member<V>().soft_value();
    }
    template <template<class> class V>
    const typename V<T>::result_type & soft_value_template() const
    {
        return soft_value< V<T> >();
    }

    // default implementations to be used by accu_context_base:
    template <class, class A> void calculate(A &) const {}
    template <template<class> class, class A> void calculate_template(A &) const
    {}
    template <class V, class ACX>
    void set_accs(ACX &)
    {
        this->template set<V>();
    }
    template <class V, class ACX>
    bool is_set_accs(const ACX &) const
    {
        return this->template is_set<V>();
    }
};

// calculators -- using just use an ordinary static type_lists::tuple, i.e.,
// they are always available, although not necessarily with meaningful results.
template <class T, class L, class EXTRACTOR = type_lists::nil>
struct calc_context : public feature_context_base<T, L, EXTRACTOR>
{
    // by-need calculation via get::val
    template <class V, class ACX>
    void calculate(ACX & x) const
    {
        this->template ref<V>().calculate(x);
    }
    template <template<class> class V, class ACX>
    void calculate_template(ACX & x) const
    {
        this->template calculate< V<T> >(x);
    }
    // set_accs(): set all accumulators that are implied by calculator class V
    template <class V>
    struct set_accs_op
    {
        template <class ACX>
        static void exec(ACX & x)
        {
            x.template set_accs<V>(x);
        }
    };
    template <class V, class ACX>
    void set_accs(ACX & x) const
    {
        type_lists::for_all_implied<V, set_accs_op>(x);
    }
    // consider a calculator as set if all implied features are (considered to
    // be) set:
    template <class V>
    struct is_set_accs_op
    {
        template <class PAIR>
        static void exec(PAIR & z)
        {
            *z.first = *z.first && z.second->template is_set<V>();
        }
    };
    template <class V, class ACX>
    bool is_set_accs(const ACX & x) const
    {
        bool ret = true; // OK if there are no implied features
        std::pair<bool*, const ACX*> z(&ret, &x); // std::pair vetoes references
        type_lists::for_all_implied<V, is_set_accs_op>(z);
        return ret;
    }
};

template <class L>
struct feature_types
{
    template <class X> struct one_pass : public passes<1>::at_least<X> {};
    template <class X> struct two_pass : public passes<2>::at_least<X> {};
    typedef typename type_lists::topological_sort_expanded<L>::type  list;
    typedef typename type_lists::remove_if_not<one_pass, list>::type accu_list;
    typedef typename type_lists::remove_if    <one_pass, list>::type calc_list;

    typedef typename type_lists::remove_if_not<two_pass, list>::type pass_list;
    typedef typename type_lists::remove_if<one_pass, typename type_lists::
                                          implies_expand<pass_list>::type>::type
        sync_list; // calculators that are used by passes >= 2
    static const unsigned max_passes
        = type_lists::max_value<passes_number, unsigned, accu_list>::value;
    // base class for extractors:
    typedef typename type_lists::map_unique<get_result_type, list>::type
        result_type_list;
    typedef extractor_base<result_type_list> extractor_type;
};

template <class T, class L, template<class, template<class> class TEST>
                            class TUPLE = type_lists::tuple>
struct accu_context_base
    : public feature_context_base<T, typename feature_types<L>::accu_list,
                                     typename feature_types<L>::extractor_type,
                                  TUPLE>,
      public UnaryAnalyserTag
{
    typedef          accu_context_base                            feature_type;
    typedef typename accu_context_base::tuple_type                base_type;
    typedef typename accu_context_base::list_type                 list_type;
    typedef typename accu_context_base::feature_context_base_type feature_cx;
    typedef          accu_context_base                            extra_passes;
    typedef          feature_types<L>                             feature_lists;
    static const unsigned max_passes = feature_lists::max_passes;

    // for Python interface and similar:
    typedef typename accu_context_base::extractor_type extractor_type;

    typedef typename feature_lists::list                          full_list;
    typedef typename feature_lists::accu_list                     accu_list;
    typedef typename feature_lists::calc_list                     calc_list;
    // features that do not require even a single pass:
    typedef typename feature_lists::sync_list                     sync_list;
    typedef calc_context<T, calc_list, extractor_type>            calc_type;

    calc_type calculators;
    unsigned  number_of_passes;
    accu_context_base() : number_of_passes(0) { update_number_of_passes(); }

    static const unsigned accu_size = type_lists::size<accu_list>::of;
    static const unsigned calc_size = type_lists::size<calc_list>::of;
    static const unsigned size = accu_size + calc_size;

    template <template<class> class V> // works around cl.exe bugs
    struct has_feature
    {
        static const bool value
            =  base_type::template has_element< V<T> >::value
            || calc_type::template has_element< V<T> >::value;
    };

    // use type feature_cx here in place of accu_context_base to
    // prevent infinite recursion.
    template <class V>
    struct chooser : public type_lists::choose_tuple<V, list_type,
                                                     feature_cx, calc_type>
    {};
    template <class V>
    typename chooser<V>::ref_type ref()
    {
        return chooser<V>::ref(*this, calculators);
    }
    template <class V>
    typename chooser<V>::const_ref_type ref() const
    {
        return chooser<V>::ref(*this, calculators);
    }
 
    // dispatch sub-accumulator member requests to calculators and collectors:
    template <class V>
    const V & member() const
    {
        return static_cast<const V &>(ref<V>());
    }
    template <class V>
    V & member()
    {
        return static_cast<V &>(ref<V>());
    }
    template <template<class> class V>
    const V<T> & member_template() const
    {
        return member< V<T> >();
    }
    template <template<class> class V>
    V<T> & member_template()
    {
        return member< V<T> >();
    }

    // dispatch value requests to calculators and collectors:
    template <class V>
    const typename V::result_type & value() const
    {
        return chooser<V>::at(*this, calculators).template value<V>();
    }
    template <template<class> class V>
    const typename V<T>::result_type & value_template() const
    {
        return this->template value< V<T> >();
    }
    template <class V>
    const typename V::result_type & soft_value()
    {
        return chooser<V>::at(*this, calculators).template soft_value<V>();
    }
    template <template<class> class V>
    const typename V<T>::result_type & soft_value_template()
    {
        return soft_value< V<T> >();
    }
    template <class V>
    const typename V::result_type & soft_value() const
    {
        return chooser<V>::at(*this, calculators).template soft_value<V>();
    }
    template <template<class> class V>
    const typename V<T>::result_type & soft_value_template() const
    {
        return soft_value< V<T> >();
    }

    // work around the workaround of using templated inheritance as a
    // surrogate for missing template aliases in C++ 2003,
    // such as used for the definition of SelectedAccumulators,
    // for which the following inline conversion operator is (silently)
    // called within Binder*
    typedef typename accu_context_base::derived_type derived_type;
    operator derived_type & ()
    {
        return static_cast_2<derived_type, base_type>(*this);
    }
    operator const derived_type & () const
    {
        return static_cast_2<const derived_type, const base_type>(*this);
    }

    template <class V>
    struct number_of_passes_op
    {
        template <class ACX>
        static void exec(ACX & x)
        {
            const unsigned number_of_passes = V::number_of_passes;
            if (x.template is_set<V>())
               updateMax(x.number_of_passes, number_of_passes);
        }
    };
    void update_number_of_passes()
    {
        type_lists::for_all<accu_list, number_of_passes_op>(*this);
    }

    template <class V>
    void set()
    {
        chooser<V>::at(*this, calculators).template set_accs<V>(*this);
        // the number of passes may have been increased:
        update_number_of_passes();
    }
    template <template<class> class V>
    void set_template()
    {
        set< V<T> >();
    }

    template <class V>
    bool is_set() const
    {
        return chooser<V>::at(*this, calculators).template
                                                  is_set_accs<V>(*this);
    }
    template <template<class> class V>
    bool is_set_template() const
    {
        return is_set< V<T> >();
    }

    template <class V, class ACX>
    void calculate(ACX & x) const
    {
        chooser<V>::at(*this, calculators).template calculate<V>(x);
    }
    template <template<class> class V, class ACX>
    void calculate_template(ACX & x) const
    {
        chooser< V<T> >::at(*this, calculators).template
                                                calculate_template<V>(x);
    }

    template <class V>
    struct calc_sync_op
    {
        static void exec(const accu_context_base & x)
        {
            x.calculators.template calculate<V>(x);
        }
    };
    template <class V>
    struct calc_sync_op<V*> : public calc_sync_op<V> {};
    void calc_sync() // run all calculators that are used by passes >= 2
    {
        type_lists::for_all<sync_list, calc_sync_op>(*this);
    }

        /** re-init
        */
    void reset()
    {
        BinderReset binder;
        base_type::exec_bound(binder);
    }
    /** update
    */
    void operator()(const T & v)
    {
        BinderPass1<accu_context_base, T> binder(*this, v);
        base_type::exec_bound(binder);
    }
    /** update using weighted input -> not yet completed.
    */
    void operator()(const T & v, double weight)
    {
        BinderPassN<accu_context_base, T, double> binder(*this, v, weight);
        base_type::exec_bound(binder);
    }

    // update: passes >= 2
    void updatePassN(const T & v, unsigned pass_number)
    {
        BinderPassN<accu_context_base, T, unsigned>
            binder(*this, v, pass_number);
        base_type::call_bound(binder);
    }
    struct pass_n_dispatch
    {
        accu_context_base & accu;
        unsigned            pass_number;
        pass_n_dispatch(accu_context_base & a, unsigned n)
            : accu(a), pass_number(n) {}
        void operator()(const T & v)
        {
            accu.updatePassN(v, pass_number);
        }
    };
    pass_n_dispatch pass_n(unsigned n)
    {
        if (n < 2 || n > max_passes)
            vigra_fail("accu_context_base::pass_n(): inconsistent use.");
        return pass_n_dispatch(*this, n);
    }

    unsigned numberOfPasses()
    {
        return number_of_passes;
    }

    /** select accumulators that are set on the operand
    */
    void selectSelected(const accu_context_base & x)
    {
        // merge set features (if any):
        x.transfer_set_to(*this);
    }

    /** merge two statistics: *this <- (x, y)
    */
    void operator()(const accu_context_base & x0, const accu_context_base & y0)
    {
        // first, merge the selection of accumulators (if any):
        selectSelected(x0);
        selectSelected(y0);
        copy_if_same_as<accu_context_base> x(x0, *this);
        copy_if_same_as<accu_context_base> y(y0, *this);

        BinderMerge<accu_context_base> binder(*this, x(), y());
        base_type::exec_bound(binder);
    }
    /** merge two statistics: *this <- (this, x)
    */
    void operator()(const accu_context_base & x)
    {
        operator()(*this, x);
    }

    /** names of all accumulators and calculators
    */
//     ... names()()
//     {
//         ...;
//     }

    // compatibility with ArrayOfRegionStatistics
    typedef T argument_type;
    typedef void result_type;
    result_type operator()() const {}
};

} // namespace detail

// The feature accessing functions get(), member(), select(), and isSelected()
// are implemented as free template functions, since an implementation as member
// functions would necessitate the use of the 'template' disambiguator in cases
// where they are used within template classes or template functions.

// enable_if to check for availability of features
template <template<class> class V,
          class T, class L, template<class, template<class> class> class P>
inline typename
detail::accu_context_base<T, L, P>::template get<V>::result_type
get(const detail::accu_context_base<T, L, P> & x, typename
    enable_if<detail::accu_context_base<T, L, P>::template
                                                  has_feature<V>::value>::type*
                                                                            = 0)
{
    return detail::use<V>::val(x);
}

template <template<class> class V,
          class T, class L, template<class, template<class> class> class P>
inline void
select(detail::accu_context_base<T, L, P> & x, typename
       enable_if<detail::accu_context_base<T, L, P>::template
                                                   has_feature<V>::value>::type*
                                                                            = 0)
{
    detail::use<V>::set(x);
}
template <template<class> class V,
          class T, class L, template<class, template<class> class> class P>
inline bool
isSelected(detail::accu_context_base<T, L, P> & x, typename
       enable_if<detail::accu_context_base<T, L, P>::template
                                                   has_feature<V>::value>::type*
                                                                            = 0)
{
    return detail::use<V>::is_set(x); // always returns true for calculators
}

template <template<class> class V,
          class T, class L, template<class, template<class> class> class P>
inline
const typename detail::accu_context_base<T, L, P>::template get<V>::type &
member(const   detail::accu_context_base<T, L, P> & x, typename
       enable_if<detail::accu_context_base<T, L, P>::template
                                                   has_feature<V>::value>::type*
                                                                            = 0)
{
    return detail::use<V>::member(x);
}

template <template<class> class V,
          class T, class L, template<class, template<class> class> class P>
inline
typename      detail::accu_context_base<T, L, P>::template get<V>::type &
member(       detail::accu_context_base<T, L, P> & x, typename
    enable_if<detail::accu_context_base<T, L, P>::template
                                                  has_feature<V>::value>::type*
                                                                            = 0)
{
    return detail::use<V>::member(x);
}

// 'type_lists::use_template_list' requires this
// (binding of the TUPLE parameter):

template <class T, class L>
struct Accumulators : public detail::accu_context_base<T, L> {};

template <class T, class L>
struct NestedAccumulators
    : public detail::accu_context_base<T, L, type_lists::bit_cond_tuple> {};

template <class T, class L>
struct BitAccumulators
    : public detail::accu_context_base<T, L, type_lists::bit_tuple> {};

template <class T, class L>
struct PlainCondAccumulators
    : public detail::accu_context_base<T, L, type_lists::cond_tuple_plain> {};

template <class T, class L>
struct CondAccumulators
    : public detail::accu_context_base<T, L, type_lists::cond_tuple> {};

template <class T, class L>
struct SelectedAccumulators
    : public detail::accu_context_base<T, L,
        type_lists::cond_virtual_tuple<SelectedAccumulators, T,
                                       detail::v_feature_dispatch>::template
                                                                    type>
{};

namespace detail {

// Some member functions may take the accumulator object as a parameter despite
// it being _currently_ unused.
template <class ACX>
struct feature_selector_base
{
    virtual void select(ACX &) const = 0;
    virtual bool isSelected(const ACX &) const = 0;
    virtual std::string name() const = 0;
    virtual void extract(const ACX &, typename ACX::extractor_type &) const = 0;
    virtual unsigned numberOfPasses(const ACX &) const = 0;
    virtual feature_selector_base* clone() const = 0;
};

// Class feature_selector is the dispatch class to use for bridging runtime
// access to a number of features (such as within the Python interface) to
// the compile-time selection of get<>() / select<>() / etc. via C++ symbol
// names.
// Class feature_selector has got value semantics and is thus safe to use in,
// e.g., arrays, or the containers of the standard library.
// Use get_feature_selectors() to create feature_selector objects corresponding
// to a certain feature class.
template <class ACX>
class feature_selector
{
private:
    typedef feature_selector_base<ACX>* ptr_type;
    ptr_type p;
public:
    feature_selector() : p(0) {}
    feature_selector(ptr_type p_) : p(p_) {}
private:
    void check_pointer() const
    {
        if (!p)
            vigra_fail("class detail::feature_selector: empty object.");
    }
public:
    void select(ACX & x) const
    {
        check_pointer();
        p->select(x);
    }
    bool isSelected(const ACX & x) const
    {
        check_pointer();
        return p->isSelected(x);
    }
    std::string name() const
    {
        check_pointer();
        return p->name();
    }
    void extract(const ACX & x, typename ACX::extractor_type & g) const
    {
        check_pointer();
        p->extract(x, g);
    }
    unsigned numberOfPasses(const ACX & x) const
    {
        check_pointer();
        return p->numberOfPasses(x);
    }

    feature_selector & operator=(const feature_selector & e)
    {
        ptr_type tmp = e.p->clone();
        delete p;
        p = tmp;
        return *this;
    }
    feature_selector(const feature_selector & e)
        : p(0)
    {
        if (e.p) // properly care for empty original
            p = e.p->clone();
    }
    ~feature_selector()
    {
        delete p;
    }
};

// These classes map run-time access via the feature_selector wrapper
// class to class-name based compile-time access to features.
template <class ACX, class V>
struct feature_selector_impl : public feature_selector_base<ACX>
{
    void select(ACX & x) const
    {
        x.template set<V>();
    }
    bool isSelected(const ACX & x) const
    {
        return x.template is_set<V>();
    }
    std::string name() const
    {
        return V::name();
    }
    void extract(const ACX & x, typename ACX::extractor_type & g) const
    {
        x.template member<V>().calculate(x);
        x.template member<V>().extract(g);
    }
    unsigned numberOfPasses(const ACX & x) const
    {
        return V::number_of_passes;
    }
    feature_selector_impl* clone() const
    {
        return new feature_selector_impl(*this);
    }
};

template <class ACX>
struct feature_selector_op
{
    template <class V>
    struct get
    {
        template <class OutputIterator>
        static void exec(OutputIterator & j)
        {
            *j = detail::feature_selector<ACX>(
                                     new detail::feature_selector_impl<ACX, V>);
            ++j;
        }
    };
};
// fill a sequence pointed to by an output iterator
// with feature_selector objects:
template <class ACX, class OutputIterator>
void get_feature_selectors(OutputIterator j)
{
    type_lists::for_all<typename ACX::full_list,
                        feature_selector_op<ACX>::template get
                       >(j);
}

// N_STRING must be at least explicitly constructible from std::string,
// and the std::map-like MAP template class must implement two additional member
// functions, the advisory 'void reserve(unsigned)' as well as
// std::string name_convert(const N_STRING &) const
template <class ACX, class N_STRING, template<class, class> class MAP>
struct selector_name_dispatch
{
    typedef feature_selector<ACX>             selector_type;
    typedef MAP<N_STRING, selector_type>      map_type;
    typedef typename map_type::const_iterator const_iterator;
    static const unsigned size = ACX::size;
    map_type selectors;
    selector_name_dispatch()
    {
        selector_type input[size];
        get_feature_selectors<ACX>(input);
        selectors.reserve(size);
        for (unsigned i = 0; i != size; ++i)
            selectors[N_STRING(input[i].name())] = input[i];
    }
    template <class OutputIterator>
    void names(OutputIterator j) const
    {
        for (const_iterator i = selectors.begin(); i != selectors.end();
                                                                       ++i, ++j)
            *j = (*i).first;
    }
    const selector_type & ref(const N_STRING & name) const
    {
        const_iterator i = selectors.find(name);
        if (i == selectors.end())
            vigra_fail("A feature named '" + selectors.name_convert(name)
                       + "' was not found.");
        return (*i).second;
    }
    void select(const N_STRING & n, ACX & x) const
    {
        ref(n).select(x);
    }
    bool isSelected(const N_STRING & n, const ACX & x) const
    {
        return ref(n).isSelected(x);
    }
    void extract(const N_STRING & n, const ACX & x, typename
                                     ACX::extractor_type & g) const
    {
        ref(n).extract(x, g);
    }
    unsigned numberOfPasses(const N_STRING & n, const ACX & x) const
    {
        return ref(n).numberOfPasses(x);
    }
};

} // namespace detail

// mask cl.exe shortcomings [end]
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

} // namespace vigra

#endif // VIGRA_OBJECT_FEATURES_HXX
