#ifndef VIGRA_OBJECT_FEATURES_HXX
#define VIGRA_OBJECT_FEATURES_HXX

#include <cmath>

#include "error.hxx"
#include "numerictraits.hxx"
#include "functortraits.hxx"
#include "type_lists.hxx"
#include "algorithm.hxx"
#include "mathutil.hxx"

namespace vigra {

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

template <class T>
struct features_real_promote
{
    typedef typename NumericTraits<T>::RealPromote TraitsRealPromote;
    typedef typename
        IfBool<
            IsSameType<
                Error_NumericTraits_not_specialized_for_this_case,
                TraitsRealPromote>::value, typename
            get_optional_real_promote<T>::type,
            TraitsRealPromote>::type
        type;
};

template <class V>
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
        /** (re-)init
        */
    void reset()
    {
        this->value = result_type(); // guaranteed to be zero for built-in types
    }
    feature_base()
    {
        reset();
    }
    feature_base(const V & v) : value(v) {}
};

// all collectors / accumulators must implement void reset()

template <class T, class V>
struct collector_base : public feature_base<V>
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

template <class T, class A = T>
struct numeric_accumulator_base
    : public collector_base<T, typename features_real_promote<A>::type>
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

template <class T>
struct calculator_base
    : public feature_base<typename features_real_promote<T>::type>
{
    // the actual calculator class templates should implement
    // template <class ACX> void calculate(const ACX & x),
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
        // return x.template value_template<V>(); // works with gcc 4.4 -- bug?
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
struct accumulable_coord_access
{
    typedef T accumulable_type;
    typedef T type;
    static const T & get(const T & v) { return v; }
};
template<class T>
struct accumulable_value_access
{
    typedef T accumulable_type;
    typedef T type;
    static const T & get(const T & v) { return v; }
};
template<class T>
struct accumulable_weighted_access
{
    typedef T accumulable_type;
    typedef T type;
    static const T & get(const T & v) { return v; }
};

} // namespace detail

// -- start feature template classes --

namespace acc {

using namespace detail;

template <class T>
struct Count : public collector_base<T, double>
{
    /** update
     */
    template <class ACX>
    void updatePass1(ACX &, const T &)
    {
        ++this->value;
    }
    /** update using weighted input
     */
    template <class ACX>
    void updatePass1(ACX & x, const T &, double weight)
    {
        this->value += weight;
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        this->value = use_d<Count>::val(x) + use_d<Count>::val(y);
    }
};

} // namespace acc

namespace detail {

template <class T, template<class> class CV,
                   template<class> class SELF,
                   template<class> class SUM_SQUARED_DIFF,
                   template<class> class COVAR = SUM_SQUARED_DIFF>
struct gen_mean_acc : public numeric_accumulator_base<T, typename CV<T>::type>,
                      public type_lists::uses_template<T, acc::Count>
                                       ::template follows<SUM_SQUARED_DIFF,
                                                          COVAR>
{
    /** update average
    */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v)
    {
        this->value
            += (CV<T>::get(v) - this->value)
               / use<acc::Count>::f_val(accumulators);
    }
    /** update average, using weighted input.
     * <tt>stats(value, 1.0)</tt> is equivalent to the unweighted
     * call <tt>stats(value)</tt>, and <tt>stats(value, 2.0)</tt>
     * is equivalent to two unweighted calls.
     */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v, double weight)
    {
        this->value += (CV<T>::get(v) - this->value) * weight
                       / use<acc::Count>::f_val(accumulators);
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        this->value = (   use<acc::Count>::val(x) * use<SELF>::val(x)
                        + use<acc::Count>::val(y) * use<SELF>::val(y)
                      ) / use<acc::Count>::val(z);
    }
};

template <class T, template<class> class CV,
                   template<class> class SELF,
                   template<class> class MEAN>
struct gen_sum_squared_diff_acc
    : public numeric_accumulator_base<T, typename CV<T>::type>,
      public type_lists::uses_template<T, acc::Count>
                       ::template implies<MEAN>
{
    typedef typename gen_sum_squared_diff_acc::result_type result_type;

    /** update sum of squared differences
    */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v)
    {
        result_type t1 = CV<T>::get(v) - use<MEAN>::f_val(accumulators);
        result_type t2 = t1 / use<acc::Count>::f_val(accumulators);
        this->value += t1 * t2 * (use<acc::Count>::f_val(accumulators) - 1.0);
    }
    /** update sum of squared differences, using weighted input.
     * <tt>stats(value, 1.0)</tt> is equivalent to the unweighted
     * call <tt>stats(value)</tt>, and <tt>stats(value, 2.0)</tt>
     * is equivalent to two unweighted calls.
     */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v, double weight)
    {
        result_type t1 = CV<T>::get(v) - use<MEAN>::f_val(accumulators);
        result_type t2 = t1 * weight / use<acc::Count>::f_val(accumulators);
        // count > weight is ensured by preceding acc::Count::operator() call
        this->value
            += t1 * t2 * (use<acc::Count>::f_val(accumulators) - weight);
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        double count_x = use<acc::Count>::val(x);
        double count_y = use<acc::Count>::val(y);
        this->value =   use<SELF>::val(x)
                      + use<SELF>::val(y)
                      + count_x * count_y / (count_x + count_y)
                        * (use<MEAN>::val(x) - use<MEAN>::val(y))
                        * (use<MEAN>::val(x) - use<MEAN>::val(y));
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
                                                           CV<T>::type>::type>,
      public type_lists::uses_template<T, acc::Count>
                       ::template implies<MEAN>
{
    typedef typename gen_covariance_acc::result_type result_type;

    /** update sum of squared differences
    */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v)
    {
//         result_type t1 = CV<T>::get(v) - use<MEAN>::f_val(accumulators);
//         result_type t2 = t1 / use<acc::Count>::f_val(accumulators);
//         this->value += t1 * t2 * (use<acc::Count>::f_val(accumulators) - 1.0);
    }
    /** update sum of squared differences, using weighted input.
     * <tt>stats(value, 1.0)</tt> is equivalent to the unweighted
     * call <tt>stats(value)</tt>, and <tt>stats(value, 2.0)</tt>
     * is equivalent to two unweighted calls.
     */
    template <class ACX>
    void updatePass1(ACX & accumulators, const T & v, double weight)
    {
//         result_type t1 = CV<T>::get(v) - use<MEAN>::f_val(accumulators);
//         result_type t2 = t1 * weight / use<acc::Count>::f_val(accumulators);
//         // count > weight is ensured by preceding acc::Count::operator() call
//         this->value
//             += t1 * t2 * (use<acc::Count>::f_val(accumulators) - weight);
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
                                     Moment2, Moment2> {};
template <class T>
struct Moment2
    : public detail::gen_sum_squared_diff_acc<T,
                                              detail::accumulable_value_access,
                                              Moment2, Mean
                                             > {};

// mean and second moment for coordinates

template <class T>
struct CoordMoment2;

template <class T>
struct CoordMean
    : public detail::gen_mean_acc<T, detail::accumulable_coord_access, 
                                     CoordMean, CoordMoment2, CoordMoment2> {};
template <class T>
struct CoordMoment2
    : public detail::gen_sum_squared_diff_acc<T,
                                              detail::accumulable_coord_access,
                                              CoordMoment2, CoordMean
                                             > {};

// mean and second moment for coordinates weighted by values

template <class T>
struct WeightedMoment2;

template <class T>
struct WeightedMean
    : public detail::gen_mean_acc<T, detail::accumulable_weighted_access,
                                     WeightedMean, WeightedMoment2,
                                                   WeightedMoment2> {};
template <class T>
struct WeightedMoment2
    : public detail::gen_sum_squared_diff_acc<T,
                                            detail::accumulable_weighted_access,
                                            WeightedMoment2, WeightedMean> {};

} // namespace acc

namespace detail {

// mean and second moment for templated things [design study]
// - to use like temp_mean<accumulable_weighted_access>::template acc
// in place of WeightedMean
// - sample use may be accumulable_weighted_grey_access / etc.

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

template <class T, template<class> class CV,
                   template<class> class SUM_SQUARED_DIFF>
struct gen_variance_calc
    : public calculator_base<typename CV<T>::type>,
      public type_lists::implies_template<T, SUM_SQUARED_DIFF>
{
    template <class ACX>
    void calculate(const ACX & x) const // calculator tuple == x.calculators
    {
        this->value = use<SUM_SQUARED_DIFF>::val(x)
                           / use<acc::Count>::val(x);
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Variance
    : public detail::gen_variance_calc<T, detail::accumulable_value_access,
                                          Moment2> {};
template <class T>
struct CoordVariance
    : public detail::gen_variance_calc<T, detail::accumulable_coord_access,
                                          CoordMoment2> {};

template <class T>
struct WeightedVariance
    : public detail::gen_variance_calc<T, detail::accumulable_weighted_access,
                                          WeightedMoment2> {};

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
        this->value = use<SUM_SQUARED_DIFF>::val(x)
                          / (use<acc::Count>::val(x) - 1.0);
    }
};

} // namespace detail

namespace acc {

template <class T>
struct UnbiasedVariance
    : public detail::gen_unbiased_variance_calc<T,
                                               detail::accumulable_value_access,
                                               Moment2> {};
template <class T>
struct CoordUnbiasedVariance
    : public detail::gen_unbiased_variance_calc<T,
                                               detail::accumulable_coord_access,
                                               CoordMoment2> {};
template <class T>
struct WeightedUnbiasedVariance
    : public detail::gen_unbiased_variance_calc<T,
                                            detail::accumulable_weighted_access,
                                            WeightedMoment2> {};

} // namespace acc

namespace detail {

template <class T, template<class> class CV, template<class> class SELF>
struct gen_sum_acc : public numeric_accumulator_base<T, typename CV<T>::type>
{
    /** update sum
    */
    template <class ACX>
    void updatePass1(ACX &, const T & v)
    {
        this->value += CV<T>::get(v);
    }
    /** update average, using weighted input.
     * <tt>stats(value, 1.0)</tt> is equivalent to the unweighted
     * call <tt>stats(value)</tt>, and <tt>stats(value, 2.0)</tt>
     * is equivalent to two unweighted calls.
     */
    template <class ACX>
    void updatePass1(ACX &, const T & v, double weight)
    {
        this->value += CV<T>::get(v) * weight;
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        this->value = use<SELF>::val(x) + use<SELF>::val(y);
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Sum
    : public detail::gen_sum_acc<T, detail::accumulable_value_access, Sum> {};

template <class T>
struct CoordSum
    : public detail::gen_sum_acc<T, detail::accumulable_coord_access, CoordSum>
{};

template <class T>
struct WeightedSum
    : public detail::gen_sum_acc<T, detail::accumulable_weighted_access,
                                    WeightedSum> {};

} // namespace acc

namespace detail {

template <unsigned n, class T, template<class> class CV,
                   template<class> class MEAN,
                   template<class> class SELF>
struct gen_central_moment_acc2
    : public numeric_accumulator_base<T, typename CV<T>::type>,
      public type_lists::implies_template<T, MEAN>
{
    static const unsigned number_of_passes = 2;
    // reset() and unimplemented variants
    using numeric_accumulator_base<T, typename CV<T>::type>::operator();

    template <class ACX> void updatePass1(ACX &, const T &) {} // pass 1
    
    template <class ACX>
    void updatePassN(ACX & x, const T & v, unsigned pass) // pass 2
    {
        if (test_pass<ACX, 2>(pass))
            this->value += power<n>(CV<T>::get(v) - use<MEAN>::f_val(x));
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void operator()(const ACX & z, const ACX & x, const ACX & y)
    {
        this->value =   use<SELF>::val(x)
                      + use<SELF>::val(y);
    }
};

template <unsigned n, class T, template<class> class SELF>
struct central_moment_acc2
    : public gen_central_moment_acc2<n, T, accumulable_value_access,
                                           acc::Mean, SELF> {};
template <unsigned n, class T, template<class> class SELF>
struct coord_central_moment_acc2
    : public gen_central_moment_acc2<n, T, accumulable_coord_access,
                                           acc::CoordMean, SELF> {};
template <unsigned n, class T, template<class> class SELF>
struct weighted_central_moment_acc2
    : public gen_central_moment_acc2<n, T, accumulable_weighted_access,
                                           acc::WeightedMean, SELF> {};

} // namespace detail

namespace acc {

template <class T>
struct Moment2_2 : public detail::central_moment_acc2<2, T, acc::Moment2_2> {};

template <class T>
struct Moment3 : public detail::central_moment_acc2<3, T, Moment3> {};

template <class T>
struct Moment4 : public detail::central_moment_acc2<4, T, Moment4> {};

template <class T>
struct CoordMoment2_2
    : public detail::coord_central_moment_acc2<2, T, CoordMoment2> {};

template <class T>
struct CoordMoment3
    : public detail::coord_central_moment_acc2<3, T, CoordMoment3> {};

template <class T>
struct CoordMoment4
    : public detail::coord_central_moment_acc2<4, T, CoordMoment4> {};

template <class T>
struct WeightedMoment2_2
    : public detail::weighted_central_moment_acc2<2, T, WeightedMoment2> {};

template <class T>
struct WeightedMoment3
    : public detail::weighted_central_moment_acc2<3, T, WeightedMoment3> {};

template <class T>
struct WeightedMoment4
    : public detail::weighted_central_moment_acc2<4, T, WeightedMoment4> {};

} // namespace acc

namespace detail {

template <class T, template<class> class CV,
                   template<class> class M2_ACC2,
                   template<class> class M3_ACC2>
struct gen_skewness_calc2
    : public calculator_base<typename CV<T>::type>,
      public type_lists::implies_template<T, M2_ACC2, M3_ACC2>
{
    template <class ACX>
    void calculate(const ACX & x) const
    {
        using std::sqrt;
        this->value = use<M3_ACC2>::val(x)
                      / sqrt(power<3>(use<M2_ACC2>::val(x)));
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Skewness
    : public detail::gen_skewness_calc2<T, detail::accumulable_value_access,
                                           Moment2_2, Moment3> {};
template <class T>
struct CoordSkewness
    : public detail::gen_skewness_calc2<T, detail::accumulable_coord_access,
                                           CoordMoment2_2, CoordMoment3> {};
template <class T>
struct WeightedSkewness
    : public detail::gen_skewness_calc2<T, detail::accumulable_weighted_access,
                                           WeightedMoment2_2, WeightedMoment3>
{};

} // namespace acc

namespace detail {

template <class T, template<class> class CV,
                   template<class> class M2_ACC2,
                   template<class> class M4_ACC2>
struct gen_kurtosis_calc2
    : public calculator_base<typename CV<T>::type>,
      public type_lists::implies_template<T, M2_ACC2, M4_ACC2>
{
    template <class ACX>
    void calculate(const ACX & x) const
    {
        this->value = use<M4_ACC2>::val(x) / power<2>(use<M2_ACC2>::val(x));
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Kurtosis
    : public detail::gen_kurtosis_calc2<T, detail::accumulable_value_access,
                                           Moment2_2, Moment4> {};
template <class T>
struct CoordKurtosis
    : public detail::gen_kurtosis_calc2<T, detail::accumulable_coord_access,
                                           CoordMoment2_2, CoordMoment4> {};
template <class T>
struct WeightedKurtosis
    : public detail::gen_kurtosis_calc2<T, detail::accumulable_weighted_access,
                                           WeightedMoment2_2, WeightedMoment4>
{};

} // namespace acc

namespace detail {

template <class T, template<class> class CV, template<class> class SELF>
struct gen_min_acc : public collector_base<T, typename CV<T>::type>
{
    typedef typename gen_min_acc::collector_base_type collector_base_type;
    void reset()
    {
        this->value = NumericTraits<typename CV<T>::type>::max();
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
        updateMin(this->value, CV<T>::get(v));
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        this->value = use<SELF>::val(x);
        updateMin(this->value, use<SELF>::val(y));
    }
};

template <class T, template<class> class CV, template<class> class SELF>
struct gen_max_acc : public collector_base<T, typename CV<T>::type>
{
    typedef typename gen_max_acc::collector_base_type collector_base_type;
    void reset()
    {
        this->value = NumericTraits<typename CV<T>::type>::min();
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
        updateMax(this->value, CV<T>::get(v));
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void merge(const ACX & z, const ACX & x, const ACX & y)
    {
        this->value = use<SELF>::val(x);
        updateMax(this->value, use<SELF>::val(y));
    }
};

} // namespace detail

namespace acc {

template <class T>
struct Min : public gen_min_acc<T, detail::accumulable_value_access, Min> {};
template <class T>
struct CoordMin
    : public gen_min_acc<T, detail::accumulable_coord_access, CoordMin> {};
template <class T>
struct WeightedMin
    : public gen_min_acc<T, detail::accumulable_weighted_access, WeightedMin>
{};

template <class T>
struct Max : public gen_max_acc<T, detail::accumulable_value_access, Max> {};
template <class T>
struct CoordMax
    : public gen_max_acc<T, detail::accumulable_coord_access, CoordMax> {};
template <class T>
struct WeightedMax
    : public gen_max_acc<T, detail::accumulable_weighted_access, WeightedMax>
{};

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

template <class T, class L,
          template<class, template<class> class TES>
          class TUPLE = type_lists::tuple>
struct feature_context_base
    : public TUPLE<L, passes<2>::template at_least>
{
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

    template <class, class A> void calculate(A &) const {}
    template <template<class> class, class A> void calculate_template(A &) const
    {}
};

template <class T, class L>
struct calc_context : public feature_context_base<T, L>
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
};

template <class T, class L, template<class, template<class> class TEST>
                            class TUPLE = type_lists::tuple>
struct accu_context_base
    : public feature_context_base<T, typename feature_types<L>::accu_list,
                                  TUPLE>,
      public UnaryAnalyserTag
{
    typedef typename accu_context_base::tuple_type                base_type;
    typedef typename accu_context_base::list_type                 list_type;
    typedef typename accu_context_base::feature_context_base_type feature_cx;
    typedef          accu_context_base                            extra_passes;
    typedef          feature_types<L>                             feature_lists;     
    static const unsigned max_passes = feature_lists::max_passes;

    // features that do not require even a single pass:
    typedef typename feature_lists::sync_list                     sync_list;
    typedef calc_context<T, typename feature_lists::calc_list>    calc_type;
    calc_type                                                     calculators;

    template <template<class> class V> // works around cl.exe bugs
    struct has_feature
    {
        static const bool value
            =  base_type::template has_element< V<T> >::value
            || calc_type::template has_element< V<T> >::value;
    };

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
    typename chooser<V>::ref_type ref() const
    {
        return chooser<V>::ref(*this, calculators);
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

    template <template<class> class V>
    void set_template()
    {
        this->template set< V<T> >();
    }

    template <template<class> class V>
    bool is_set_template() const
    {
        return chooser<V<T> >::at(*this, calculators).template is_set<V<T> >();
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
    /** update using weighted input
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

    /** merge two statistics: *this <- (x, y)
    */
    void operator()(const accu_context_base & x0, const accu_context_base & y0)
    {
        // first, merge set features (if any):
        x0.transfer_set_to(*this);
        y0.transfer_set_to(*this);
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

    // compatibility with ArrayOfRegionStatistics
    typedef T argument_type;
    typedef void result_type;
    result_type operator()() const {}
};

} // namespace detail

// The feature accessing functions get(), select(), and member() are
// implemented as free template functions, since an implementation as member
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
struct SelectedAccumulators
    : public detail::accu_context_base<T, L,
        type_lists::cond_virtual_tuple<SelectedAccumulators, T,
                                       detail::v_feature_dispatch>::template
                                                                    type>
{};


} // namespace vigra

#endif // VIGRA_OBJECT_FEATURES_HXX
