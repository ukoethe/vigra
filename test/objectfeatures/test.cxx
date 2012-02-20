#include <iostream>
#include <sstream>
#include <map>
#include <set>

#include <unittest.hxx>

#include <vigra/accessor.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/rgbvalue.hxx>

#include <vigra/coordinate_iterator.hxx>
#include <vigra/object_features.hxx>

#include <vigra/multi_pointoperators.hxx>

#include <vigra/basicimage.hxx>
#include <vigra/stdimage.hxx> // BImage
#include <vigra/inspectimage.hxx> // FindAverageAndVariance

#include <vigra/multi_array.hxx>
#include <vigra/multi_convolution.hxx>
#include <vigra/impex.hxx>
#include <vigra/imageinfo.hxx>
#include <vigra/functorexpression.hxx>

#include <vigra/algorithm.hxx>
#include <vigra/histogram.hxx>
#include <vigra/random.hxx>
#include <vigra/convolution.hxx>
#include <vigra/accumulator.hxx>

namespace std {

template <unsigned int N, class T, class Stride>
ostream & operator<<(ostream & o, vigra::MultiArrayView<N, T, Stride> const & m)
{
    for(MultiArrayIndex k=0; k<m.size(); ++k)
        o << m[k] << " ";
    return o;
}

} // namespace std

using namespace vigra;

// mask cl.exe shortcomings
#if defined(_MSC_VER)
#pragma warning( disable : 4503 )
#endif

#if 0
template <class X, unsigned N>
TinyVector<X, N> & tab_assign(TinyVector<X, N> & x, const X y[N])
{
    for (unsigned k = 0; k != N; ++k)
        x[k] = y[k];
    return x;
}

template <class X>
X nan_squash(const X & y)
{
    if (y == y)
        return y;
    else
        return X(-0.0);
}
template <class X, unsigned N>
TinyVector<X, N> nan_squash(const TinyVector<X, N> & y)
{
    TinyVector<X, N> x;
    for (unsigned k = 0; k != N; ++k)
        x[k] = y[k] == y[k] ? y[k] : -0.0;
    return x;
}

using vigra::type_lists::use_template_list;
using vigra::Accumulators;

typedef double pixel_type;

typedef vigra::MultiArrayShape<2>::type  shape_2d;
typedef vigra::MultiArray<2, pixel_type> array_2d;

typedef vigra::RGBValue<double> rgb_type;

template <class T>
struct gen_accumulators
    : public use_template_list<Accumulators, T,
                               acc::Variance,
                               acc::Skewness,
                               acc::Kurtosis
                              >
{};

typedef gen_accumulators<pixel_type> pixel_accumulators;

typedef use_template_list<Accumulators, StridePairPointer<2, pixel_type>,
                          acc::Count,
                          acc::Sum,
                          acc::Moment2,
                          acc::CoordSum,
                          acc::CoordMean,
                          acc::CoordMoment2,
                          acc::WeightedSum,
                          acc::WeightedMean,
                          acc::WeightedMoment2
                         >
    cv_test_accumulators;

typedef use_template_list<Accumulators, pixel_type,
                          acc::Count,
                          acc::Sum,
                          acc::Mean,
                          acc::Moment2
                         >
yy_test_accumulators;

typedef use_template_list<Accumulators, pixel_type,
                          acc::Mean,
                          acc::Moment2
                         >
yyy_test_accumulators;

struct count : public InitializerTag
{
    mutable unsigned long t;
    count() : t(0) {}
    unsigned long operator()() const
    {
        return t++;
    }
};

void test1()
{
    pixel_accumulators pix_acc;
    yyy_test_accumulators yy_acc;
    
    array_2d image(array_2d::size_type(10, 10));

    for (unsigned i = 0; i != 100; ++i)
        image(i / 10, i % 10) = i;

    inspectMultiArray(srcMultiArrayRange(image), pix_acc);

    shouldEqual(get<acc::Count>(pix_acc), 100);
    shouldEqualTolerance(get<acc::Mean>(pix_acc), 49.5, 2e-15);
    shouldEqualTolerance(get<acc::Moment2>(pix_acc), 83325, 2e-15);
    shouldEqualTolerance(get<acc::Moment2_2Pass>(pix_acc), 83325, 2e-15);
    shouldEqualTolerance(get<acc::Variance>(pix_acc), 833.25, 2e-15);
    shouldEqualTolerance(get<acc::Skewness>(pix_acc), 0, 2e-15);
    shouldEqualTolerance(get<acc::Kurtosis>(pix_acc), 0.017997599759976, 2e-15);
}

void test2()
{
    array_2d stride_image2(array_2d::size_type(2, 3));
    initMultiArray(destMultiArrayRange(stride_image2), count());

    cv_test_accumulators cv_acc;

    inspectMultiArray(srcCoordinateMultiArrayRange(stride_image2), cv_acc);

    shouldEqual(get<acc::Count>(cv_acc), 6);
    shouldEqual(get<acc::Sum>(cv_acc), 15);
    { TinyVector<double, 2> x(3, 6); shouldEqualSequence(x.begin(), x.end(), get<acc::CoordSum>(cv_acc).begin()); };
    { TinyVector<double, 2> x(0.5, 1); shouldEqualSequenceTolerance(x.begin(), x.end(), get<acc::CoordMean>(cv_acc).begin(), 2e-15); };
    { TinyVector<double, 2> x(1.5, 4); shouldEqualSequenceTolerance(x.begin(), x.end(), get<acc::CoordMoment2>(cv_acc).begin(), 2e-15); };
    shouldEqualTolerance(get<acc::Mean>(cv_acc), 2.5, 2e-15);
    shouldEqualTolerance(get<acc::Moment2>(cv_acc), 17.5, 2e-15);
    { TinyVector<double, 2> x(9, 23); shouldEqualSequence(x.begin(), x.end(), get<acc::WeightedSum>(cv_acc).begin()); };
    { TinyVector<double, 2> x(1.5, 3.83333333333333); shouldEqualSequenceTolerance(x.begin(), x.end(), get<acc::WeightedMean>(cv_acc).begin(), 2e-15); };
    { TinyVector<double, 2> x(21.5, 88.8333333333333); shouldEqualSequenceTolerance(x.begin(), x.end(), get<acc::WeightedMoment2>(cv_acc).begin(), 2e-15); };

    rgb_type epstest = 1e-5;
    { TinyVector<double, 3> x(1e-05, 1e-05, 1e-05); shouldEqualSequenceTolerance(x.begin(), x.end(), epstest.begin(), 2e-15); };
}


using vigra::NestedAccumulators;
using vigra::SelectedAccumulators;

template <class T, template<class> class CV,
                   template<class> class SUM>
struct gen_mean_from_sum_calc
    : public detail::calculator_base<typename CV<T>::type>,
      public type_lists::implies_template<T, SUM, acc::Count>
{
    template <class ACX>
    void calculate(const ACX & x) const // calculator tuple == x.calculators
    {
        this->set_val(detail::use<SUM>::val(x)
                      / get<acc::Count>(x));
    }
};
template <class T>
struct mean_from_sum_calc
    : public gen_mean_from_sum_calc<T, detail::accumulable_value_access,
                                       acc::Sum> {};
template <class T>
struct coord_mean_from_sum_calc
    : public gen_mean_from_sum_calc<T, detail::accumulable_coord_access,
                                       acc::CoordSum>
{};

// testing contexts:

template <class T, class L>
struct cond_plain_accu_context
    : public detail::accu_context_base<T, L, type_lists::cond_tuple_plain> {};

template <class T, class L>
struct cond_accu_context
    : public detail::accu_context_base<T, L, type_lists::cond_tuple> {};

template <class T, class L>
struct virtual_accu_context
    : public detail::accu_context_base<T, L,
              type_lists::virtual_tuple<virtual_accu_context, T,
                                        detail::v_feature_dispatch>::template
                                                                     type>
{};

template <class Z> const Z & make_const(const Z & z) { return z; }

template <class T>
struct mean_context
    : public use_template_list<Accumulators, T, acc::Count, acc::Mean>
{
};

template <class T>
struct nested_accumulators_mean
    : public use_template_list<NestedAccumulators, T,
                                    acc::Count, acc::Mean>
{
};

template <class T>
struct selected_accumulators_mean
    : public use_template_list<SelectedAccumulators, T, acc::Mean>
{
};

// sum of squares

template <class T>
struct accumulators_sq
    : public use_template_list<Accumulators, T, acc::Moment2>
{
};

template <class T>
struct accumulators_sq2
    : public use_template_list<Accumulators, T, acc::Variance,
                                                     acc::UnbiasedVariance>
{
};

template <class T>
struct nested_accumulators_sq
    : public use_template_list<NestedAccumulators, T,
                                    acc::Moment2>
{
};

template <class T>
struct selected_accumulators_sq
    : public use_template_list<SelectedAccumulators, T,
                                    acc::Moment2>
{
};

template <class T>
struct cond_plain_mean_context
    : public use_template_list<cond_plain_accu_context, T,
                                    acc::Count, acc::Mean>
{
};

template <class T>
struct cond_mean_context
    : public use_template_list<cond_accu_context, T,
                                    acc::Count, acc::Mean>
{
};

template <class T>
struct accumulators_mean
    : public use_template_list<Accumulators, T, acc::Mean>
{
};

template <class T>
struct virtual_mean_context
    : public use_template_list<virtual_accu_context, T, acc::Mean>
{
};

template <class T>
struct accumulators_mean_sum
    : public use_template_list<Accumulators, T, acc::Sum, mean_from_sum_calc,
                               acc::Moment2_2Pass, acc::Variance>
{};

template <class key, class data>
struct gen_store {
    typedef std::map<key, data> gen_map;
    typedef typename gen_map::iterator    iterator;
    typedef typename gen_map::size_type    size_type;
    typedef typename gen_map::mapped_type    data_type; 
    typedef typename gen_map::key_type    key_type;
    typedef typename gen_map::value_type    value_type;
    typedef typename gen_map::const_iterator    const_iterator;
    typedef pair<data, bool> join_type;

    gen_map store_map;

    bool find(const key & x, data & id) const {
        const_iterator i = store_map.find(x);
        bool found = (i != store_map.end());
        if (found)
            id = (*i).second;
        return found;
    }
    bool exists(const key & x) const {
        return store_map.find(x) != store_map.end();
    }
};

template <class search, class id_type = unsigned>
class label_map
{
private:
    typedef std::map<search, id_type>        gen_map;
    typedef typename gen_map::iterator       iterator;
    typedef typename gen_map::size_type      size_type;
    typedef typename gen_map::mapped_type    data_type; 
    typedef typename gen_map::key_type       key_type;
    typedef typename gen_map::value_type     value_type;
    typedef typename gen_map::const_iterator const_iterator;
    typedef          search                  reverse_type;
    typedef pair<id_type, bool> join_type;

    gen_map                   store_map;
    id_type                   id_count;
    std::vector<reverse_type> reverse_map;

public:
    bool find(const search & x, id_type & id) const
    {
        const_iterator i = store_map.find(x);
        bool found = (i != store_map.end());
        if (found)
            id = (*i).second;
        return found;
    }
    bool exists(const search & x) const
    {
        return store_map.find(x) != store_map.end();
    }
private:
    id_type push_back(const search & x = search())
    {
        reverse_map.push_back(x);
        return id_count++;
    }
    join_type join(const search & x)
    {
        std::pair<iterator, bool> inserted =
            store_map.insert(value_type(x, id_count));
        if (inserted.second)
            push_back(x);
        // return (id number, was it new?)
        return join_type((*(inserted.first)).second, inserted.second);
    }
public:
    label_map() : id_count(0) {}
    void clear()
    {
        id_count = 0;
        store_map.clear();
        reverse_map.clear();
    }
    id_type count() const { return id_count; }
    // adds x into symbol map if new, returns id #
    id_type operator()(const search & x)
    {
        return join(x).first;
    }
    // transform to colour label
    unsigned operator()(const search & v) const
    {
        unsigned id;
        if (find(v, id))
        {
            return id;
        }
        else
            return ~(unsigned(0));
    }
    // adds x into symbol map if new, returns "if new" and sets id
    bool subjoin(const search & x, id_type & id)
    {
        join_type j = join(x);
        id = j.first;
        return j.second;
    }
    bool subjoin(const search & x)
    {
        join_type j = join(x);
        return j.second;
    }
    const reverse_type & operator[](id_type id) const
    {
        return reverse_map[id];
    }
    reverse_type & operator[](id_type id)
    {
        return reverse_map[id];
    }
};

template <class T>
struct colour_label_acc : public detail::collector_base<T, label_map<T> >
{
    bool test_member;
    colour_label_acc() : test_member(false) {}
    
    typedef typename colour_label_acc::collector_base_type collector_base_type;
    void reset()
    {
        this->value.clear();
    }
    // reset() and unimplemented variants
    using collector_base_type::operator();
    /** update
    */
    template <class ACX>
    void updatePass1(ACX &, const T & v)
    {
        this->value(v);
    }
    /** merge two statistics: z === *this <- (x, y)
    */
    template <class ACX>
    void operator()(const ACX & z, const ACX & x, const ACX & y);
};

typedef vigra::RGBValue<double> vpixel_type;

typedef vigra::MultiArrayShape<2>::type                          shape_v2d;
typedef vigra::MultiArray<2, vpixel_type>                        array_v2d;
typedef vigra::MultiArray<2, int>                                labels_v2d;
typedef vigra::MultiArray<2, vigra::TinyVector<vpixel_type, 2> > grad_v2d;
typedef vigra::MultiArray<2, vigra::TinyVector<vpixel_type, 3> > symm_v2d;
typedef vigra::BasicImageView<vpixel_type>                       image_v2d;
typedef vigra::ConvolutionOptions<2>                             options_v2d;

template <class T>
struct accumulators_colour_label
    : public use_template_list<Accumulators, T, colour_label_acc>
{};

void test3()
{
    BImage img(10, 10);
    for (unsigned i = 0; i != 100; ++i)
        img(i / 10, i % 10) = i;

    FindAverageAndVariance<BImage::PixelType> averageAndVariance; // init functor
    mean_context<BImage::PixelType> my_mean;

    cond_plain_mean_context<BImage::PixelType>    my_c_mean;
    cond_mean_context<BImage::PixelType>          my_c_mean_new;
    nested_accumulators_mean<BImage::PixelType>   my_c_mean_bit;
    selected_accumulators_mean<BImage::PixelType> my_v_mean; 
    accumulators_mean<BImage::PixelType>          test_mean;
    cond_plain_mean_context<BImage::PixelType>    test_c_mean;
    cond_mean_context<BImage::PixelType>          test_c_mean_new;
    nested_accumulators_mean<BImage::PixelType>   test_c_mean_bit;
    selected_accumulators_mean<BImage::PixelType> test_v_mean;

    accumulators_mean_sum<BImage::PixelType>      mean_sum_0;
    accumulators_colour_label<BImage::PixelType>  colour_0;

    virtual_mean_context<BImage::PixelType> virtual_meanest;
    select<acc::Mean>(virtual_meanest);

    virtual_mean_context<BImage::PixelType> virtual_empty;

    accumulators_sq<BImage::PixelType> my_sq;
    accumulators_sq2<BImage::PixelType> my_sq2;
    // at least 1 inpect pass:
    shouldEqual(detail::passes<1>::at_least<acc::Count<BImage::PixelType> >::value, true);
    shouldEqualTolerance(detail::passes<1>::at_least<acc::Variance<BImage::PixelType> >::value, 0, 2e-15);
    
    shouldEqualTolerance(detail::use<acc::Count>::f_val(my_sq), 0, 2e-15);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_sq), 0, 2e-15);
    shouldEqualTolerance(detail::use<acc::Moment2>::f_val(my_sq), 0, 2e-15);

    shouldEqualTolerance(detail::use<acc::Count>::f_val(my_c_mean), 0, 2e-15);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_c_mean), 0, 2e-15);
    shouldEqualTolerance(detail::use<acc::Count>::f_val(my_c_mean_new), 0, 2e-15);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_c_mean_new), 0, 2e-15);
    shouldEqualTolerance(detail::use<acc::Count>::f_val(my_c_mean_bit), 0, 2e-15);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_c_mean_bit), 0, 2e-15);
    shouldEqualTolerance(get<acc::Count>(my_v_mean), 0, 2e-15); //getting virtual stuff
    shouldEqualTolerance(get<acc::Mean>(my_v_mean), 0, 2e-15);  //sets them, thus val

    select<acc::Count>(my_c_mean); // only use count...
    select<acc::Count>(my_c_mean_new); // only use count...
    select<acc::Count>(my_c_mean_bit); // only use count...
    select<acc::Count>(my_v_mean); // only use count...
    
    select<acc::Mean>(test_c_mean);
    select<acc::Mean>(test_v_mean);
    select<acc::Count>(test_v_mean); // test idempotence
    select<acc::Mean>(test_c_mean_new);
    select<acc::Mean>(test_c_mean_bit);

    inspectImage(srcImageRange(img), averageAndVariance);
    inspectImage(srcImageRange(img), my_sq);
    inspectImage(srcImageRange(img), my_sq2);
    inspectImage(srcImageRange(img), my_mean);
    inspectImage(srcImageRange(img), my_c_mean);
    inspectImage(srcImageRange(img), my_c_mean_new);
    inspectImage(srcImageRange(img), my_c_mean_bit);
    inspectImage(srcImageRange(img), my_v_mean);
    inspectImage(srcImageRange(img), test_mean);
    inspectImage(srcImageRange(img), test_c_mean);
    inspectImage(srcImageRange(img), test_c_mean_new);
    inspectImage(srcImageRange(img), test_c_mean_bit);
    inspectImage(srcImageRange(img), test_v_mean);

    inspectImage(srcImageRange(img), mean_sum_0);
    inspectImage(srcImageRange(img), colour_0);

    inspectImage(srcImageRange(img), virtual_meanest);
    inspectImage(srcImageRange(img), virtual_empty);

    shouldEqual(averageAndVariance.count_, 100);
    shouldEqualTolerance(averageAndVariance.mean_, 49.5, 2e-15);
    shouldEqualTolerance(averageAndVariance.sumOfSquaredDifferences_, 83325, 2e-15);
 
    shouldEqual(detail::use<acc::Count>::f_val(my_sq2), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_sq2), 49.5, 2e-15);
    shouldEqualTolerance(detail::use<acc::Moment2>::f_val(my_sq2), 83325, 2e-15);
    shouldEqualTolerance(get<acc::Variance>(my_sq2), 833.25, 2e-15);
    shouldEqualTolerance(get<acc::UnbiasedVariance>(my_sq2), 841.666666666666, 2e-15);

    shouldEqual(detail::use<acc::Count>::f_val(my_sq), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_sq), 49.5, 2e-15);
    shouldEqualTolerance(detail::use<acc::Moment2>::f_val(my_sq), 83325, 2e-15);
    shouldEqualTolerance(std::sqrt(detail::use<acc::Moment2>::f_val(my_sq) / detail::use<acc::Count>::f_val(my_sq)), 28.8660700477221, 2e-15);

    shouldEqual(detail::use<acc::Count>::f_val(my_mean), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_mean), 49.5, 2e-15);

    shouldEqual(detail::use<acc::Count>::f_val(my_c_mean), 100);
    // mean must be zero as only count is selected
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_c_mean), 0, 2e-15);
    shouldEqual(detail::use<acc::Count>::f_val(my_c_mean_new), 100);
    // mean must be zero as only count is selected
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_c_mean_new), 0, 2e-15);
    shouldEqual(detail::use<acc::Count>::f_val(my_c_mean_bit), 100);
    // mean must be zero as only count is selected
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_c_mean_bit), 0, 2e-15);
    shouldEqual(get<acc::Count>(my_v_mean), 100);
    // mean must be zero as only count is selected
    shouldEqualTolerance(get<acc::Mean>(my_v_mean), 0, 2e-15);

    shouldEqual(detail::use<acc::Count>::f_val(test_mean), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(test_mean), 49.5, 2e-15);
    
    shouldEqual(detail::use<acc::Count>::f_val(test_c_mean), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(test_c_mean), 49.5, 2e-15);
    shouldEqual(detail::use<acc::Count>::f_val(test_c_mean_new), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(test_c_mean_new), 49.5, 2e-15);
    shouldEqual(detail::use<acc::Count>::f_val(test_c_mean_bit), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(test_c_mean_bit), 49.5, 2e-15);
    shouldEqual(detail::use<acc::Count>::f_val(test_v_mean), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(test_v_mean), 49.5, 2e-15);

    shouldEqual(detail::use<acc::Count>::f_val(virtual_meanest), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(virtual_meanest), 49.5, 2e-15);

    // both must be zero as nothing is selected
    shouldEqualTolerance(get<acc::Count>(virtual_empty), 0, 2e-15);
    shouldEqualTolerance(get<acc::Mean>(virtual_empty), 0, 2e-15);

    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean), 100);

    mean_context<BImage::PixelType> my_mean2 = my_mean;
    cond_plain_mean_context<BImage::PixelType> my_c_mean2 = my_c_mean;
    cond_mean_context<BImage::PixelType> my_c_mean_new2 = my_c_mean_new;
    nested_accumulators_mean<BImage::PixelType> my_c_mean_bit2 = my_c_mean_bit;

    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean), 100);

    my_mean2(my_mean2, my_mean);
    shouldEqual(detail::use<acc::Count>::f_val(my_mean2), 200);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_mean2), 49.5, 2e-15);
    mean_context<BImage::PixelType> my_mean3;
    my_mean3(my_mean, my_mean3);
    shouldEqual(detail::use<acc::Count>::f_val(my_mean3), 100);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_mean3), 49.5, 2e-15);
    mean_context<BImage::PixelType> my_mean4(my_mean3);
    my_mean4(my_mean2);
    shouldEqual(detail::use<acc::Count>::f_val(my_mean4), 300);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_mean4), 49.5, 2e-15);
    my_mean4.reset();
    // both must be zero as everything is reset
    shouldEqualTolerance(detail::use<acc::Count>::f_val(my_mean4), 0, 2e-15);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_mean4), 0, 2e-15);

    // begin copy constructor tests.
    mean_context<BImage::PixelType> my_mean5;
    my_mean4.transfer_set_to(my_mean5);

    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean), 100);
    shouldEqualTolerance(get<acc::Mean>(my_v_mean), 0, 2e-15);

    selected_accumulators_mean<BImage::PixelType> my_v_mean2 = my_v_mean;
    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean2), 100);
    shouldEqualTolerance(get<acc::Mean>(my_v_mean2), 0, 2e-15);
    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean), 100);
    shouldEqualTolerance(get<acc::Mean>(my_v_mean), 0, 2e-15);

    shouldEqual(detail::use<acc::Count>::f_val(test_v_mean), 100);
    shouldEqualTolerance(get<acc::Mean>(test_v_mean), 49.5, 2e-15);

    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean2), 100);
    shouldEqualTolerance(get<acc::Mean>(my_v_mean2), 0, 2e-15);
    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean), 100);
    shouldEqualTolerance(get<acc::Mean>(my_v_mean), 0, 2e-15);
    my_v_mean2(test_v_mean);
    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean), 100);
    shouldEqualTolerance(get<acc::Mean>(my_v_mean), 0, 2e-15);
    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean2), 200);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_v_mean2), 24.75, 2e-15); // not a bug, but a feature..

    shouldEqual(detail::use<acc::Count>::f_val(test_v_mean), 100);
    shouldEqualTolerance(get<acc::Mean>(test_v_mean), 49.5, 2e-15);

    shouldEqual(get<acc::Count>(my_v_mean), 100);
    shouldEqualTolerance(get<acc::Mean>(my_v_mean), 0, 2e-15);
    my_v_mean(test_v_mean);
    shouldEqual(detail::use<acc::Count>::f_val(my_v_mean), 200);
    shouldEqualTolerance(detail::use<acc::Mean>::f_val(my_v_mean), 24.75, 2e-15); // not a bug, but a feature..

    typedef selected_accumulators_mean<BImage::PixelType> cvmc;
    cvmc* my_v_mean3 = new cvmc(my_v_mean);
    shouldEqual(get<acc::Count>(*my_v_mean3), 200);
    shouldEqualTolerance(get<acc::Mean>(*my_v_mean3), 24.75, 2e-15);
    my_v_mean3->reset();
    shouldEqualTolerance(get<acc::Count>(*my_v_mean3), 0, 2e-15);
    shouldEqualTolerance(get<acc::Mean>(*my_v_mean3), 0, 2e-15);
    (*my_v_mean3)(test_v_mean);
    shouldEqual(get<acc::Count>(*my_v_mean3), 100);
    shouldEqualTolerance(get<acc::Mean>(*my_v_mean3), 49.5, 2e-15);
    cvmc* my_v_mean4 = new cvmc;
    shouldEqualTolerance(get<acc::Count>(*my_v_mean4), 0, 2e-15);
    shouldEqualTolerance(get<acc::Mean>(*my_v_mean4), 0, 2e-15);
    cvmc* my_v_mean5 = new cvmc(*my_v_mean3);
    delete my_v_mean3;
    shouldEqual(get<acc::Count>(*my_v_mean5), 100);
    shouldEqualTolerance(get<acc::Mean>(*my_v_mean5), 49.5, 2e-15);
    (*my_v_mean4)(*my_v_mean5, *my_v_mean5);
    shouldEqual(get<acc::Count>(*my_v_mean4), 200);
    shouldEqualTolerance(get<acc::Mean>(*my_v_mean4), 49.5, 2e-15);
    // end copy constructor tests.

    // calculator and member<> access function tests.
    shouldEqual(get<acc::Count>(mean_sum_0), 100);
    shouldEqual(get<acc::Sum>(mean_sum_0), 4950);
    shouldEqualTolerance(detail::use<mean_from_sum_calc>::c_val(mean_sum_0.calculators), 0, 2e-15);
    shouldEqualTolerance(detail::use<mean_from_sum_calc>::val(mean_sum_0), 49.5, 2e-15);
    shouldEqualTolerance(detail::use<mean_from_sum_calc>::c_val(mean_sum_0.calculators), 49.5, 2e-15);
    shouldEqual(accumulators_mean_sum<BImage::PixelType>::max_passes, 2);
    shouldEqual(get<acc::Moment2_2Pass>(mean_sum_0), 83325);

    shouldEqualTolerance(detail::use<acc::Variance>::c_val(mean_sum_0.calculators), 0, 2e-15);
    shouldEqualTolerance(get<acc::Variance>(mean_sum_0), 833.25, 2e-15);
    shouldEqualTolerance(detail::use<acc::Variance>::c_val(mean_sum_0.calculators), 833.25, 2e-15);

    shouldEqual(detail::use<colour_label_acc>::val(colour_0).count(), 100);

    shouldEqual(get<colour_label_acc>(colour_0).count(), 100);
    shouldEqualTolerance(member<colour_label_acc>(colour_0).test_member, 0, 2e-15);
    member<colour_label_acc>(colour_0).test_member = true;
    shouldEqualTolerance(member<colour_label_acc>(make_const(colour_0)).test_member, 1, 2e-15);
}

void test4()
{
    typedef MultiArray<2, float> Array;

    typedef use_template_list<Accumulators, float,
                              mean_from_sum_calc,
                              acc::Moment2_2Pass, acc::Variance,
                              acc::Skewness, acc::Kurtosis>
        accumulators_mean_sum;

    typedef use_template_list<Accumulators, float, acc::Variance>
        accumulators_var;

    typedef use_template_list<Accumulators, float, acc::Min, acc::Max>
        accumulators_minmax;

    accumulators_var test_acc;
    test_acc(static_cast<float>(1.2));
    test_acc(static_cast<float>(100.2));
    test_acc(1023);
    shouldEqual(get<acc::Count>(test_acc), 3);
    shouldEqualTolerance(get<acc::Mean>(test_acc), 374.799987792969, 3e-6);
    shouldEqualTolerance(get<acc::Variance>(test_acc), 211715.125, 2e-6);
    
    unsigned h = 19;
    unsigned w = 10;
    unsigned max_label = 2;

    Array image(Array::size_type(h, w));
    Array labels(Array::size_type(h, w));

    for (unsigned i = 0; i < h; ++i)
    {
        for (unsigned j = 0; j < w; ++j)
        {
            if (i <= h / 2)
            {
                image(i, j) = static_cast<float>(10 * i + j);
                labels(i, j) = 0;
            }
            else
            {
                image(i, j) = static_cast<float>(1.234);
                labels(i, j) = 1;

            }
        }

    }
    ArrayOfRegionStatistics<accumulators_minmax, float> minmax(max_label - 1);
    ArrayOfRegionStatistics<accumulators_var, float>    sum_v(max_label - 1);
    ArrayOfRegionStatistics<accumulators_mean_sum, float>
                                                        sum_aors(max_label - 1);
    ArrayOfRegionStatistics<accumulators_colour_label<float>, float>
        colour_aors(max_label - 1);

    inspectTwoMultiArrays(srcMultiArrayRange(image),
                                 srcMultiArray(labels), minmax);
    inspectTwoMultiArrays(srcMultiArrayRange(image),
                                 srcMultiArray(labels), sum_aors);
    inspectTwoMultiArrays(srcMultiArrayRange(image),
                                 srcMultiArray(labels), sum_v);
    inspectTwoMultiArrays(srcMultiArrayRange(image),
                                 srcMultiArray(labels), colour_aors);

    shouldEqualTolerance(get<acc::Min>(minmax[0]), 0, 2e-6);
    shouldEqualTolerance(get<acc::Max>(minmax[0]), 99, 2e-6);
    shouldEqual(get<acc::Count>(sum_aors[0]), 100);
    shouldEqualTolerance(get<acc::Mean>(sum_aors[0]), 49.5, 2e-6);
    shouldEqualTolerance(get<acc::Moment2>(sum_aors[0]), 83324.9921875, 2e-6);
    shouldEqualTolerance(get<acc::Variance>(sum_aors[0]), 833.249938964844, 2e-6);
    shouldEqual(get<acc::Sum>(sum_aors[0]), 4950);
    shouldEqualTolerance(detail::use<mean_from_sum_calc>::val(sum_aors[0]), 49.5, 2e-6);
    shouldEqual(get<acc::Moment2_2Pass>(sum_aors[0]), 83325);
    shouldEqualTolerance(get<acc::Skewness>(sum_aors[0]), 0, 2e-6);
    shouldEqualTolerance(get<acc::Kurtosis>(sum_aors[0]), 0.0179976038634777, 2e-6);
    shouldEqualTolerance(get<acc::Mean>(sum_v[0]), 49.5, 2e-6);
    shouldEqualTolerance(get<acc::Variance>(sum_v[0]), 833.249938964844, 2e-6);
    shouldEqual(detail::use<colour_label_acc>::val(colour_aors[0]).count(), 100);

    shouldEqualTolerance(get<acc::Min>(minmax[1]), 1.23399996757507, 3e-6);
    shouldEqualTolerance(get<acc::Max>(minmax[1]), 1.23399996757507, 3e-6);
    shouldEqual(get<acc::Count>(sum_aors[1]), 90);
    shouldEqualTolerance(get<acc::Mean>(sum_aors[1]), 1.23399996757507, 3e-6);
    shouldEqualTolerance(get<acc::Moment2>(sum_aors[1]), 0, 2e-6);
    shouldEqualTolerance(get<acc::Variance>(sum_aors[1]), 0, 2e-6);
    shouldEqualTolerance(get<acc::Sum>(sum_aors[1]), 111.060066223145, 5e-6);
    shouldEqualTolerance(detail::use<mean_from_sum_calc>::val(sum_aors[1]), 1.23400068283081, 2e-6);
    shouldEqualTolerance(get<acc::Moment2_2Pass>(sum_aors[1]), 0, 2e-6);
    shouldEqualTolerance(get<acc::Mean>(sum_v[1]), 1.23399996757507, 3e-6);
    shouldEqualTolerance(get<acc::Variance>(sum_v[1]), 0, 2e-6);
    shouldEqual(detail::use<colour_label_acc>::val(colour_aors[1]).count(), 1);
}

void test5()
{
    vigra::ImageImportInfo import_info("of.gif");
    array_v2d test_image(shape_v2d(import_info.width(), import_info.height()));
    vigra::importImage(import_info, destImage(test_image));

    accumulators_colour_label<vpixel_type> colour_labels;
    inspectMultiArray(srcMultiArrayRange(test_image), colour_labels);
    shouldEqual(detail::use<colour_label_acc>::val(colour_labels).count(), 12);

    labels_v2d test_labels(shape_v2d(import_info.width(), import_info.height()));
    vigra::transformMultiArray(srcMultiArrayRange(test_image),
                               destMultiArray(test_labels),
                               detail::use<colour_label_acc>::val(colour_labels));
    shouldEqual(detail::use<colour_label_acc>::val(colour_labels).count(), 12);
    unsigned c_count = detail::use<colour_label_acc>::val(colour_labels).count();
    vigra::exportImage(srcImageRange(test_labels),
                       vigra::ImageExportInfo("of_labels.gif"));

    array_v2d test_image_gs(test_image);
    vigra::gaussianSmoothMultiArray(vigra::srcMultiArrayRange(test_image),
                                    vigra::destMultiArray(test_image_gs),
                                    options_v2d().stdDev(3, 5));

    using namespace vigra::functor;
    vigra::combineThreeMultiArrays(srcMultiArrayRange(test_image),
                                   srcMultiArray(test_image_gs),
                                   srcMultiArray(test_labels),
                                   destMultiArray(test_image),
                                   ifThenElse(!(Arg3() % Param(4)),
                                              Arg1(), Arg2()));

    vigra::exportImage(srcImageRange(test_image),
                       vigra::ImageExportInfo("of_gs.gif"));


    typedef use_template_list<Accumulators, StridePairPointer<2, vpixel_type>,
            acc::Variance,
            acc::UnbiasedVariance,
            acc::Skewness,
            acc::Kurtosis,
            acc::Sum,
            acc::Min,
            acc::Max,
            acc::CoordVariance,
            acc::CoordUnbiasedVariance,
            acc::CoordSkewness,
            acc::CoordKurtosis,
            acc::CoordSum,
            acc::CoordMin,
            acc::CoordMax,
            acc::WeightedVariance,
            acc::WeightedUnbiasedVariance,
            acc::WeightedSkewness,
            acc::WeightedKurtosis,
            acc::WeightedSum,
            acc::WeightedMin,
            acc::WeightedMax>
        accumulators_v2d;

    vigra::ArrayOfRegionStatistics<accumulators_v2d> test_aors(c_count);

    shouldEqual(accumulators_v2d::max_passes, 2);
    shouldEqual(acc::Moment2_2Pass<vpixel_type>::number_of_passes, 2);
    shouldEqual((type_lists::max_value<detail::passes_number, unsigned, accumulators_v2d::list_type>::value), 2);
    
    shouldEqual(vigra::ArrayOfRegionStatistics<accumulators_v2d>::max_passes, 2);

    inspectTwoMultiArrays(srcCoordinateMultiArrayRange(test_image),
                                 srcMultiArray(test_labels), test_aors);

    for (unsigned k = 0; k != c_count; ++k)
    {
        {
            double x_tab[] =
            {
                44251,
                467,
                4716,
                2061,
                1321,
                491,
                2331,
                2531,
                1386,
                1935,
                1209,
                1301
            };
            double x = x_tab[(k)];
            shouldEqual(nan_squash(get<acc::Count>(test_aors[(k)])), x);
        }
        {
            double x_tab[][3] =
            {
                { 251, 253, 250 },
                { 190.07965721932, 209.987419616725, 181.582565021942 },
                { 67.913326469405, 201.435784493757, 156.365797422477 },
                { 241.492757985878, 193.318654606482, 120.5791435498 },
                { 187, 39, 36 },
                { 165.729286754482, 163.820392363662, 228.224842303106 },
                { 229.136746492872, 194.131506079963, 233.688057454825 },
                { 226.464699878818, 143.660482209742, 224.56605401776 },
                { 47, 43, 27 },
                { 87.5563448435715, 206.046916072375, 86.9893029678825 },
                { 140.884544460222, 62.4018669461537, 72.2334115190627 },
                { 77.6163845348157, 154.276412827855, 87.4082836188324 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Mean>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][3] =
            {
                { 0, 0, 0 },
                { 151357.879015652, 53205.9839771735, 126078.50459435 },
                { 7027068.56229722, 395452.71116911, 1423540.62712511 },
                { 1095784.88644502, 378691.884028378, 1775553.3384559 },
                { 0, 0, 0 },
                { 385703.899603018, 426965.92005703, 26580.7561779026 },
                { 284311.862464607, 737198.877904507, 925147.288415809 },
                { 686241.437583964, 1421453.13838044, 2126896.11300357 },
                { 0, 0, 0 },
                { 1817091.39576743, 148856.520141535, 1785720.47479914 },
                { 814282.127559123, 2132526.15661061, 2264676.62617193 },
                { 1586333.78340258, 526576.835171952, 1421048.36393696 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Moment2>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][3] =
            {
                { 0, 0, 0 },
                { 151357.879015652, 53205.9839771733, 126078.50459435 },
                { 7027068.56229725, 395452.711169113, 1423540.6271251 },
                { 1095784.88644502, 378691.884028384, 1775553.33845591 },
                { 0, 0, 0 },
                { 385703.899603019, 426965.920057029, 26580.7561779027 },
                { 284311.86246461, 737198.877904507, 925147.288415813 },
                { 686241.43758397, 1421453.13838044, 2126896.11300356 },
                { 0, 0, 0 },
                { 1817091.39576743, 148856.520141536, 1785720.47479914 },
                { 814282.127559124, 2132526.1566106, 2264676.62617193 },
                { 1586333.78340258, 526576.83517195, 1421048.36393696 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Moment2_2Pass>(test_aors[(k)])).begin(), 5e-15);
        }

        {
            double x_tab[][3] =
            {
                { 0, 0, 0 },
                { 324.106807314029, 113.931443205939, 269.975384570342 },
                { 1490.04846528779, 83.8534162784373, 301.853398457402 },
                { 531.67631559681, 183.741816607656, 861.500892021301 },
                { 0, 0, 0 },
                { 785.547657032624, 869.584358568289, 54.1359596291296 },
                { 121.969910967227, 316.258634879668, 396.888583618966 },
                { 271.134507144988, 561.617202046794, 840.338250890387 },
                { 0, 0, 0 },
                { 939.06532081004, 76.9284341816721, 922.852958552526 },
                { 673.517061670077, 1763.8760600584, 1873.18165936471 },
                { 1219.31881891052, 404.747759548003, 1092.2739154012 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Variance>(const_cast<const accumulators_v2d &>(test_aors[(k)]))).begin(), 5e-15);
        }

        {
            double x_tab[][3] =
            {
                { 0, 0, 0 },
                { 324.106807314029, 113.931443205939, 269.975384570342 },
                { 1490.04846528779, 83.8534162784373, 301.853398457402 },
                { 531.67631559681, 183.741816607656, 861.500892021301 },
                { 0, 0, 0 },
                { 785.547657032624, 869.584358568289, 54.1359596291296 },
                { 121.969910967227, 316.258634879668, 396.888583618966 },
                { 271.134507144988, 561.617202046794, 840.338250890387 },
                { 0, 0, 0 },
                { 939.06532081004, 76.9284341816721, 922.852958552526 },
                { 673.517061670077, 1763.8760600584, 1873.18165936471 },
                { 1219.31881891052, 404.747759548003, 1092.2739154012 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Variance>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][3] =
            {
                { 0, -0, -0 },
                { -0.0471999476074436, 0.0225100820944068, 0.0237907216783446 },
                { 0.0203748046142414, 0.0194368293242087, 0.0187743028374562 },
                { -0.0568030216004321, 0.0393764466733232, 0.0385930655364708 },
                { 0, -0, -0 },
                { 0.0129160274209086, 0.0125150661442668, -0.0118529891497778 },
                { -0.0941694199004041, -0.0546479844548244, -0.097100401212693 },
                { -0.0872381883732133, -0.0259197775770145, -0.0883967204365516 },
                { 0, -0, -0 },
                { 0.0251091853946739, 0.0256996935947648, 0.0254673230130096 },
                { 0.0198939243150143, 0.0251674700870824, 0.0188723958474796 },
                { 0.0145664507906367, 0.0122393102871076, 0.0124156595502772 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Skewness>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][3] =
            {
                { 0, -0, -0 },
                { 0.0130922311877067, 0.00542906726593727, 0.0057287568787488 },
                { 0.000995372849455257, 0.000992068018240971, 0.00108301695110664 },
                { 0.0043777958589122, 0.00286214787763822, 0.00284567786675031 },
                { 0, -0, -0 },
                { 0.00473657112866356, 0.00471809290711345, 0.0116046931055875 },
                { 0.0114994795412696, 0.00647975768808986, 0.0113138535470433 },
                { 0.00974007485620935, 0.00343407472637769, 0.0098084392352915 },
                { 0, -0, -0 },
                { 0.00152394272486234, 0.00160790709278664, 0.00154341215192386 },
                { 0.00211927272212755, 0.00240580104758013, 0.00206332315300248 },
                { 0.00201510610202909, 0.00194267614266703, 0.00195582295690013 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Kurtosis>(test_aors[(k)])).begin(), 5e-15);
        }
         
        {
            double x_tab[] =
            {
                44251,
                467,
                4716,
                2061,
                1321,
                491,
                2331,
                2531,
                1386,
                1935,
                1209,
                1301
            };
            double x = x_tab[(k)];
            shouldEqual(nan_squash(get<acc::Count>(test_aors[(k)])), x);
        }

        {
            double x_tab[][3] =
            {
                { 0, 0, 0 },
                { 324.106807314029, 113.931443205939, 269.975384570342 },
                { 1490.04846528779, 83.8534162784373, 301.853398457402 },
                { 531.67631559681, 183.741816607656, 861.500892021301 },
                { 0, 0, 0 },
                { 785.547657032624, 869.584358568289, 54.1359596291296 },
                { 121.969910967227, 316.258634879668, 396.888583618966 },
                { 271.134507144988, 561.617202046794, 840.338250890387 },
                { 0, 0, 0 },
                { 939.06532081004, 76.9284341816721, 922.852958552526 },
                { 673.517061670077, 1763.8760600584, 1873.18165936471 },
                { 1219.31881891052, 404.747759548003, 1092.2739154012 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Variance>(test_aors[(k)])).begin(), 5e-15);
        }

        {
            double x_tab[][3] =
            {
                { 251, 253, 250 },
                { 90.6309234812139, 208.945350778528, 170.528369020468 },
                { 36, 195, 145 },
                { 113.839833780299, 210.580726983075, 109.144489317492 },
                { 187, 39, 36 },
                { 118.593487117183, 114.212450351746, 216.499556981456 },
                { 134.229093937077, 53.3957699485034, 60.9830225692155 },
                { 110.446729613724, 16.4310826625539, 20.8175998381671 },
                { 47, 43, 27 },
                { 62.5690639680764, 199.390452546177, 63.5103920170893 },
                { 109.854720662883, 16.2678275526212, 19.7875100817635 },
                { 22.9135458480419, 122.145619928757, 34.6354994660794 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Min>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][3] =
            {
                { 251, 253, 250 },
                { 231.339611982162, 243.60624313095, 234.452271513179 },
                { 242.654567477069, 246.121415567509, 247.323393948054 },
                { 251.997664863313, 183.163459568119, 97.3572759131751 },
                { 187, 39, 36 },
                { 238.105293397333, 240.32538730135, 246.48304522072 },
                { 241.08146673097, 205.096154262385, 242.043304857592 },
                { 239.323715428994, 190.259623102833, 239.897336643791 },
                { 47, 43, 27 },
                { 189.353806582618, 235.869789293751, 188.672697989038 },
                { 220.890291212331, 190.023294397645, 208.934745452046 },
                { 197.472899093947, 222.2912249782, 199.458474443288 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::Max>(test_aors[(k)])).begin(), 5e-15);
        }
        
        {
            double x_tab[][3] =
            {
                { 0, 0, 0 },
                { 324.802315484231, 114.175931281488, 270.554730889163 },
                { 1490.36448829209, 83.8712006721337, 301.917418266195 },
                { 531.934410895643, 183.831011664261, 861.919096337816 },
                { 0, 0, 0 },
                { 787.150815516364, 871.359020524551, 54.2464411793931 },
                { 122.022258568501, 316.394368199359, 397.058922066871 },
                { 271.241674934373, 561.839185130607, 840.670400396668 },
                { 0, 0, 0 },
                { 939.550876818731, 76.9682110349201, 923.330131747227 },
                { 674.074608906559, 1765.33622235977, 1874.73230643372 },
                { 1220.25675646352, 405.059103978424, 1093.11412610535 }
            };
            TinyVector<double, 3> x;
            x = tab_assign<double, 3>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::UnbiasedVariance>(test_aors[(k)])).begin(), 5e-15);
        }

        {
            double x_tab[][2] =
            {
                { 157.017875302254, 101.195746988769 },
                { 53.0085653104925, 20.1006423982869 },
                { 254.780534351145, 83.0949957591179 },
                { 221.383794274624, 45.8515283842795 },
                { 71.5972747918243, 42.551854655564 },
                { 154.187372708758, 41.1059063136456 },
                { 74.5285285285287, 95.2711282711283 },
                { 91.3512445673648, 110.925325958119 },
                { 146.313131313132, 106.463924963925 },
                { 190.585012919896, 149.702842377261 },
                { 59.1927212572374, 147.971877584781 },
                { 275.843197540353, 152.858570330515 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordMean>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 9231.25355044262, 4028.66535083009 },
                { 3013.64874890526, 358.557322927795 },
                { 1211.00124110548, 582.493096246803 },
                { 341.711508359424, 139.698480196793 },
                { 256.401022096329, 338.642466280162 },
                { 151.036174563735, 119.357419290612 },
                { 321.123488997364, 395.902465436571 },
                { 367.54671812707, 369.70244433271 },
                { 57.6912705700585, 461.225610568686 },
                { 138.877398126448, 179.542188303319 },
                { 136.759384571599, 105.327579683935 },
                { 185.465036355286, 217.52188850178 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordVariance>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 9231.46216634207, 4028.75639411485 },
                { 3020.11580630635, 359.326759243091 },
                { 1211.25808124145, 582.616636670185 },
                { 341.877387732414, 139.766294993005 },
                { 256.595265294887, 338.899013603101 },
                { 151.344411654682, 119.601005860593 },
                { 321.261310237277, 396.072380657788 },
                { 367.691993509729, 369.848571781063 },
                { 57.7329249170405, 461.558625449963 },
                { 138.949206501901, 179.635022940497 },
                { 136.872595982669, 105.414771388971 },
                { 185.607701767867, 217.689213031397 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordUnbiasedVariance>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 9.35529626378307e-05, -0.000201033506253422 },
                { 0.148321279819828, 0.168281345855187 },
                { -0.047472423299276, 0.00356103018462776 },
                { 0.000192179578835091, -0.00750899734840381 },
                { 0.0447830864485579, 0.0337109734387289 },
                { -0.0464672218640941, 0.103707338527541 },
                { 0.0210872294593684, 0.0238223404097032 },
                { -0.0103137773071384, -0.00733595713606302 },
                { 0.0538155226222176, 0.0117063181947805 },
                { 0.000516439404180237, 0.000189488750831951 },
                { -8.63052630006994e-05, -0.00556197725328459 },
                { -0.00551260750647344, 0.0104159760632156 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordSkewness>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 4.13914053003152e-05, 3.62621615042355e-05 },
                { 0.0267309430128425, 0.0327149072314834 },
                { 0.00458178237455199, 0.000650670524587656 },
                { 0.000872110163342214, 0.00223852449558276 },
                { 0.0119946290399348, 0.00609280241121389 },
                { 0.0187822434433972, 0.0752720469824458 },
                { 0.00329468474121094, 0.00224127819493772 },
                { 0.000890674282739669, 0.000934912350687164 },
                { 0.0094492246661146, 0.00179198734443597 },
                { 0.0010008280421884, 0.00101530261174722 },
                { 0.00168461297160918, 0.00164348097662748 },
                { 0.00148495992856451, 0.00186664162262748 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordKurtosis>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 6948198, 4478013 },
                { 24755, 9387 },
                { 1201545, 391876 },
                { 456272, 94500 },
                { 94580, 56211 },
                { 75706, 20183 },
                { 173726, 222077 },
                { 231210, 280752 },
                { 202790, 147559 },
                { 368782, 289675 },
                { 71564, 178898 },
                { 358872, 198869 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordSum>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 408492200.860636, 178272470.439582 },
                { 1407373.96573876, 167446.269807281 },
                { 5711081.85305343, 2747037.44189992 },
                { 704267.418728772, 287918.567685591 },
                { 338705.750189251, 447346.697956093 },
                { 74158.761710794, 58604.4928716905 },
                { 748538.852852855, 922848.646932647 },
                { 930260.743579614, 935716.88660609 },
                { 79960.101010101, 639258.696248199 },
                { 268727.765374677, 347414.134366921 },
                { 165342.095947064, 127341.043837877 },
                { 241290.012298227, 282995.976940816 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordMoment2>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 408492200.860634, 178272470.439586 },
                { 1407373.96573876, 167446.269807281 },
                { 5711081.85305344, 2747037.44189991 },
                { 704267.418728773, 287918.567685588 },
                { 338705.750189251, 447346.697956094 },
                { 74158.7617107943, 58604.4928716903 },
                { 748538.852852854, 922848.646932648 },
                { 930260.74357961, 935716.886606083 },
                { 79960.1010101013, 639258.696248195 },
                { 268727.765374677, 347414.134366924 },
                { 165342.095947064, 127341.043837884 },
                { 241290.012298232, 282995.976940816 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordMoment2_2Pass>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 0, 0 },
                { 6, 18 },
                { 4, 19 },
                { 190, 20 },
                { 40, 77 },
                { 58, 166 },
                { 35, 153 },
                { 36, 153 },
                { 133, 96 },
                { 168, 138 },
                { 36, 131 },
                { 249, 139 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordMin>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 319, 199 },
                { 297, 91 },
                { 303, 51 },
                { 253, 55 },
                { 203, 132 },
                { 174, 45 },
                { 200, 174 },
                { 121, 109 },
                { 196, 130 },
                { 213, 147 },
                { 82, 164 },
                { 302, 130 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::CoordMax>(test_aors[(k)])).begin(), 5e-15);
        }

        {
            double x_tab[][2] =
            {
                { 68354.1932934027, 44053.287798127 },
                { 17650.5390022635, 6733.96508211027 },
                { 67452.4810544542, 22014.5881686907 },
                { 73645.4646687636, 15026.2391354058 },
                { 13917.5231439103, 8271.49390402137 },
                { 50325.7227311568, 13260.1582318924 },
                { 28441.5925476105, 35945.9121540738 },
                { 32279.7100861791, 38548.874140171 },
                { 10123.1347700802, 7366.04193272973 },
                { 46154.185046169, 36248.804902978 },
                { 10116.8921831957, 25714.188380754 },
                { 53814.890913385, 29995.8838902959 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedMean>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 1749414860.3444, 763472370.635817 },
                { 326816056.977875, 38873811.4838026 },
                { 122852815.587156, 47126378.2650026 },
                { 52514126.5941582, 10359137.818378 },
                { 9688369.02093186, 12795944.2308622 },
                { 33663807.1910754, 11384575.1296847 },
                { 50276374.7822094, 36209519.9562075 },
                { 66048655.9780736, 43693265.4591786 },
                { 276168.11221887, 2207886.99779231 },
                { 45253778.3866511, 33989613.8109884 },
                { 12813910.8513701, 73219896.5950605 },
                { 157016412.855393, 67947868.8844765 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedVariance>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 1749454395.1435, 763489624.248713 },
                { 327517378.988557, 38957231.6801198 },
                { 122878871.327472, 47136373.2550907 },
                { 52539618.8886214, 10364166.5260568 },
                { 9695708.69443257, 12805638.1280068 },
                { 33732508.8384041, 11407808.95648 },
                { 50297952.6254635, 36225060.5227123 },
                { 66074762.1662073, 43710535.5245775 },
                { 276367.511577873, 2209481.14002899 },
                { 45277177.4447621, 34007188.5854511 },
                { 12824518.3934656, 73280509.0922419 },
                { 157137194.711436, 68000136.4759261 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedUnbiasedVariance>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 9.35529626380483e-05, -0.000201033506253035 },
                { 0.148551326440454, 0.167819316544869 },
                { -0.012172839505341, 0.0101825665227328 },
                { 0.0151071780592265, -0.0143354144765889 },
                { 0.0447830864485582, 0.0337109734387289 },
                { -0.00366248808457122, 0.144933259013347 },
                { 0.00710677088608147, 0.00910741980275403 },
                { -0.0150819941773359, -0.00726164736997879 },
                { 0.053815522622225, 0.0117063181947808 },
                { 0.0290359694175935, 0.0300278064121328 },
                { 0.047977179924778, 0.028284213762116 },
                { 0.024983386216469, 0.0285154514435182 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedSkewness>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 4.13914053003146e-05, 3.62621615042379e-05 },
                { 0.0270052352269365, 0.0329491265476083 },
                { 0.00221650780929717, 0.000748318147749532 },
                { 0.00158060502569846, 0.00226865403942461 },
                { 0.0119946290399349, 0.00609280241121388 },
                { 0.00715408561564861, 0.0893013218488034 },
                { 0.00179024990145288, 0.00139340430671746 },
                { 0.00118828742649687, 0.00105762916711132 },
                { 0.0094492246661151, 0.00179198734443594 },
                { 0.00209607572663346, 0.00219066599302646 },
                { 0.00519164918415026, 0.00267080136228542 },
                { 0.00286334295164008, 0.00305355056574576 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedKurtosis>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 3024741407.42644, 1949402038.35477 },
                { 8242801.71405707, 3144761.6933455 },
                { 318105900.652807, 103820797.803545 },
                { 151783302.682322, 30969078.8580714 },
                { 18385048.0731055, 10926643.4472122 },
                { 24709929.860998, 6510737.69185915 },
                { 66297352.2284803, 83789921.2311459 },
                { 81699946.2281191, 97567200.4487727 },
                { 14030664.7913311, 10209334.1187635 },
                { 89308348.064337, 70141437.4872623 },
                { 12231322.6494836, 31088453.7523315 },
                { 70013173.0783138, 39024644.941275 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedSum>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 77413356985099.9, 33784415873005.6 },
                { 152623098608.668, 18154069962.9358 },
                { 579373878309.029, 222247999897.752 },
                { 108231614910.56, 21350183043.6771 },
                { 12798335476.651, 16903442328.9689 },
                { 16528929330.818, 5589826388.67521 },
                { 117194229617.33, 84404391017.9197 },
                { 167169148280.504, 110587654877.181 },
                { 382769003.535354, 3060131378.94015 },
                { 87566061178.1699, 65769902724.2625 },
                { 15492018219.3064, 88522854983.4282 },
                { 204278353124.867, 88400177418.7039 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedMoment2>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 77413356985098.7, 33784415873006.5 },
                { 152623098608.668, 18154069962.9358 },
                { 579373878309.03, 222247999897.754 },
                { 108231614910.561, 21350183043.6769 },
                { 12798335476.651, 16903442328.969 },
                { 16528929330.818, 5589826388.67521 },
                { 117194229617.33, 84404391017.9201 },
                { 167169148280.504, 110587654877.182 },
                { 382769003.535355, 3060131378.94014 },
                { 87566061178.1692, 65769902724.2626 },
                { 15492018219.3064, 88522854983.428 },
                { 204278353124.867, 88400177418.7039 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedMoment2_2Pass>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 0, 0 },
                { 2458.02275696471, 7374.06827089413 },
                { 1674.77292261546, 7955.17138242342 },
                { 53670.3296741816, 34727.8603774116 },
                { 7775.44854011651, 14967.7384397243 },
                { 18768.6802019356, 53717.2571296779 },
                { 7683.39483334837, 22893.3805238543 },
                { 6394.48642125357, 16671.3395982682 },
                { 9202.02385347919, 6642.06233033084 },
                { 38093.4371650261, 32682.4375677212 },
                { 6197.52778034227, 16604.6970718604 },
                { 35865.4488383351, 20900.0097547133 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedMin>(test_aors[(k)])).begin(), 5e-15);
        }
        {
            double x_tab[][2] =
            {
                { 138869.460681606, 86630.1651273966 },
                { 101431.63654982, 31078.3802223355 },
                { 107416.762486961, 15598.4737604828 },
                { 101258.478967449, 9605.54741193195 },
                { 39460.4013410913, 25658.9801823845 },
                { 62806.1516333851, 18622.7542633875 },
                { 64973.710444769, 56527.128086949 },
                { 46941.1768383563, 36466.6993620289 },
                { 13560.8772577588, 8994.45940565635 },
                { 70727.0854784769, 54379.7874197723 },
                { 28869.8161401725, 57739.6322803449 },
                { 108128.069880807, 45471.0757445777 }
            };
            TinyVector<double, 2> x;
            x = tab_assign<double, 2>(x, x_tab[(k)]);
            shouldEqualSequenceTolerance(x.begin(), x.end(), nan_squash(get<acc::WeightedMax>(test_aors[(k)])).begin(), 5e-15);
        }
        
    }
}

void test6()
{
    const unsigned m_size = 1200;
    unsigned m_test[m_size];
    *m_test = 2;
    for (unsigned i = 1; i != m_size; ++i)
    {
        for (unsigned n = 1 + m_test[i - 1];; ++n)
        {
            unsigned j = 0;
            for (; j != i; ++j)
                if (n % m_test[j] == 0)
                    break;
            if (j == i)
            {
                m_test[i] = n;
                break;
            }
        }
    }
    
    typedef use_template_list<Accumulators, double,
                              acc::Count,
                              acc::Sum,
                              acc::Mean,
                              acc::Moment2,
                              acc::Moment2_2Pass,
                              acc::Moment3,
                              acc::Moment4,
                              acc::Variance,
                              acc::UnbiasedVariance,
                              acc::Skewness,
                              acc::Kurtosis,
                              acc::Min,
                              acc::Max
                             >
        m_test_accumulators;

    m_test_accumulators mt_0;
    m_test_accumulators mt_1;
    m_test_accumulators mt_2;
    m_test_accumulators mt_a;
    m_test_accumulators mt_b;
    m_test_accumulators mt_c;
    
    inspectSequence(m_test, m_test + m_size, mt_0);

    inspectSequence(m_test, m_test + m_size / 2, mt_1);
    inspectSequence(m_test + m_size / 2, m_test + m_size, mt_2);
    mt_1(mt_2); // merge 2 in 1

    inspectSequence(m_test, m_test + m_size / 6, mt_a);
    inspectSequence(m_test + m_size / 6, m_test + m_size / 5 * 4, mt_b);
    inspectSequence(m_test + m_size / 5 * 4, m_test + m_size, mt_c);
    // merge 3 in 1
    mt_b(mt_a);
    mt_a(mt_b, mt_c);

    shouldEqualTolerance(get<acc::Count>(mt_0), get<acc::Count>(mt_1), 3e-15);
    shouldEqualTolerance(get<acc::Count>(mt_0), get<acc::Count>(mt_a), 3e-15);

    shouldEqualTolerance(get<acc::Sum>(mt_0), get<acc::Sum>(mt_1), 3e-15);
    shouldEqualTolerance(get<acc::Sum>(mt_0), get<acc::Sum>(mt_a), 3e-15);

    shouldEqualTolerance(get<acc::Mean>(mt_0), get<acc::Mean>(mt_1), 3e-15);
    shouldEqualTolerance(get<acc::Mean>(mt_0), get<acc::Mean>(mt_a), 3e-15);

    shouldEqualTolerance(get<acc::Moment2>(mt_0), get<acc::Moment2>(mt_1), 3e-15);
    shouldEqualTolerance(get<acc::Moment2>(mt_0), get<acc::Moment2>(mt_a), 3e-15);

    shouldEqualTolerance(get<acc::Moment2_2Pass>(mt_0), get<acc::Moment2_2Pass>(mt_1), 3e-15);
    shouldEqualTolerance(get<acc::Moment2_2Pass>(mt_0), get<acc::Moment2_2Pass>(mt_a), 3e-15);

    shouldEqualTolerance(get<acc::Moment3>(mt_0), get<acc::Moment3>(mt_1), 5e-14);
    shouldEqualTolerance(get<acc::Moment3>(mt_0), get<acc::Moment3>(mt_a), 5e-14);

    shouldEqualTolerance(get<acc::Moment4>(mt_0), get<acc::Moment4>(mt_1), 5e-14);
    shouldEqualTolerance(get<acc::Moment4>(mt_0), get<acc::Moment4>(mt_a), 5e-14);

    shouldEqualTolerance(get<acc::Variance>(mt_0), get<acc::Variance>(mt_1), 3e-15);
    shouldEqualTolerance(get<acc::Variance>(mt_0), get<acc::Variance>(mt_a), 3e-15);

    shouldEqualTolerance(get<acc::UnbiasedVariance>(mt_0), get<acc::UnbiasedVariance>(mt_1), 3e-15);
    shouldEqualTolerance(get<acc::UnbiasedVariance>(mt_0), get<acc::UnbiasedVariance>(mt_a), 3e-15);

    shouldEqualTolerance(get<acc::Skewness>(mt_0), get<acc::Skewness>(mt_1), 5e-14);
    shouldEqualTolerance(get<acc::Skewness>(mt_0), get<acc::Skewness>(mt_a), 5e-14);

    shouldEqualTolerance(get<acc::Kurtosis>(mt_0), get<acc::Kurtosis>(mt_1), 5e-14);
    shouldEqualTolerance(get<acc::Kurtosis>(mt_0), get<acc::Kurtosis>(mt_a), 5e-14);

    shouldEqualTolerance(get<acc::Min>(mt_0), get<acc::Min>(mt_1), 3e-15);
    shouldEqualTolerance(get<acc::Min>(mt_0), get<acc::Min>(mt_a), 3e-15);

    shouldEqualTolerance(get<acc::Max>(mt_0), get<acc::Max>(mt_1), 3e-15);
    shouldEqualTolerance(get<acc::Max>(mt_0), get<acc::Max>(mt_a), 3e-15);
}


struct ObjectFeaturesTest1
{
    void run1() { test1(); }
    void run2() { test2(); }
    void run3() { test3(); }
    void run4() { test4(); }
    void run5() { test5(); }
    void run6() { test6(); }
};

struct HistogramTest
{
    ArrayVector<double> data;
    
    HistogramTest()
    {
        double d[] = { 
            5.73,  5.45,  0.05,  8.52,  9.04,  4.29,  4.68,  9.65,  5.71,
            6.11,  1.98,  9.13,  0.78,  6.82,  0.41,  9.77,  7.33,  2.89,
            3.27,  4.85,  5.01,  1.86,  0.2 ,  0.62,  0.91,  5.19,  4.7 ,
            5.63,  0.2 ,  1.78,  8.27,  4.3 ,  7.46,  1.12,  8.09,  8.67,
            9.7 ,  2.94,  8.71,  2.07,  0.57,  6.71,  3.21,  8.23,  4.72,
            5.49,  8.  ,  8.86,  8.59,  2.97,  7.65,  6.46,  4.13,  3.64,
            8.74,  0.24,  9.55,  2.83,  1.69,  0.7 ,  2.27,  7.8 ,  3.04,
            6.84,  1.15,  4.29,  6.5 ,  8.76,  8.17,  4.7 ,  7.76,  9.54,
            4.58,  4.73,  2.37,  5.24,  1.23,  7.33,  4.54,  3.72,  2.41,
            3.24,  6.05,  5.25,  7.09,  0.26,  8.8 ,  8.47,  2.39,  8.35,
            5.57,  2.45,  4.83,  1.41,  1.46,  4.6 ,  3.34,  8.66,  3.5 ,  1.33
        };
       
        data.insert(data.begin(), d, d + sizeof(d) / sizeof(double));
    }
    
    void testHistogram()
    {
        const int binCount = 10;
        Histogram<double, int> hist(0.0, 10.0, binCount);
        for(unsigned int k=0; k < data.size(); ++k)
            hist.add(data[k]);
        
        int ref[binCount] = { 11, 10, 10,  8, 14, 10,  7,  7, 16,  7 };
        shouldEqualSequence(ref, ref + binCount, hist.data_.begin());

        hist.reset();
        shouldEqual(hist.data_, ArrayVector<double>(binCount));
        
        MultiArray<2, double> m(Shape2(1, binCount));
        HistogramView<double, double> hv(0.0, 10.0, binCount, &m(0), m.stride(1));

        for(unsigned int k=0; k < data.size(); ++k)
            hv.add(data[k]);
        
        shouldEqualSequence(ref, ref + binCount, &m(0));
    }
    
    void testKernelHistogram()
    {
        TrapezoidKernel<double> kernel;
        double kernelRef[] = { 
            0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375,
            0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
            0.4375, 0.375, 0.3125, 0.25, 0.1875, 0.125, 0.0625, 0.0 };
            
        int c = 0;
        for(double k = -kernel.radius(); k <= kernel.radius(); k+=0.125, ++c)
            shouldEqual(kernel[k], kernelRef[c]);
 
        const int binCount = 10;
        ArrayVector<double> hdata(binCount);
        KernelHistogramView<double, TrapezoidKernel<double> > hist(0.0, 10.0, binCount, &hdata[0]);

        hist.add(2.4);
        hist.add(1.8);
        
        //hist.add(6.3);
        hist.add(6.8);
        
        std::cerr << hdata << "\n";
        
        ArrayVector<std::pair<double, double> > modes;
        hist.findModes(&modes);
        std::cerr << modes << ", highest: " << hist.findMode() << "\n";
    }
    
    void testChannelSmoothing()
    {
        ImageImportInfo info("../impex/lenna.xv");
        MultiArray<2, UInt8> in(Shape2(info.width(), info.height())), out(in.shape());
        importImage(info, destImage(in));
        
        MultiArray<2, TinyVector<float, 20> > hdata(in.shape());
        KernelHistogramView<UInt8, TrapezoidKernel<float> > hist(0, 255, 20);
        
        MersenneTwister rand;
        for(int k=0; k<in.size(); ++k)
        {
            if(rand.uniform() < 0.9)
                continue;
            hist.setData(hdata[k].begin());
            hist.add(in[k]);
        }
        
        gaussianSmoothing(srcImageRange(hdata), destImage(hdata), 2.0);
        
        for(int k=0; k<in.size(); ++k)
        {
            hist.setData(hdata[k].begin());
            out[k] = hist.findMode();
        }
        
        exportImage(srcImageRange(out), ImageExportInfo("res.png"));
    }
};

struct FeaturesTestSuite : public vigra::test_suite
{
    FeaturesTestSuite()
        : vigra::test_suite("FeaturesTestSuite")
    {
        add(testCase(&ObjectFeaturesTest1::run1));
        add(testCase(&ObjectFeaturesTest1::run2));
        add(testCase(&ObjectFeaturesTest1::run3));
        add(testCase(&ObjectFeaturesTest1::run4));
        add(testCase(&ObjectFeaturesTest1::run5));
        add(testCase(&ObjectFeaturesTest1::run6));
        
        add(testCase(&HistogramTest::testHistogram));
        add(testCase(&HistogramTest::testKernelHistogram));
        add(testCase(&HistogramTest::testChannelSmoothing));
    }
};

#endif

struct AccumulatorTest
{
    AccumulatorTest()
    {}
    
    void test1()
    {
#if 0
        using namespace vigra::acc1;
        
        typedef Select<int, float>::type Selected;
        
        std::cerr << "bool: " << Contains<Selected, bool>::type::value << "\n";
        std::cerr << "int: " << Contains<Selected, int>::type::value << "\n";
        std::cerr << "float: " << Contains<Selected, float>::type::value << "\n";
        std::cerr << "double: " << Contains<Selected, double>::type::value << "\n";

        typedef AccumulatorList<double, Selected> Selected1;
        
        std::cerr << "bool: " << Contains<Selected1, bool>::type::value << "\n";
        std::cerr << "int: " << Contains<Selected1, int>::type::value << "\n";
        std::cerr << "float: " << Contains<Selected1, float>::type::value << "\n";
        std::cerr << "double: " << Contains<Selected1, double>::type::value << "\n";
#endif
    }
    
    void testScalar()
    {
        using namespace vigra::acc1;
        
        { 
            Accumulator<double, Select<Count> > a;

            shouldEqual(1, a.passesRequired());

            a(1.0);
            a(2.0);
            a(3.0);

            shouldEqual(get<Count>(a), 3.0);

            try 
            {
                get<Mean>(a);
                failTest("get<Mean>() failed to throw exception");
            }
            catch(ContractViolation & c) 
            {
                std::string expected("\nPrecondition violation!\nget(accumulator): attempt to access inactive statistic");
                std::string message(c.what());
                should(0 == expected.compare(message.substr(0,expected.size())));
            }
        }

        { 
            Accumulator<double, Select<Covariance, UnbiasedStdDev, StdDev, Minimum, Maximum, Skewness, Kurtosis> > a;

            shouldEqual(2, a.passesRequired());

            double data[] = { 1.0, 2.0, 3.0, 5.0 };

            for(int k=0; k<4; ++k)
                a(data[k]);

            shouldEqual(get<Count>(a), 4.0);
            shouldEqual(get<Minimum>(a), 1.0);
            shouldEqual(get<Maximum>(a), 5.0);
            shouldEqual(get<Sum>(a), 11.0);
            shouldEqual(get<Mean>(a), 2.75);
            shouldEqualTolerance(get<UnbiasedVariance>(a), 2.9166666666666665, 1e-15);
            shouldEqualTolerance(get<UnbiasedStdDev>(a), sqrt(2.9166666666666665), 1e-15);
            shouldEqualTolerance(get<Variance>(a), 2.1875, 1e-15);
            shouldEqualTolerance(get<StdDev>(a), sqrt(2.1875), 1e-15);
            shouldEqualTolerance(get<Covariance>(a), 2.1875, 1e-15);

            for(int k=0; k<4; ++k)
                a.updatePass2(data[k]);

            shouldEqual(get<Count>(a), 4.0);
            shouldEqualTolerance(get<CentralMoment<2> >(a),  2.1875, 1e-15);
            shouldEqualTolerance(get<Skewness>(a), 0.43465075957466565, 1e-15);
            shouldEqualTolerance(get<Kurtosis>(a), 1.8457142857142856, 1e-15);
        }

        { 
            DynamicAccumulator<double, Select<Covariance, StdDev, Minimum, CentralMoment<2> > > a;
            activate<Count>(a);

            shouldEqual(1, a.passesRequired());

            a(1.0);
            a(2.0);
            a(3.0);

            shouldEqual(get<Count>(a), 3.0);

            try 
            {
                get<Mean>(a);
                failTest("get<Mean>() failed to throw exception");
            }
            catch(ContractViolation & c) 
            {
                std::string expected("\nPrecondition violation!\nget(accumulator): attempt to access inactive statistic");
                std::string message(c.what());
                should(0 == expected.compare(message.substr(0,expected.size())));
            }

            a.reset();
            activate<StdDev>(a);
            activate<Minimum>(a);
            activate<Covariance>(a);
            activate<CentralMoment<2> >(a);

            shouldEqual(2, a.passesRequired());

            a(1.0);
            a(2.0);
            a(3.0);

            shouldEqual(get<Count>(a), 3.0);
            shouldEqual(get<Minimum>(a), 1.0);
            shouldEqual(get<Sum>(a), 6.0);
            shouldEqual(get<Mean>(a), 2.0);
            shouldEqual(get<Variance>(a), 2.0/3.0);
            shouldEqual(get<StdDev>(a), sqrt(2.0/3.0));
            shouldEqual(get<Covariance>(a), 2.0/3.0);

            a.updatePass2(1.0);
            a.updatePass2(2.0);
            a.updatePass2(3.0);

            shouldEqual(get<Count>(a), 3.0);
            shouldEqual(get<CentralMoment<2> >(a), 2.0/3.0);
        }
    }

    void testVector()
    {
        using namespace vigra::acc1;

        {
            typedef TinyVector<int, 3> V;
            typedef Accumulator<V, Select<StdDev, Covariance, Minimum, Maximum, CentralMoment<2> > > A;
            typedef LookupTag<Mean, A>::result_type W;
            typedef LookupTag<Covariance, A>::result_type Var;

            A a;

            a(V(1,2,3));
            a(V(2,3,1));
            a(V(3,1,2));

            a.updatePass2(V(1,2,3));
            a.updatePass2(V(2,3,1));
            a.updatePass2(V(3,1,2));

            shouldEqual(get<Count>(a), 3.0);
            shouldEqual(get<Minimum>(a), V(1));
            shouldEqual(get<Maximum>(a), V(3));
            shouldEqual(get<Sum>(a), W(6.0));
            shouldEqual(get<Mean>(a), W(2.0));
            shouldEqual(get<CentralMoment<2> >(a), W(2.0/3.0));
            shouldEqual(get<Variance>(a), W(2.0/3.0));

            W stddev = sqrt( W(2.0/3.0));
            shouldEqualSequenceTolerance(stddev.data(), stddev.data()+stddev.size(), get<StdDev>(a).data(), 1e-15);

            double covarianceData[] = { 
                2.0/3.0, -1.0/3.0, -1.0/3.0,
               -1.0/3.0,  2.0/3.0, -1.0/3.0,
               -1.0/3.0, -1.0/3.0,  2.0/3.0 };
            Var covariance(3,3, covarianceData);
            shouldEqual(get<Covariance>(a), covariance);
        }

        {
            using namespace vigra::multi_math;
            
            typedef MultiArray<1, int> V;
            typedef TinyVector<int, 3> T;
            typedef Accumulator<V::view_type, Select<Covariance, StdDev, Minimum, Maximum, CentralMoment<2> > > A;
            typedef LookupTag<Mean, A>::result_type W;
            typedef LookupTag<Covariance, A>::result_type Var;

            A a;

            Shape1 s(3);

            V data[] = { V(s, 1), V(s, 2), V(s, 3) };

            a(data[0]);
            a(data[1]);
            a(data[2]);

            a.updatePass2(data[0]);
            a.updatePass2(data[1]);
            a.updatePass2(data[2]);

            shouldEqual(get<Count>(a), 3.0);
            shouldEqual(get<Minimum>(a), V(s, 1));
            shouldEqual(get<Maximum>(a), V(s, 3));
            shouldEqual(get<Sum>(a), W(s, 6.0));
            shouldEqual(get<Mean>(a),  W(s, 2.0));
            shouldEqual(get<CentralMoment<2> >(a),  W(s, 2.0 / 3.0));

            Var covariance(3,3);
            covariance.init(2.0/3.0);
            shouldEqual(get<Covariance>(a), covariance);

            W variance(s, 2.0/3.0);
            shouldEqual(get<Variance>(a), variance);

            W stddev = sqrt(variance);
            shouldEqualSequenceTolerance(stddev.data(), stddev.data()+stddev.size(), get<StdDev>(a).data(), 1e-15);

            a(V(s, T(0, 2, 4).begin()));
            shouldEqual(get<Minimum>(a), V(s, T(0,1,1).begin()));
            shouldEqual(get<Maximum>(a), V(s, T(3,3,4).begin()));
        }
    }

    void testMerge()
    {
        using namespace vigra::acc1;
        
        typedef Accumulator<double, Select<Covariance, StdDev, Minimum, Maximum, CentralMoment<2> > > A;

        A a, b;

        a(1.0);
        a(2.0);
        a(3.0);

        b(4.0);
        b(5.0);

        a.updatePass2(1.0);
        a.updatePass2(2.0);
        a.updatePass2(3.0);

        b.updatePass2(4.0);
        b.updatePass2(5.0);

        a += b;

        shouldEqual(get<Count>(a), 5.0);
        shouldEqual(get<Minimum>(a), 1.0);
        shouldEqual(get<Maximum>(a), 5.0);
        shouldEqual(get<Sum>(a), 15.0);
        shouldEqual(get<Mean>(a), 3.0);
        shouldEqual(get<SSD>(a), 10.0);
        shouldEqual(get<Variance>(a), 2.0);
        shouldEqual(get<Covariance>(a), 2.0);
        shouldEqual(get<StdDev>(a), sqrt(2.0));
        shouldEqual(get<CentralMoment<2> >(a), 2.0);
    }
};

struct FeaturesTestSuite : public vigra::test_suite
{
    FeaturesTestSuite()
        : vigra::test_suite("FeaturesTestSuite")
    {
        add(testCase(&AccumulatorTest::testScalar));
        add(testCase(&AccumulatorTest::testVector));
        add(testCase(&AccumulatorTest::testMerge));
    }
};

int main(int argc, char** argv)
{
    FeaturesTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
