#ifndef MULTI_RANGE_TRAITS_HXX_
#define MULTI_RANGE_TRAITS_HXX_

#include <vigra/multi_shape.hxx>

namespace vigra
{

namespace multi_range_traits_detail
{

template <class Range, unsigned int N>
struct value_type_impl
{
    typedef typename value_type_impl<typename Range::value_type, N - 1>::type type;
};

template <class ValueType>
struct value_type_impl<ValueType, 0>
{
    typedef ValueType type;
};

}

namespace multi_range_traits
{

template <class Range>
struct value_type
{
    typedef typename multi_range_traits_detail::value_type_impl<Range, Range::dimension>::type type;
};

template <class Range>
struct decayed_value_type
{
    typedef typename std::decay<typename value_type<Range>::type>::type type;
};

template <class Range>
struct multi_size_type
{
    typedef typename MultiArrayShape<Range::dimension>::type type; // TODO
};

}

}

#endif

