#ifndef MULTI_RANGE_HXX_
#define MULTI_RANGE_HXX_

#include <vigra/range/range.hxx>
#include <vigra/range/multi_array_range.hxx>

#include <vigra/range/multi_zip.hxx>
#include <vigra/range/multi_map.hxx>
#include <vigra/range/multi_at.hxx>

namespace vigra
{

template <class MultiRange, class Visitor>
void for_each_flattened(MultiRange&& range, Visitor& visitor);

namespace multi_range_detail
{
template <unsigned int N>
struct for_each_flattened
{   
    // for N >= 1
    template <class MultiRange_, class Visitor>
    static void make(MultiRange_ range, Visitor& visitor)
    {
        typedef typename MultiRange_::value_type NextRange;
        for_each(range, [&visitor](NextRange sub_range)
        {
            for_each_flattened<N -1>::make(sub_range, visitor);
        });
    }
};
template <>
struct for_each_flattened<0>
{
    template <class Value, class Visitor>
    static void make(Value& val, Visitor& visitor)
    {
        visitor(val);
    }
};

}

template <class MultiRange, class Visitor>
void for_each_flattened(MultiRange range, Visitor&& visitor)
{
    multi_range_detail::for_each_flattened<MultiRange::dimension>::make(range, visitor);
}

}

#endif

