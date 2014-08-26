#ifndef MULTI_MAP_HXX_
#define MULTI_MAP_HXX_

#include <vigra/range/map.hxx>

namespace vigra
{

namespace multi_map_detail
{

template <class Range, class Functor, unsigned int N>
class MappedMultiRangeImpl
  : public Range
{
private:
    Functor functor;
public:
    typedef MappedMultiRangeImpl<typename Range::value_type, Functor, N - 1> value_type;
    
    MappedMultiRangeImpl(Range range, Functor functor)
      : Range(std::move(range)),
        functor(std::move(functor))
    {}
    using Range::Range;

    value_type front() const
    {
        return value_type(Range::front(), functor);
    }
    value_type back() const
    {
        return value_type(Range::back(), functor);
    }
};

template <class Range, class Functor>
class MappedMultiRangeImpl<Range, Functor, 1>
  : public MappedRange<Range, Functor>
{
public:
    enum{ Dimension = 1 };
    using MappedRange<Range, Functor>::MappedRange;
};

}

template <class Range, class Functor>
using MultiMappedRange = multi_map_detail::MappedMultiRangeImpl<Range, Functor, Range::dimension>;

template <class Range, class Functor>
MultiMappedRange<Range, Functor> multi_map(Range range, Functor functor)
{
    return {std::move(range), std::move(functor)};
}

}

#endif

