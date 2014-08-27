#ifndef MULTI_MAP_HXX_
#define MULTI_MAP_HXX_

#include <vigra/range/map.hxx>

namespace vigra
{

//namespace multi_map_detail
//{
//
//template <class Range, class Functor, unsigned int N>
//class MappedMultiRangeImpl
//  : public Range
//{
//private:
//    Functor functor;
//public:
//    typedef MappedMultiRangeImpl<typename Range::value_type, Functor, N - 1> value_type;
//    
//    MappedMultiRangeImpl(Range range, Functor functor)
//      : Range(std::move(range)),
//        functor(std::move(functor))
//    {}
//    using Range::Range;
//
//    value_type front() const
//    {
//        return value_type(Range::front(), functor);
//    }
//    value_type back() const
//    {
//        return value_type(Range::back(), functor);
//    }
//};
//
//template <class Range, class Functor>
//class MappedMultiRangeImpl<Range, Functor, 1>
//  : public MappedRange<Range, Functor>
//{
//public:
//    enum{ Dimension = 1 };
//    //using MappedRange<Range, Functor>::MappedRange;
//    MappedMultiRangeImpl(Range const & range, Functor const & functor)
//    : MappedRange<Range, Functor>(range, functor)
//    {}
//
//};
//
//}

//template <class Range, class Functor>
//using MultiMappedRange = multi_map_detail::MappedMultiRangeImpl<Range, Functor, Range::dimension>;

template <class Range, class Functor, unsigned int N=Range::dimension>
class MultiMappedRange
  : public Range
{
private:
    Functor functor;
public:
    typedef MultiMappedRange<typename Range::value_type, Functor, N - 1> value_type;
    
    MultiMappedRange(Range range, Functor functor)
      : Range(std::move(range)),
        functor(std::move(functor))
    {}

    //using Range::Range;

    template <class T1, class T2, class T3>
    MultiMappedRange(T1 t1, T2 t2, T3 t3)
    : base_type(t1, t2, t3)
    {}


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
class MultiMappedRange<Range, Functor, 1>
  : public MappedRange<Range, Functor>
{
public:
    enum{ Dimension = 1 };
    //using MappedRange<Range, Functor>::MappedRange;
    MultiMappedRange(Range const & range, Functor const & functor)
    : MappedRange<Range, Functor>(range, functor)
    {}

};

template <class Range, class Functor>
MultiMappedRange<Range, Functor> multi_map(Range range, Functor functor)
{
    return MultiMappedRange<Range, Functor>(std::move(range), std::move(functor));
}

}

#endif

