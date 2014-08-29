#ifndef MAP_HXX_
#define MAP_HXX_

#include <utility>

namespace vigra
{

template <class Range, class Functor>
class MappedRange
  : public Range
{
private:
    Functor functor;
public:
    typedef decltype(std::declval<Functor>()(std::declval<typename Range::value_type>())) value_type;

    MappedRange(Range range, Functor functor)
      : Range(std::move(range)),
        functor(std::move(functor))
    {}

    value_type front() const
    {
        return functor(Range::front());
    }
    value_type back() const
    {
        return functor(Range::back());
    }
};

template<class Range, class Functor>
MappedRange<Range, Functor> map(Range range, Functor functor)
{
    return MappedRange<Range, Functor>(std::move(range), std::move(functor));
}

}

#endif

