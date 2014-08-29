#ifndef ZIPPED_RANGE_HXX_
#define ZIPPED_RANGE_HXX_

#include <vigra/error.hxx>

#include <utility>

namespace vigra
{

template <class First, class Second>
class ZippedRange
{
public:
    typedef std::pair<typename First::value_type, typename Second::value_type> value_type;
    typedef typename First::size_type size_type;
private:
    First f;
    Second s;
public:
    ZippedRange(First f, Second s)
      : f(f),
        s(s)
    {
        vigra_precondition(f.length() >= s.length(), "second range has lower length than first range");
    }
    ZippedRange(std::pair<First, Second> p)
      : f(p.first),
        s(p.second)
    {
        vigra_precondition(f.length() >= s.length(), "second range has lower length than first range");
    }

    bool empty() const
    {
        return f.empty();
    }
    size_type length() const
    {
        return f.length();
    }
    value_type front() const
    {
        return value_type(f.front(), s.front());
    }
    void pop_front()
    {
        f.pop_front();
        s.pop_front();
    }
    void pop_front(size_type i)
    {
        f.pop_front(i);
        s.pop_front(i);
    }
    void take_front(size_type i)
    {
        f.take_front(i);
        s.take_front(i);
    }

    value_type back() const
    {
        return value_type(f.back(), s.back());
    }
    void pop_back()
    {
        f.pop_back();
        s.pop_back();
    }
    void pop_back(size_type i)
    {
        f.pop_back(i);
        s.pop_back(i);
    }
    void take_back(size_type i)
    {
        f.take_back(i);
        s.take_back(i);
    }
};

template <class First, class Second>
ZippedRange<First, Second> zip(First f, Second s)
{
    return ZippedRange<First, Second>(std::move(f), std::move(s));
}

}

#endif

