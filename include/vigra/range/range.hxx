#ifndef RANGE_HXX_
#define RANGE_HXX_

#include <vigra/range/zip.hxx>
#include <vigra/range/map.hxx>

namespace vigra
{

template <class Range, class Visitor>
void for_each(Range range, Visitor visitor)
{
    while(!range.empty())
    {
        visitor(range.front());
        range.pop_front();
    }
}
template <class Range>
size_t length(const Range r)
{
    return r.length();
}
template <class Range>
typename Range::value_type at(Range r, size_t i)
{
    r.pop_front(i);
    return r.front();
}


template <class Iterator>
class IteratorRange
{
  private:
    Iterator pos;
    Iterator end;
  public:
    IteratorRange(Iterator begin, Iterator end)
    : pos(begin),
      end(end)
    {}
    
    typedef decltype(*pos) value_type;
    typedef unsigned long size_type;
    
    bool empty() const
    {
        return pos == end;
    }
    size_type length() const
    {
        return end - pos;
    }

    value_type front() const
    {
        return *pos;
    }
    void pop_front()
    {
        ++pos;
    }

    value_type back() const
    {
        return *(end - 1);
    }
    void pop_back()
    {
        --end;
    }

    void pop_front(size_type number)
    {
        pos += number;
    }
    void pop_back(size_type number)
    {
        end -= number;
    }
    void take_front(size_type number)
    {
        end = pos + number;
    }
    void take_back(size_type number)
    {
        pos = end - number;
    }
};

template <class Iterator>
IteratorRange<Iterator> range(Iterator begin, Iterator end)
{
    return IteratorRange<Iterator>(begin, end);
}

}

#endif

