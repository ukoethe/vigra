#ifndef SPLICES_HXX
#define SPLICES_HXX
#include <vigra/multi_array.hxx>
#include <iterator>

namespace vigra
{
/** Idea for a class to use for easy splicing
 *
 * usage (with factory function _spl)
 *
 * // copy every even indexed element to a 5 x 5 
 * // sized matrix
 * Matrix<double> a(10, 10)
 * MultiArrayView<2,double> b(_spl_shp(_spl(0,2,10),
 *                                     _spl(0,2,10)));
 * copy_splice(_spl(0,2,10),_spl(0,2,10), a, b);
 * 
 * it is also possible to supply iterator ranges
 * std::vector<int> indices;
 * indices.push_back(3) (...)
 * 
 * copy_splice(_spl(indices.begin(), indices.end()),
 *             _spl(a.shape(1)),
 *             a, b)
 *
 * if you only have a forward iterator then you must
 * specify the size of the splice with 
 * _spl(set.begin(), set.end(), set.size());
 *
 * ok.. what we actually need is a decent iota iterator
 * or something like xrange but for now it should suffice.
 */
template<class T>
class Splice
{
    int size_;
    T begin_;
    T end_;
public:
    Splice(T  &begin, T &end)
        : size_(std::distance(begin, end)),
          begin_(begin),
          end_(end)
    {}

    int operator[](int index)
    {
        T ii = begin_;
        std::advance(ii, index);
        return *ii;
    }

    int size()
    {
        return size_;
    }
};

template<>
class Splice<int>
{
    int begin_;
    int interval_;
    int end_;
    int size_;
    public:
    Splice(int begin, int end)
    : begin_(begin), 
      interval_(1),
      end_(end),
      size_(end - begin)
    {}

    Splice(int begin, int interval, int end) 
    : begin_(begin), 
      interval_(interval),
      end_(end),
      size_(int(std::floor((double(end) -double(begin))/interval)))
    {}
    
    int operator[](int index)
    {
        int ii = begin_ + index * interval_;
        return ii;
    }

    int size()
    {
        return size_;
    }
};

template<class T>
Splice<T> _spl(T b, T e)
{
    return Splice<T>(b, e);
}
template<class T>
Splice<T> _spl(T b, int size, T e)
{
    return Splice<T>(b, size, e);
}

inline Splice<int> _spl(int size)
{
    return Splice<int>(0, size);
}



template<class T, class G>
inline MultiArrayShape<2>::type _spl_shp(Splice<T> f,
                                                  Splice<G> h)
{
    return MultiArrayShape<2>::type(f.size(), h.size());
}




template<   class R, class F, 
            class T, class C,
            class T2, class C2 >
void copy_splice(  Splice<R> _first,
                   Splice<F> _second,
                   MultiArrayView<2, T, C>  src,
                   MultiArrayView<2, T2, C2> dest)
{
    for(int jj = 0 ; jj < _second.size(); ++jj)
    {
        for(int ii = 0 ; ii < _first.size(); ++ii)
        {
            dest(ii, jj) = src(_first[ii], _second[jj]);
        }
    }
}

};
#endif //SPLICES_HXX
