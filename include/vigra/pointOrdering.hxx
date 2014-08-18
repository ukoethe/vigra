/************************************************************************/
/*                                                                      */
/*               Copyright 2011-2014 by                                 */
/*               Ullrich Koethe                                         */
/*               Esteban Pardo                                          */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_POINTORDERING_HXX
#define VIGRA_POINTORDERING_HXX

namespace vigra {

namespace detail {

static const bool ASCENDING = false;
static const bool DESCENDING = true;

/**
 * functor to order Points by x-coordinate
 */
template<bool DESCENDING>
class PointXOrdering {
public:
    inline
    bool operator()(Diff2D const & a, Diff2D const & b) const {
        return (DESCENDING) ? !(a.x < b.x) : (a.x < b.x);
    }
    
    template<typename T>
    bool operator()(TinyVector<T, 2> const &a,
            TinyVector<T, 2> const &b) const {
        return (DESCENDING) ? !(a[0] < b[0]) : (a[0] < b[0]);
    }
};

/**
 * functor to order Points by y-coordinate
 */
template<bool DESCENDING>
class PointYOrdering {
public:
    inline
    bool operator()(Diff2D const & a, Diff2D const & b) const {
        return (DESCENDING) ? !(a.y < b.y) : (a.y < b.y);
    }
    
    template<typename T>
    bool operator()(TinyVector<T, 2> const &a,
            TinyVector<T, 2> const &b) const {
        return (DESCENDING) ? !(a[1] < b[1]) : (a[1] < b[1]);
    }
};

/**
 * functor to order Points by yx-coordinates (y first, x after)
 */
template<bool DESCENDING>
class PointYXOrdering {
public:
    inline
    bool operator()(Diff2D const & a, Diff2D const & b) const {
        if (a.y == b.y)
            return PointXOrdering<DESCENDING>()(a, b);
        else
            return PointYOrdering<DESCENDING>()(a, b);
    }
    
    template<typename T>
    bool operator()(TinyVector<T, 2> const &a,
            TinyVector<T, 2> const &b) const {
        if (a[1] == b[1])
            return PointXOrdering<DESCENDING>()(a, b);
        else
            return PointYOrdering<DESCENDING>()(a, b);
    }
};

/**
 * functor to order Points by xy-coordinates (x first, y after)
 */
template<bool DESCENDING>
class PointXYOrdering {
public:
    inline
    bool operator()(Diff2D const & a, Diff2D const & b) const {
        if (a.x == b.x)
            return PointYOrdering<DESCENDING>()(a, b);
        else
            return PointXOrdering<DESCENDING>()(a, b);
    }
    
    template<typename T>
    bool operator()(TinyVector<T, 2> const &a,
            TinyVector<T, 2> const &b) const {
        if (a[0] == b[0])
            return PointYOrdering<DESCENDING>()(a, b);
        else
            return PointXOrdering<DESCENDING>()(a, b);
    }    
};

} // namespace detail

} // namespace vigra

#endif
