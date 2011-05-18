/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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

#ifndef VIGRA_INTERPOLATING_ACCESSOR_HXX
#define VIGRA_INTERPOLATING_ACCESSOR_HXX


#include "accessor.hxx"
#include "diff2d.hxx"

namespace vigra {

/** \addtogroup DataAccessors
*/
//@{

/********************************************************/
/*                                                      */
/*            BilinearInterpolatingAccessor             */
/*                                                      */
/********************************************************/

/** \brief Bilinear interpolation at non-integer positions.

    This accessor allows an image be accessed at arbitrary non-integer
    coordinates and performs an bi-linear interpolation to
    obtain a pixel value.
    It uses the given ACCESSOR (which is usually the
    accessor originally associated with the iterator)
    to access data.

    <b>\#include</b> \<vigra/accessor.hxx\>
    Namespace: vigra

    <b> Required Interface:</b>

    \code
    ITERATOR iter;
    ACCESSOR a;
    VALUETYPE destvalue;
    float s;
    int x, y;

    destvalue = s * a(iter, x, y) + s * a(iter, x, y);

    \endcode
*/
template <class ACCESSOR, class VALUETYPE>
class BilinearInterpolatingAccessor
{
  public:
    /** the iterators' pixel type
    */
    typedef VALUETYPE value_type;

    /** init from given accessor
    */
    BilinearInterpolatingAccessor(ACCESSOR a)
    : a_(a)
    {}

    /** Interpolate the data item at a non-integer position.
        Ensure that no outside pixels are accessed if
        (x, y) is near the image border (as long as
        0 <= x <= width-1, 0 <= y <= height-1).
    */
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i, float x, float y) const
    {
        int ix = int(x);
        int iy = int(y);
        float dx = x - ix;
        float dy = y - iy;

        value_type ret;

        // avoid dereferencing the iterator outside its range
        if(dx == 0.0)
        {
            if(dy == 0.0)
            {
                ret = a_(i, Diff2D(ix, iy));
            }
            else
            {
                ret = detail::RequiresExplicitCast<value_type>::cast(
                  (1.0 - dy) * a_(i, Diff2D(ix, iy)) +
                  dy * a_(i, Diff2D(ix, iy + 1)));
            }
        }
        else
        {
            if(dy == 0.0)
            {
                ret = detail::RequiresExplicitCast<value_type>::cast(
                  (1.0 - dx) * a_(i, Diff2D(ix, iy)) +
                  dx * a_(i, Diff2D(ix + 1, iy)));
            }
            else
            {
                ret = detail::RequiresExplicitCast<value_type>::cast(
                  (1.0 - dx) * (1.0 - dy) * a_(i, Diff2D(ix, iy)) +
                  dx * (1.0 - dy) * a_(i, Diff2D(ix + 1, iy)) +
                  (1.0 - dx) * dy * a_(i, Diff2D(ix, iy + 1)) +
                  dx * dy * a_(i, Diff2D(ix + 1, iy + 1)));
            }
        }

        return ret;
    }

    /** Interpolate the data item at a non-integer position.
        This function works as long as 0 <= x < width-1,
        0 <= y < height-1. It is slightly faster than <TT>operator()</TT>.
    */
    template <class ITERATOR>
    value_type unchecked(ITERATOR const & i, float x, float y) const
    {
        int ix = int(x);
        int iy = int(y);
        float dx = x - ix;
        float dy = y - iy;
        return detail::RequiresExplicitCast<value_type>::cast(
               (1.0 - dx) * (1.0 - dy) * a_(i, Diff2D(ix, iy)) +
               dx * (1.0 - dy) * a_(i, Diff2D(ix + 1, iy)) +
               (1.0 - dx) * dy * a_(i, Diff2D(ix, iy + 1)) +
               dx * dy * a_(i, Diff2D(ix + 1, iy + 1)));
    }

  private:
    ACCESSOR a_;
};

//@}

} // namespace vigra

#endif /* VIGRA_INTERPOLATING_ACCESSOR_HXX */
