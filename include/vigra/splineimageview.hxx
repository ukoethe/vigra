/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#ifndef VIGRA_SPLINEIMAGEVIEW_HXX
#define VIGRA_SPLINEIMAGEVIEW_HXX

#include "vigra/mathutil.hxx"
#include "vigra/recursiveconvolution.hxx"
#include "vigra/splines.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/basicimage.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/tinyvector.hxx"
#include "vigra/fixedpoint.hxx"
#include "vigra/multi_array.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*                    SplineImageView                   */
/*                                                      */
/********************************************************/
/** \brief Create a continuous view onto a discrete image using splines.

    This class is very useful if image values or derivatives at arbitrary
    real-valued coordinates are needed. Access at such coordinates is implemented by 
    interpolating the given discrete image values with a spline of the 
    specified <tt>ORDER</TT>. Continuous derivatives are available up to 
    degree <tt>ORDER-1</TT>. If the requested coordinates are near the image border, 
    reflective boundary conditions are applied. In principle, this class can also be used 
    for image resizing, but here the functions from the <tt>resize...</tt> family are 
    more efficient. 
    
*/
template <int ORDER, class VALUETYPE>
class SplineImageView
{
    typedef typename NumericTraits<VALUETYPE>::RealPromote InternalValue;

  public:
  
        /** The view's value type (return type of access and derivative functions).
        */
    typedef VALUETYPE value_type;
    
        /** The view's size type.
        */
    typedef Size2D size_type;

        /** The view's difference type.
        */
    typedef TinyVector<double, 2> difference_type;

        /** The order of the spline used.
        */
    enum StaticOrder { order = ORDER };
    
        /** The type of the internal image holding the spline coefficients.
        */
    typedef BasicImage<InternalValue> InternalImage;
  
  private:
    typedef typename InternalImage::traverser InternalTraverser;
    typedef typename InternalTraverser::row_iterator InternalRowIterator;
    typedef typename InternalTraverser::column_iterator InternalColumnIterator;
    typedef BSpline<ORDER, double> Spline;
    
    enum { ksize_ = ORDER + 1, kcenter_ = ORDER / 2 };
 
  public:   
        /** Construct SplineImageView for the given image.
        
            If <tt>skipPrefiltering = true</tt> (default: <tt>false</tt>), the recursive
            prefilter of the cardinal spline function is not applied, resulting
            in an approximating (smoothing) rather than interpolating spline. This is
            especially useful if customized prefilters are to be applied.
        */
    template <class SrcIterator, class SrcAccessor>
    SplineImageView(SrcIterator is, SrcIterator iend, SrcAccessor sa, bool skipPrefiltering = false)
    : w_(iend.x - is.x), h_(iend.y - is.y), w1_(w_-1), h1_(h_-1),
      x0_(kcenter_), x1_(w_ - kcenter_ - 2), y0_(kcenter_), y1_(h_ - kcenter_ - 2),
      image_(w_, h_),
      x_(-1.0), y_(-1.0),
      u_(-1.0), v_(-1.0)
    {
        copyImage(srcIterRange(is, iend, sa), destImage(image_));
        if(!skipPrefiltering)
            init();
    }
    
        /** Construct SplineImageView for the given image.
        
            If <tt>skipPrefiltering = true</tt> (default: <tt>false</tt>), the recursive
            prefilter of the cardinal spline function is not applied, resulting
            in an approximating (smoothing) rather than interpolating spline. This is
            especially useful if customized prefilters are to be applied.
        */
    template <class SrcIterator, class SrcAccessor>
    SplineImageView(triple<SrcIterator, SrcIterator, SrcAccessor> s, bool skipPrefiltering = false)
    : w_(s.second.x - s.first.x), h_(s.second.y - s.first.y), w1_(w_-1), h1_(h_-1),
      x0_(kcenter_), x1_(w_ - kcenter_ - 2), y0_(kcenter_), y1_(h_ - kcenter_ - 2),
      image_(w_, h_),
      x_(-1.0), y_(-1.0),
      u_(-1.0), v_(-1.0)
    {
        copyImage(srcIterRange(s.first, s.second, s.third), destImage(image_));
        if(!skipPrefiltering)
            init();
    }
    
        /** Access interpolated function at real-valued coordinate <tt>(x, y)</tt>.
            If <tt>(x, y)</tt> is near the image border or outside the image, the value
            is calculated with reflective boundary conditions. An exception is thrown if the 
            coordinate is outside the first reflection. 
        */
    value_type operator()(double x, double y) const;
    
        /** Access derivative of order <tt>(dx, dy)</tt> at real-valued coordinate <tt>(x, y)</tt>.
            If <tt>(x, y)</tt> is near the image border or outside the image, the value
            is calculated with reflective boundary conditions. An exception is thrown if the 
            coordinate is outside the first reflection. 
        */
    value_type operator()(double x, double y, unsigned int dx, unsigned int dy) const;
    
        /** Access 1st derivative in x-direction at real-valued coordinate <tt>(x, y)</tt>.
            Equivalent to <tt>splineView(x, y, 1, 0)</tt>.
        */
    value_type dx(double x, double y) const
        { return operator()(x, y, 1, 0); }
    
        /** Access 1st derivative in y-direction at real-valued coordinate <tt>(x, y)</tt>.
            Equivalent to <tt>splineView(x, y, 0, 1)</tt>.
        */
    value_type dy(double x, double y) const
        { return operator()(x, y, 0, 1); }
    
        /** Access 2nd derivative in x-direction at real-valued coordinate <tt>(x, y)</tt>.
            Equivalent to <tt>splineView(x, y, 2, 0)</tt>.
        */
    value_type dxx(double x, double y) const
        { return operator()(x, y, 2, 0); }
    
        /** Access mixed 2nd derivative at real-valued coordinate <tt>(x, y)</tt>.
            Equivalent to <tt>splineView(x, y, 1, 1)</tt>.
        */
    value_type dxy(double x, double y) const
        { return operator()(x, y, 1, 1); }
    
        /** Access 2nd derivative in y-direction at real-valued coordinate <tt>(x, y)</tt>.
            Equivalent to <tt>splineView(x, y, 0, 2)</tt>.
        */
    value_type dyy(double x, double y) const
        { return operator()(x, y, 0, 2); }
    
        /** Access 3rd derivative in x-direction at real-valued coordinate <tt>(x, y)</tt>.
            Equivalent to <tt>splineView(x, y, 3, 0)</tt>.
        */
    value_type dx3(double x, double y) const
        { return operator()(x, y, 3, 0); }
    
        /** Access 3rd derivative in y-direction at real-valued coordinate <tt>(x, y)</tt>.
            Equivalent to <tt>splineView(x, y, 0, 3)</tt>.
        */
    value_type dy3(double x, double y) const
        { return operator()(x, y, 0, 3); }
    
        /** Access mixed 3rd derivative dxxy at real-valued coordinate <tt>(x, y)</tt>.
            Equivalent to <tt>splineView(x, y, 2, 1)</tt>.
        */
    value_type dxxy(double x, double y) const
        { return operator()(x, y, 2, 1); }
    
        /** Access mixed 3rd derivative dxyy at real-valued coordinate <tt>(x, y)</tt>.
            Equivalent to <tt>splineView(x, y, 1, 2)</tt>.
        */
    value_type dxyy(double x, double y) const
        { return operator()(x, y, 1, 2); }
        
        /** Access interpolated function at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView(d[0], d[1])</tt>.
        */
    value_type operator()(difference_type const & d) const
        { return operator()(d[0], d[1]); }
        
        /** Access derivative of order <tt>(dx, dy)</tt> at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView(d[0], d[1], dx, dy)</tt>.
        */
    value_type operator()(difference_type const & d, unsigned int dx, unsigned int dy) const
        { return operator()(d[0], d[1], dx, dy); }
        
        /** Access 1st derivative in x-direction at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView.dx(d[0], d[1])</tt>.
        */
    value_type dx(difference_type const & d) const
        { return dx(d[0], d[1]); }
        
        /** Access 1st derivative in y-direction at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView.dy(d[0], d[1])</tt>.
        */
    value_type dy(difference_type const & d) const
        { return dy(d[0], d[1]); }
        
        /** Access 2nd derivative in x-direction at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView.dxx(d[0], d[1])</tt>.
        */
    value_type dxx(difference_type const & d) const
        { return dxx(d[0], d[1]); }
        
        /** Access mixed 2nd derivative at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView.dxy(d[0], d[1])</tt>.
        */
    value_type dxy(difference_type const & d) const
        { return dxy(d[0], d[1]); }
        
        /** Access 2nd derivative in y-direction at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView.dyy(d[0], d[1])</tt>.
        */
    value_type dyy(difference_type const & d) const
        { return dyy(d[0], d[1]); }
        
        /** Access 3rd derivative in x-direction at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView.dx3(d[0], d[1])</tt>.
        */
    value_type dx3(difference_type const & d) const
        { return dx3(d[0], d[1]); }
        
        /** Access 3rd derivative in y-direction at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView.dy3(d[0], d[1])</tt>.
        */
    value_type dy3(difference_type const & d) const
        { return dy3(d[0], d[1]); }
        
        /** Access mixed 3rd derivative dxxy at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView.dxxy(d[0], d[1])</tt>.
        */
    value_type dxxy(difference_type const & d) const
        { return dxxy(d[0], d[1]); }
        
        /** Access mixed 3rd derivative dxyy at real-valued coordinate <tt>d</tt>.
            Equivalent to <tt>splineView.dxyy(d[0], d[1])</tt>.
        */
    value_type dxyy(difference_type const & d) const
        { return dxyy(d[0], d[1]); }
        
        /** Access gradient squared magnitude at real-valued coordinate <tt>(x, y)</tt>.
        */
    value_type g2(double x, double y) const;
        
        /** Access 1st derivative in x-direction of gradient squared magnitude 
            at real-valued coordinate <tt>(x, y)</tt>.
        */
    value_type g2x(double x, double y) const;
        
        /** Access 1st derivative in y-direction of gradient squared magnitude 
            at real-valued coordinate <tt>(x, y)</tt>.
        */
    value_type g2y(double x, double y) const;
        
        /** Access 2nd derivative in x-direction of gradient squared magnitude 
            at real-valued coordinate <tt>(x, y)</tt>.
        */
    value_type g2xx(double x, double y) const;
        
        /** Access mixed 2nd derivative of gradient squared magnitude 
            at real-valued coordinate <tt>(x, y)</tt>.
        */
    value_type g2xy(double x, double y) const;
        
        /** Access 2nd derivative in y-direction of gradient squared magnitude 
            at real-valued coordinate <tt>(x, y)</tt>.
        */
    value_type g2yy(double x, double y) const;
        
        /** Access gradient squared magnitude at real-valued coordinate <tt>d</tt>.
        */
    value_type g2(difference_type const & d) const
        { return g2(d[0], d[1]); }
        
        /** Access 1st derivative in x-direction of gradient squared magnitude 
            at real-valued coordinate <tt>d</tt>.
        */
    value_type g2x(difference_type const & d) const
        { return g2x(d[0], d[1]); }
        
        /** Access 1st derivative in y-direction of gradient squared magnitude 
            at real-valued coordinate <tt>d</tt>.
        */
    value_type g2y(difference_type const & d) const
        { return g2y(d[0], d[1]); }
        
        /** Access 2nd derivative in x-direction of gradient squared magnitude 
            at real-valued coordinate <tt>d</tt>.
        */
    value_type g2xx(difference_type const & d) const
        { return g2xx(d[0], d[1]); }
        
        /** Access mixed 2nd derivative of gradient squared magnitude 
            at real-valued coordinate <tt>d</tt>.
        */
    value_type g2xy(difference_type const & d) const
        { return g2xy(d[0], d[1]); }
        
        /** Access 2nd derivative in y-direction of gradient squared magnitude 
            at real-valued coordinate <tt>d</tt>.
        */
    value_type g2yy(difference_type const & d) const
        { return g2yy(d[0], d[1]); }
    
        /** The width of the image.
            <tt>0 &lt;= x &lt;= width()-1</tt> is required for all access functions.
        */
    unsigned int width() const
        { return w_; }
    
        /** The height of the image.
            <tt>0 &lt;= y &lt;= height()-1</tt> is required for all access functions.
        */
    unsigned int height() const
        { return h_; }
    
        /** The size of the image.
            <tt>0 &lt;= x &lt;= size().x-1</tt> and <tt>0 &lt;= y &lt;= size().y-1</tt> 
            are required for all access functions.
        */
    size_type size() const
        { return size_type(w_, h_); }
        
        /** The internal image holding the spline coefficients.
        */
    InternalImage const & image() const
    {
        return image_;
    }
    
        /** Get the array of polynomial coefficients for the facet containing 
            the point <tt>(x, y)</tt>. The array <tt>res</tt> will be resized to
            dimension <tt>(ORDER+1)x(ORDER+1)</tt>. From these coefficients, the
            value of the interpolated function can be calculated by the following
            algorithm
            
            \code
            SplineImageView<ORDER, float> view(...);
            double x = ..., y = ...;
            double dx, dy;
            
            // calculate the local facet coordinates of x and y
            if(ORDER % 2)
            {
                // odd order => facet coordinates between 0 and 1
                dx = x - floor(x);
                dy = y - floor(y);
            }
            else
            {
                // even order => facet coordinates between -0.5 and 0.5
                dx = x - floor(x + 0.5);
                dy = y - floor(y + 0.5);
            }
            
            BasicImage<float> coefficients;
            view.coefficientArray(x, y, coefficients);
            
            float f_x_y = 0.0;
            for(int ny = 0; ny < ORDER + 1; ++ny)
                for(int nx = 0; nx < ORDER + 1; ++nx)
                    f_x_y += pow(dx, nx) * pow(dy, ny) * coefficients(nx, ny);
                    
            assert(abs(f_x_y - view(x, y)) < 1e-6);
            \endcode
        */
    template <class Array>
    void coefficientArray(double x, double y, Array & res) const;
    
        /** Check if x is in the original image range.
            Equivalent to <tt>0 &lt;= x &lt;= width()-1</tt>.
        */
    bool isInsideX(double x) const
    {
        return x >= 0.0 && x <= width()-1.0;
    }
        
        /** Check if y is in the original image range.
            Equivalent to <tt>0 &lt;= y &lt;= height()-1</tt>.
        */
    bool isInsideY(double y) const
    {
        return y >= 0.0 && y <= height()-1.0;
    }
        
        /** Check if x and y are in the original image range.
            Equivalent to <tt>0 &lt;= x &lt;= width()-1</tt> and <tt>0 &lt;= y &lt;= height()-1</tt>.
        */
    bool isInside(double x, double y) const
    {
        return isInsideX(x) && isInsideY(y);
    }
    
        /** Check if x and y are in the valid range. Points outside the original image range are computed
            by reflcective boundary conditions, but only within the first reflection.
            Equivalent to <tt>-width() + ORDER/2 + 2 &lt; x &lt; 2*width() - ORDER/2 - 2</tt> and 
            <tt>-height() + ORDER/2 + 2 &lt; y &lt; 2*height() - ORDER/2 - 2</tt>.
        */
    bool isValid(double x, double y) const
    {
        return x < w1_ + x1_ && x > -x1_ && y < h1_ + y1_ && y > -y1_;
    }
    
        /** Check whether the points <tt>(x0, y0)</tt> and <tt>(x1, y1)</tt> are in
            the same spline facet. For odd order splines, facets span the range
            <tt>(floor(x), floor(x)+1) x (floor(y), floor(y)+1)</tt> (i.e. we have 
            integer facet borders), whereas even order splines have facet between
            half integer values 
            <tt>(floor(x)-0.5, floor(x)+0.5) x (floor(y)-0.5, floor(y)+0.5)</tt>.
        */
    bool sameFacet(double x0, double y0, double x1, double y1) const
    {
         x0 = VIGRA_CSTD::floor((ORDER % 2) ? x0 : x0 + 0.5);
         y0 = VIGRA_CSTD::floor((ORDER % 2) ? y0 : y0 + 0.5);
         x1 = VIGRA_CSTD::floor((ORDER % 2) ? x1 : x1 + 0.5);
         y1 = VIGRA_CSTD::floor((ORDER % 2) ? y1 : y1 + 0.5);
         return x0 == x1 && y0 == y1;
    }
        
  protected:
  
    void init();
    void calculateIndices(double x, double y) const;
    void coefficients(double t, double * const & c) const;
    void derivCoefficients(double t, unsigned int d, double * const & c) const;
    value_type convolve() const;
  
    unsigned int w_, h_;
    int w1_, h1_;
    double x0_, x1_, y0_, y1_;
    InternalImage image_;
    Spline k_;
    mutable double x_, y_, u_, v_, kx_[ksize_], ky_[ksize_];
    mutable int ix_[ksize_], iy_[ksize_];
};

template <int ORDER, class VALUETYPE>
void SplineImageView<ORDER, VALUETYPE>::init()
{
    ArrayVector<double> const & b = k_.prefilterCoefficients();
    
    for(unsigned int i=0; i<b.size(); ++i)
    {
        recursiveFilterX(srcImageRange(image_), destImage(image_), b[i], BORDER_TREATMENT_REFLECT);
        recursiveFilterY(srcImageRange(image_), destImage(image_), b[i], BORDER_TREATMENT_REFLECT);
    }
}

namespace detail
{

template <int i>
struct SplineImageViewUnrollLoop1
{
    template <class Array>
    static void exec(int c0, Array c)
    {
        SplineImageViewUnrollLoop1<i-1>::exec(c0, c);
        c[i] = c0 + i;
    }
};

template <>
struct SplineImageViewUnrollLoop1<0>
{
    template <class Array>
    static void exec(int c0, Array c)
    {
        c[0] = c0;
    }
};

template <int i, class ValueType>
struct SplineImageViewUnrollLoop2
{
    template <class Array1, class RowIterator, class Array2>
    static ValueType
    exec(Array1 k, RowIterator r, Array2 x)
    {
        return k[i] * r[x[i]] + SplineImageViewUnrollLoop2<i-1, ValueType>::exec(k, r, x);
    }
};

template <class ValueType>
struct SplineImageViewUnrollLoop2<0, ValueType>
{
    template <class Array1, class RowIterator, class Array2>
    static ValueType
    exec(Array1 k, RowIterator r, Array2 x)
    {
        return k[0] * r[x[0]];
    }
};

} // namespace detail

template <int ORDER, class VALUETYPE>
void 
SplineImageView<ORDER, VALUETYPE>::calculateIndices(double x, double y) const
{
    if(x == x_ && y == y_)
        return;   // still in cache
    
    if(x > x0_ && x < x1_ && y > y0_ && y < y1_)
    {
        detail::SplineImageViewUnrollLoop1<ORDER>::exec(
                                (ORDER % 2) ? int(x - kcenter_) : int(x + 0.5 - kcenter_), ix_);
        detail::SplineImageViewUnrollLoop1<ORDER>::exec(
                                (ORDER % 2) ? int(y - kcenter_) : int(y + 0.5 - kcenter_), iy_);

        u_ = x - ix_[kcenter_];
        v_ = y - iy_[kcenter_];
    }
    else
    {
        vigra_precondition(isValid(x,y),
                    "SplineImageView::calculateIndices(): coordinates out of range.");
        
        int xCenter = (ORDER % 2) ?
                      (int)VIGRA_CSTD::floor(x) :
                      (int)VIGRA_CSTD::floor(x + 0.5);
        int yCenter = (ORDER % 2) ?
                      (int)VIGRA_CSTD::floor(y) :
                      (int)VIGRA_CSTD::floor(y + 0.5);
        
        if(x >= x1_)
        {
            for(int i = 0; i < ksize_; ++i)
                ix_[i] = w1_ - vigra::abs(w1_ - xCenter - (i - kcenter_));
        }
        else
        {
            for(int i = 0; i < ksize_; ++i)
                ix_[i] = vigra::abs(xCenter - (kcenter_ - i));
        }
        if(y >= y1_)
        {
            for(int i = 0; i < ksize_; ++i)
                iy_[i] = h1_ - vigra::abs(h1_ - yCenter - (i - kcenter_));
        }
        else
        {
            for(int i = 0; i < ksize_; ++i)
                iy_[i] = vigra::abs(yCenter - (kcenter_ - i));
        }
        u_ = x - xCenter;
        v_ = y - yCenter;
    }
    x_ = x;
    y_ = y;
}

template <int ORDER, class VALUETYPE>
void SplineImageView<ORDER, VALUETYPE>::coefficients(double t, double * const & c) const
{
    t += kcenter_;
    for(int i = 0; i<ksize_; ++i)
        c[i] = k_(t-i);
}

template <int ORDER, class VALUETYPE>
void SplineImageView<ORDER, VALUETYPE>::derivCoefficients(double t, 
                                               unsigned int d, double * const & c) const
{
    t += kcenter_;
    for(int i = 0; i<ksize_; ++i)
        c[i] = k_(t-i, d);
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineImageView<ORDER, VALUETYPE>::convolve() const
{
    InternalValue sum;
    sum = ky_[0]*detail::SplineImageViewUnrollLoop2<ORDER, InternalValue>::exec(kx_, image_.rowBegin(iy_[0]), ix_);
    
    for(int j=1; j<ksize_; ++j)
    {
        sum += ky_[j]*detail::SplineImageViewUnrollLoop2<ORDER, InternalValue>::exec(kx_, image_.rowBegin(iy_[j]), ix_);
    }
    return NumericTraits<VALUETYPE>::fromRealPromote(sum);
}

template <int ORDER, class VALUETYPE>
template <class Array>
void 
SplineImageView<ORDER, VALUETYPE>::coefficientArray(double x, double y, Array & res) const
{
    typename Spline::WeightMatrix & weights = Spline::weights();
    InternalValue tmp[ksize_][ksize_]; 
    
    calculateIndices(x, y);
    for(int j=0; j<ksize_; ++j)
    {
        for(int i=0; i<ksize_; ++i)
        {
            tmp[i][j] = 0.0;
            for(int k=0; k<ksize_; ++k)
            {
                tmp[i][j] += weights[i][k]*image_(ix_[k], iy_[j]);
            }
       }       
    }
    res.resize(ksize_, ksize_);
    for(int j=0; j<ksize_; ++j)
    {
        for(int i=0; i<ksize_; ++i)
        {
            res(i,j) = 0.0;
            for(int k=0; k<ksize_; ++k)
            {
                res(i,j) += weights[j][k]*tmp[i][k];
            }
        }       
    }
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineImageView<ORDER, VALUETYPE>::operator()(double x, double y) const
{
    calculateIndices(x, y);
    coefficients(u_, kx_);
    coefficients(v_, ky_);
    return convolve();
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineImageView<ORDER, VALUETYPE>::operator()(double x, double y,
                                                 unsigned int dx, unsigned int dy) const
{
    calculateIndices(x, y);
    derivCoefficients(u_, dx, kx_);
    derivCoefficients(v_, dy, ky_);
    return convolve();
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineImageView<ORDER, VALUETYPE>::g2(double x, double y) const
{
    return sq(dx(x,y)) + sq(dy(x,y));
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineImageView<ORDER, VALUETYPE>::g2x(double x, double y) const
{
    return 2.0*(dx(x,y) * dxx(x,y) + dy(x,y) * dxy(x,y));
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineImageView<ORDER, VALUETYPE>::g2y(double x, double y) const
{
    return 2.0*(dx(x,y) * dxy(x,y) + dy(x,y) * dyy(x,y));
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineImageView<ORDER, VALUETYPE>::g2xx(double x, double y) const
{
    return 2.0*(sq(dxx(x,y)) + dx(x,y) * dx3(x,y) + sq(dxy(x,y)) + dy(x,y) * dxxy(x,y));
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineImageView<ORDER, VALUETYPE>::g2yy(double x, double y) const
{
    return 2.0*(sq(dxy(x,y)) + dx(x,y) * dxyy(x,y) + sq(dyy(x,y)) + dy(x,y) * dy3(x,y));
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineImageView<ORDER, VALUETYPE>::g2xy(double x, double y) const
{
    return 2.0*(dx(x,y) * dxxy(x,y) + dy(x,y) * dxyy(x,y) + dxy(x,y) * (dxx(x,y) + dyy(x,y)));
}

/********************************************************/
/*                                                      */
/*                    SplineImageView0                  */
/*                                                      */
/********************************************************/
template <class VALUETYPE, class INTERNAL_INDEXER>
class SplineImageView0Base
{
    typedef typename INTERNAL_INDEXER::value_type InternalValue;
  public:
    typedef VALUETYPE value_type;
    typedef Size2D size_type;
    typedef TinyVector<double, 2> difference_type;
    enum StaticOrder { order = 0 };
  
  public:

    SplineImageView0Base(unsigned int w, unsigned int h)
    : w_(w), h_(h)
    {}

    SplineImageView0Base(int w, int h, INTERNAL_INDEXER i)
    : w_(w), h_(h), internalIndexer_(i)
    {}

    template <unsigned IntBits1, unsigned FractionalBits1,
              unsigned IntBits2, unsigned FractionalBits2>
    value_type unchecked(FixedPoint<IntBits1, FractionalBits1> x, 
                         FixedPoint<IntBits2, FractionalBits2> y) const
    {
        return internalIndexer_(round(x), round(y));
    }

    template <unsigned IntBits1, unsigned FractionalBits1,
              unsigned IntBits2, unsigned FractionalBits2>
    value_type unchecked(FixedPoint<IntBits1, FractionalBits1> x, 
                         FixedPoint<IntBits2, FractionalBits2> y, 
                         unsigned int dx, unsigned int dy) const
    {
        if((dx != 0) || (dy != 0))
            return NumericTraits<VALUETYPE>::zero();
        return unchecked(x, y);
    }

    value_type unchecked(double x, double y) const
    {
        return internalIndexer_((int)(x + 0.5), (int)(y + 0.5));
    }

    value_type unchecked(double x, double y, unsigned int dx, unsigned int dy) const
    {
        if((dx != 0) || (dy != 0))
            return NumericTraits<VALUETYPE>::zero();
        return unchecked(x, y);
    }

    value_type operator()(double x, double y) const
    {
        int ix, iy;
        if(x < 0.0)
        {
            ix = (int)(-x + 0.5);
            vigra_precondition(ix <= (int)w_ - 1,
                    "SplineImageView::operator(): coordinates out of range.");
        }
        else
        {
            ix = (int)(x + 0.5);
            if(ix >= (int)w_)
            {
                ix = 2*w_-2-ix;
                vigra_precondition(ix >= 0,
                        "SplineImageView::operator(): coordinates out of range.");
            }
        }
        if(y < 0.0)
        {
            iy = (int)(-y + 0.5);
            vigra_precondition(iy <= (int)h_ - 1,
                    "SplineImageView::operator(): coordinates out of range.");
        }
        else 
        {
            iy = (int)(y + 0.5);
            if(iy >= (int)h_)
            {
                iy = 2*h_-2-iy;
                vigra_precondition(iy >= 0,
                        "SplineImageView::operator(): coordinates out of range.");
            }
        }
        return unchecked(ix, iy);
    }

    value_type operator()(double x, double y, unsigned int dx, unsigned int dy) const
    {
        if((dx != 0) || (dy != 0))
            return NumericTraits<VALUETYPE>::zero();
        return operator()(x, y);
    }

    value_type dx(double x, double y) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dy(double x, double y) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dxx(double x, double y) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dxy(double x, double y) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dyy(double x, double y) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dx3(double x, double y) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dy3(double x, double y) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dxxy(double x, double y) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dxyy(double x, double y) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type operator()(difference_type const & d) const
        { return operator()(d[0], d[1]); }

    value_type operator()(difference_type const & d, unsigned int dx, unsigned int dy) const
        { return operator()(d[0], d[1], dx, dy); }

    value_type dx(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dy(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dxx(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dxy(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dyy(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dx3(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dy3(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dxxy(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type dxyy(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type g2(double x, double y) const;

    value_type g2x(double x, double y) const;

    value_type g2y(double x, double y) const;

    value_type g2xx(double x, double y) const;

    value_type g2xy(double x, double y) const;

    value_type g2yy(double x, double y) const;

    value_type g2(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type g2x(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type g2y(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type g2xx(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type g2xy(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    value_type g2yy(difference_type const & d) const
        { return NumericTraits<VALUETYPE>::zero(); }

    unsigned int width() const
        { return w_; }

    unsigned int height() const
        { return h_; }

    size_type size() const
        { return size_type(w_, h_); }

    template <class Array>
    void coefficientArray(double x, double y, Array & res) const
    {
        res.resize(1, 1);
        res(0, 0) = operator()(x,y);
    }

    bool isInsideX(double x) const
    {
        return x >= 0.0 && x < width() - 0.5; // FIXME
    }

    bool isInsideY(double y) const
    {
        return y >= 0.0 && y < height() - 0.5;
    }

    bool isInside(double x, double y) const
    {
        return isInsideX(x) && isInsideY(y);
    }

    bool isValid(double x, double y) const
    {
        return x < 2*w_ && x > -w_ && y < 2*h_ && y > -h_;
    }

    bool sameFacet(double x0, double y0, double x1, double y1) const
    {
         x0 = VIGRA_CSTD::floor((ORDER % 2) ? x0 : x0 + 0.5);
         y0 = VIGRA_CSTD::floor((ORDER % 2) ? y0 : y0 + 0.5);
         x1 = VIGRA_CSTD::floor((ORDER % 2) ? x1 : x1 + 0.5);
         y1 = VIGRA_CSTD::floor((ORDER % 2) ? y1 : y1 + 0.5);
         return x0 == x1 && y0 == y1;
    }

  protected:
    unsigned int w_, h_;
    INTERNAL_INDEXER internalIndexer_;
};

template <class VALUETYPE, class INTERNAL_TRAVERSER = typename BasicImage<VALUETYPE>::const_traverser>
class SplineImageView0
: public SplineImageView0Base<VALUETYPE, INTERNAL_TRAVERSER>
{
    typedef SplineImageView0Base<VALUETYPE, INTERNAL_TRAVERSER> Base;
  public:
    typedef typename Base::value_type value_type;
    typedef typename Base::size_type size_type;
    typedef typename Base::difference_type difference_type;
    enum StaticOrder { order = 0 };
    typedef BasicImage<VALUETYPE> InternalImage;
  
  protected:
    typedef typename IteratorTraits<INTERNAL_TRAVERSER>::mutable_iterator InternalTraverser;
    typedef typename IteratorTraits<InternalTraverser>::DefaultAccessor InternalAccessor;
    typedef typename IteratorTraits<INTERNAL_TRAVERSER>::const_iterator InternalConstTraverser;
    typedef typename IteratorTraits<InternalConstTraverser>::DefaultAccessor InternalConstAccessor;

  public:

        /* when traverser and accessor types passed to the constructor are the same as the corresponding
           internal types, we need not copy the image (speed up)
        */
    SplineImageView0(InternalTraverser is, InternalTraverser iend, InternalAccessor sa)
    : Base(iend.x - is.x, iend.y - is.y, is)
    {}

    SplineImageView0(triple<InternalTraverser, InternalTraverser, InternalAccessor> s)
    : Base(s.second.x - s.first.x, s.second.y - s.first.y, s.first)
    {}

    SplineImageView0(InternalConstTraverser is, InternalConstTraverser iend, InternalConstAccessor sa)
    : Base(iend.x - is.x, iend.y - is.y, is)
    {}

    SplineImageView0(triple<InternalConstTraverser, InternalConstTraverser, InternalConstAccessor> s)
    : Base(s.second.x - s.first.x, s.second.y - s.first.y, s.first)
    {}

    template<class T, class SU>
    SplineImageView0(MultiArrayView<2, T, SU> const & i)
    : Base(i.shape(0), i.shape(1)),
      image_(i.shape(0), i.shape(1))
    {
        for(unsigned int y=0; y<this->height(); ++y)
            for(unsigned int x=0; x<this->width(); ++x)
                image_(x,y) = detail::RequiresExplicitCast<VALUETYPE>::cast(i(x,y));
        this->internalIndexer_ = image_->upperLeft();
    }

    template <class SrcIterator, class SrcAccessor>
    SplineImageView0(SrcIterator is, SrcIterator iend, SrcAccessor sa)
    : Base(iend.x - is.x, iend.y - is.y),
      image_(iend - is)
    {
        copyImage(srcIterRange(is, iend, sa), destImage(image_));
        this->internalIndexer_ = image_->upperLeft();
    }

    template <class SrcIterator, class SrcAccessor>
    SplineImageView0(triple<SrcIterator, SrcIterator, SrcAccessor> s)
    : Base(s.second.x - s.first.x, s.second.y - s.first.y),
      image_(s.second - s.first)
    {
        copyImage(s, destImage(image_));
        this->internalIndexer_ = image_->upperLeft();
    }
    
    InternalImage const & image() const
        { return image_; }

  private:
    InternalImage image_;
};

template <class VALUETYPE, class StridedOrUnstrided>
class SplineImageView0<VALUETYPE, MultiArrayView<2, VALUETYPE, StridedOrUnstrided> >
: public SplineImageView0Base<VALUETYPE, MultiArrayView<2, VALUETYPE, StridedOrUnstrided> >
{
    typedef SplineImageView0Base<VALUETYPE, MultiArrayView<2, VALUETYPE, StridedOrUnstrided> > Base;
  public:
    typedef typename Base::value_type value_type;
    typedef typename Base::size_type size_type;
    typedef typename Base::difference_type difference_type;
    enum StaticOrder { order = 0 };
    typedef BasicImage<VALUETYPE> InternalImage;

  protected:
    typedef MultiArrayView<2, VALUETYPE, StridedOrUnstrided> InternalIndexer;
  
  public:

        /* when traverser and accessor types passed to the constructor are the same as the corresponding
           internal types, we need not copy the image (speed up)
        */
    SplineImageView0(InternalIndexer const & i)
    : Base(i.shape(0), i.shape(1), i)
    {}

    template<class T, class SU>
    SplineImageView0(MultiArrayView<2, T, SU> const & i)
    : Base(i.shape(0), i.shape(1)),
      image_(i.shape(0), i.shape(1))
    {
        for(unsigned int y=0; y<this->height(); ++y)
            for(unsigned int x=0; x<this->width(); ++x)
                image_(x,y) = detail::RequiresExplicitCast<VALUETYPE>::cast(i(x,y));
        this->internalIndexer_ = InternalIndexer(typename InternalIndexer::difference_type(this->width(), this->height()),
                                                 image_.data());
    }

    template <class SrcIterator, class SrcAccessor>
    SplineImageView0(SrcIterator is, SrcIterator iend, SrcAccessor sa)
    : Base(iend.x - is.x, iend.y - is.y),
      image_(iend-is)
    {
        copyImage(srcIterRange(is, iend, sa), destImage(image_));
        this->internalIndexer_ = InternalIndexer(typename InternalIndexer::difference_type(this->width(), this->height()),
                                                 image_.data());
    }

    template <class SrcIterator, class SrcAccessor>
    SplineImageView0(triple<SrcIterator, SrcIterator, SrcAccessor> s)
    : Base(s.second.x - s.first.x, s.second.y - s.first.y),
      image_(s.second - s.first)
    {
        copyImage(s, destImage(image_));
        this->internalIndexer_ = InternalIndexer(typename InternalIndexer::difference_type(this->width(), this->height()),
                                                 image_.data());
    }
    
    InternalImage const & image() const
        { return image_; }
    
  private:
    InternalImage image_;
};

template <class VALUETYPE>
class SplineImageView<0, VALUETYPE>
: public SplineImageView0<VALUETYPE>
{
    typedef SplineImageView0<VALUETYPE> Base;
  public:
    typedef typename Base::value_type value_type;
    typedef typename Base::size_type size_type;
    typedef typename Base::difference_type difference_type;
    enum StaticOrder { order = 0 };
    typedef typename Base::InternalImage InternalImage;
  
  protected:
    typedef typename Base::InternalTraverser InternalTraverser;
    typedef typename Base::InternalAccessor InternalAccessor;
    typedef typename Base::InternalConstTraverser InternalConstTraverser;
    typedef typename Base::InternalConstAccessor InternalConstAccessor;

public:

        /* when traverser and accessor types passed to the constructor are the same as the corresponding
           internal types, we need not copy the image (speed up)
        */
    SplineImageView(InternalTraverser is, InternalTraverser iend, InternalAccessor sa, bool /* unused */ = false)
    : Base(is, iend, sa)
    {}

    SplineImageView(triple<InternalTraverser, InternalTraverser, InternalAccessor> s, bool /* unused */ = false)
    : Base(s)
    {}

    SplineImageView(InternalConstTraverser is, InternalConstTraverser iend, InternalConstAccessor sa, bool /* unused */ = false)
    : Base(is, iend, sa)
    {}

    SplineImageView(triple<InternalConstTraverser, InternalConstTraverser, InternalConstAccessor> s, bool /* unused */ = false)
    : Base(s)
    {}

    template <class SrcIterator, class SrcAccessor>
    SplineImageView(SrcIterator is, SrcIterator iend, SrcAccessor sa, bool /* unused */ = false)
    : Base(is, iend, sa)
    {
        copyImage(srcIterRange(is, iend, sa), destImage(image_));
    }

    template <class SrcIterator, class SrcAccessor>
    SplineImageView(triple<SrcIterator, SrcIterator, SrcAccessor> s, bool /* unused */ = false)
    : Base(s)
    {
        copyImage(s, destImage(image_));
    }
};


} // namespace vigra


#endif /* VIGRA_SPLINEIMAGEVIEW_HXX */
