/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_SPLINEIMAGEVIEW_HXX
#define VIGRA_SPLINEIMAGEVIEW_HXX

#include "vigra/mathutil.hxx"
#include "vigra/recursiveconvolution.hxx"
#include "vigra/splines.hxx"
#include "vigra/array_vector.hxx"

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
        */
    value_type operator()(double x, double y) const;
    
        /** Access derivative of order <tt>(dx, dy)</tt> at real-valued coordinate <tt>(x, y)</tt>.
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
    
        /** Check if x is in the valid range.
            Equivalent to <tt>0 &lt;= x &lt;= width()-1</tt>.
        */
    bool isInsideX(double x) const
    {
        return x >= 0.0 && x <= width()-1.0;
    }
        
        /** Check if y is in the valid range.
            Equivalent to <tt>0 &lt;= y &lt;= height()-1</tt>.
        */
    bool isInsideY(double y) const
    {
        return y >= 0.0 && y <= height()-1.0;
    }
        
        /** Check if x and y are in the valid range.
            Equivalent to <tt>0 &lt;= x &lt;= width()-1</tt> and <tt>0 &lt;= y &lt;= height()-1</tt>.
        */
    bool isInside(double x, double y) const
    {
        return isInsideX(x) && isInsideY(y);
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

template <int i>
struct SplineImageViewUnrollLoop2
{
    template <class Array1, class Image, class Array2>
    static typename Image::value_type
    exec(Array1 k, Image const & img, Array2 x, int y)
    {
        return k[i] * img(x[i], y) + SplineImageViewUnrollLoop2<i-1>::exec(k, img, x, y);
    }
};

template <>
struct SplineImageViewUnrollLoop2<0>
{
    template <class Array1, class Image, class Array2>
    static typename Image::value_type
    exec(Array1 k, Image const & img, Array2 x, int y)
    {
        return k[0] * img(x[0], y);
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
    }
    else
    {
        vigra_precondition(isInside(x, y),
             "SplineImageView<ORDER, VALUETYPE>::calculateIndices(): index out of bounds.");

        ix_[kcenter_] = (ORDER % 2) ?
                             (int)x :
                             (int)(x + 0.5);
        iy_[kcenter_] = (ORDER % 2) ?
                             (int)y :
                             (int)(y + 0.5);
        
        for(int i=0; i<kcenter_; ++i)
        {
            ix_[i] = vigra::abs(ix_[kcenter_] - (kcenter_ - i));
            iy_[i] = vigra::abs(iy_[kcenter_] - (kcenter_ - i));
        }
        for(int i=kcenter_+1; i<ksize_; ++i)
        {
            ix_[i] = w1_ - vigra::abs(w1_ - ix_[kcenter_] - (i - kcenter_));
            iy_[i] = h1_ - vigra::abs(h1_ - iy_[kcenter_] - (i - kcenter_));
        }
    }
    x_ = x;
    y_ = y;
    u_ = x - ix_[kcenter_];
    v_ = y - iy_[kcenter_];
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
    sum = ky_[0]*detail::SplineImageViewUnrollLoop2<ORDER>::exec(kx_, image_, ix_, iy_[0]);
    
    for(int j=1; j<ksize_; ++j)
    {
        sum += ky_[j]*detail::SplineImageViewUnrollLoop2<ORDER>::exec(kx_, image_, ix_, iy_[j]);
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

} // namespace vigra


#endif /* VIGRA_SPLINEIMAGEVIEW_HXX */
