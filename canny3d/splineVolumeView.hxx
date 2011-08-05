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

#ifndef VIGRA_SPLINEVOLUMEVIEW_HXX
#define VIGRA_SPLINEVOLUMEVIEW_HXX

#include "vigra/mathutil.hxx"
#include "vigra/recursiveconvolution.hxx"
#include "vigra/splines.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/tinyvector.hxx"
#include "vigra/fixedpoint.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_pointoperators.hxx"
#include "vigra/navigator.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*                    SplineImageView                   */
/*                                                      */
/********************************************************/
  
  template <int ORDER, class VALUETYPE>
  class SplineVolumeView
  {
    typedef typename NumericTraits<VALUETYPE>::RealPromote InternalValue;
    
  public:
    
    typedef VALUETYPE value_type;    
    typedef TinyVector<double, 3> size_type;
    typedef TinyVector<double, 3> difference_type;

    typedef MultiArray<3, value_type> volume_type;    
    typedef MultiArray<3, double> InternalVolume;

    typedef typename AccessorTraits<value_type>::default_accessor default_accessor;

    enum StaticOrder { order = ORDER };
  
  private:
    //fix
    typedef typename InternalVolume::traverser InternalTraverser;
    //     typedef typename InternalTraverser::row_iterator InternalRowIterator;
    //     typedef typename InternalTraverser::column_iterator InternalColumnIterator;
    typedef BSpline<ORDER, double> Spline;
  
    enum { ksize_ = ORDER + 1, kcenter_ = ORDER / 2 };
    
  public:   
    
    //IMPL:d_,d1_,z0_,z1_,z_,vz_
    //IMPL:volume_
    template<typename SrcIterator, typename SrcShape, typename SrcAccessor>
    SplineVolumeView(SrcIterator iter, SrcShape shape, SrcAccessor acc, 
             bool skipPrefiltering = false)
    
      :  w_(shape[0]), h_(shape[1]), d_(shape[2]),
     w1_(w_-1), h1_(h_-1), d1_(d_-1),
     x0_(kcenter_), x1_(w_ - kcenter_ - 2), 
     y0_(kcenter_), y1_(h_ - kcenter_ - 2),
     z0_(kcenter_), z1_(d_ - kcenter_ - 2),
     volume_(shape),
     x_(-1.0), y_(-1.0), z_(-1.0),
     cx_(-1.0), cy_(-1.0), cz_(-1.0)
    {
      copyMultiArray(iter, shape, acc, 
             volume_.traverser_begin(), 
             AccessorTraits<value_type>::default_accessor());
    
      if(!skipPrefiltering)
    init();
    }

    SplineVolumeView(const volume_type & volume,  bool skipPrefiltering = false)
      :  w_(volume.shape()[0]), h_(volume.shape()[1]), d_(volume.shape()[2]),
     w1_(w_-1), h1_(h_-1), d1_(d_-1),
     x0_(kcenter_), x1_(w_ - kcenter_ - 2), 
     y0_(kcenter_), y1_(h_ - kcenter_ - 2),
     z0_(kcenter_), z1_(d_ - kcenter_ - 2),
     volume_(volume.shape()),
     x_(-1.0), y_(-1.0), z_(-1.0),
     cx_(-1.0), cy_(-1.0), cz_(-1.0)
    {
      copyMultiArray(volume.traverser_begin(), volume.shape(), default_accessor(), 
             volume_.traverser_begin(), 
             default_accessor());
      
      if(!skipPrefiltering)
    init();
    }
    
    value_type operator()(double x, double y, double z, unsigned int dx, unsigned int dy, unsigned int dz) const;

    value_type operator()(double x, double y, double z) const;
  
    value_type operator()(difference_type const & d) const
    { return operator()(d[0], d[1], d[2]); }

    value_type operator[](difference_type const & d) const
    { return operator()(d[0], d[1], d[2]); }

    difference_type d(difference_type const & d) const
    {  return difference_type(dx(d),dy(d),dz(d)); }

    value_type d(double x, double y, double z) const
    {  return difference_type(dx(d),dy(d),dz(d)); }

    value_type dx(double x, double y, double z) const
    { return operator()(x, y, z, 1, 0, 0); }
    
    value_type dy(double x, double y, double z) const
    { return operator()(x, y, z, 0, 1, 0); }
    
    value_type dz(double x, double y, double z) const
    { return operator()(x, y, z, 0, 0, 1); }

    value_type dxx(double x, double y, double z) const
    { return operator()(x, y, z, 2, 0, 0); }
    
    value_type dxy(double x, double y, double z) const
    { return operator()(x, y, z, 1, 1, 0); }
    
    value_type dxz(double x, double y, double z) const
    { return operator()(x, y, z, 1, 0, 1); }
    
    value_type dyy(double x, double y, double z) const
    { return operator()(x, y, z, 0, 2, 0); }
    
    value_type dyz(double x, double y, double z) const
    { return operator()(x, y, z, 0, 1, 1); }
    
    value_type dzz(double x, double y, double z) const
    { return operator()(x, y, z, 0, 0, 2); }
    
    value_type g2(double x, double y, double z) const
    { return sq(dx(x,y,z)) + sq(dy(x,y,z)) + sq(dz(x,y,z)); }
    
    value_type g2(difference_type const & d) const
    { return g2(d[0], d[1], d[2]); }

    value_type dx(difference_type const & d) const
    { return dx(d[0], d[1], d[2]); }

    value_type dy(difference_type const & d) const
    { return dy(d[0], d[1], d[2]); }

    value_type dz(difference_type const & d) const
    { return dz(d[0], d[1], d[2]); }

    value_type dxx(difference_type const & d) const
    { return operator()(d[0], d[1], d[2], 2, 0, 0); }
    
    value_type dyy(difference_type const & d) const
    { return operator()(d[0], d[1], d[2], 0, 2, 0); }

    value_type dzz(difference_type const & d) const
    { return operator()(d[0], d[1], d[2], 0, 0, 2); }
    
    value_type dxy(difference_type const & d) const
    { return operator()(d[0], d[1], d[2], 1, 1, 0); }
    
    value_type dxz(difference_type const & d) const
    { return operator()(d[0], d[1], d[2], 1, 0, 1); }
    
    value_type dyz(difference_type const & d) const
    { return operator()(d[0], d[1], d[2], 0, 1, 1); }
    
    unsigned int width() const
    { return w_; }
  
    unsigned int height() const
    { return h_; }

    unsigned int depth() const
    { return d_; }

    size_type size() const
    { return size_type(w_, h_, d_); }

    size_type shape() const
    { return size_type(w_, h_, d_); }
  
    InternalVolume const & volume() const
    {
      return volume_;
    }
    
    bool isInsideX(double x) const
    {
      return x >= 0.0 && x <= width()-1.0;
    }
  
    bool isInsideY(double y) const
    {
      return y >= 0.0 && y <= height()-1.0;
    }
  
    bool isInsideZ(double z) const
    {
      return z >= 0.0 && z <= depth()-1.0;
    }

    bool isInside(double x, double y, double z) const
    {
      return isInsideX(x)  && isInsideY(y) && isInsideZ(z);
    }

    bool isInside(difference_type const & d) const
    {
      return isInside(d[0],d[1],d[2]);
    }
    
    bool isValid(double x, double y, double z) const
    {
      return x < w1_ + x1_ && x > -x1_ && y < h1_ + y1_ && y > -y1_ && z < d1_ + z1_ && z > -z1_;
    }
    
    bool sameFacet(double x0, double y0, double z0,
           double x1, double y1, double z1) const
    {
      x0 = VIGRA_CSTD::floor((ORDER % 2) ? x0 : x0 + 0.5);
      y0 = VIGRA_CSTD::floor((ORDER % 2) ? y0 : y0 + 0.5);
      z0 = VIGRA_CSTD::floor((ORDER % 2) ? z0 : z0 + 0.5);
      x1 = VIGRA_CSTD::floor((ORDER % 2) ? x1 : x1 + 0.5);
      y1 = VIGRA_CSTD::floor((ORDER % 2) ? y1 : y1 + 0.5);
      z1 = VIGRA_CSTD::floor((ORDER % 2) ? z1 : z1 + 0.5);
    
      return x0 == x1 && y0 == y1 && z0 == z1;
    }
        
  protected:
  
    void init();
    void calculateIndices(double x, double y, double z) const;
    void coefficients(double t, double * const & c) const;
    void derivCoefficients(double t, unsigned int d, double * const & c) const;
    value_type convolve() const;
  
    unsigned int w_, h_, d_;
    int w1_, h1_, d1_;
    double x0_, x1_, y0_, y1_, z0_, z1_;
    InternalVolume volume_;
    Spline k_;
    mutable double x_, y_, z_, cx_, cy_, cz_, kx_[ksize_], ky_[ksize_], kz_[ksize_];
    mutable int ix_[ksize_], iy_[ksize_], iz_[ksize_];
  };
  
  template <int ORDER, class VALUETYPE>
  void SplineVolumeView<ORDER, VALUETYPE>::init()
  {


    typedef typename AccessorTraits<value_type>::default_accessor Accessor;
    typedef vigra::MultiArrayNavigator<InternalVolume::traverser, 3> Navigator;

    ArrayVector<double> const & b = k_.prefilterCoefficients();
    Accessor as;  
    
    typename InternalVolume::difference_type shape = volume_.shape();
    
    for(unsigned int i=0; i<b.size(); ++i)
      {       
    for(int d=0; d<3; ++d)
      {
        int s = shape[d];
        Navigator nav(volume_.traverser_begin(), shape, d);
        for(; nav.hasMore(); ++nav)
          {
        Navigator::iterator begin = nav.begin();
        
        recursiveFilterLine(begin, begin+s, as,
                    begin, as,
                    b[i], BORDER_TREATMENT_REFLECT);	      
          }
      }
      }
}
  
namespace detail
{

  // performs c[i] = c0 + i;
template <int i>
struct SplineVolumeViewUnrollLoop1
{
    template <class Array>
    static void exec(int c0, Array c)
    {
        SplineVolumeViewUnrollLoop1<i-1>::exec(c0, c);
        c[i] = c0 + i;
    }
};

template <>
struct SplineVolumeViewUnrollLoop1<0>
{
    template <class Array>
    static void exec(int c0, Array c)
    {
        c[0] = c0;
    }
};

  // performs sum_i(k[i] * r[x[i]])
template <int i, class ValueType>
struct SplineVolumeViewUnrollLoop2
{
    template <class Array1, class RowIterator, class Array2>
    static ValueType
    exec(Array1 k, RowIterator r, Array2 x)
    {
        return k[i] * r[x[i]] + SplineVolumeViewUnrollLoop2<i-1, ValueType>::exec(k, r, x);
    }
};

template <class ValueType>
struct SplineVolumeViewUnrollLoop2<0, ValueType>
{
    template <class Array1, class RowIterator, class Array2>
    static ValueType
    exec(Array1 k, RowIterator r, Array2 x)
    {
        return k[0] * r[x[0]];
    }
};

} // namespace detail

  // calculate u_,v_, ix_, iy_
template <int ORDER, class VALUETYPE>
void 
SplineVolumeView<ORDER, VALUETYPE>::calculateIndices(double x, double y, double z) const
{
    if(x == x_ && y == y_ && z == z_)
        return;   // still in cache
    
    if(x > x0_ && x < x1_ && y > y0_ && y < y1_ && z >z0_ && z < z1_)
    {
        detail::SplineVolumeViewUnrollLoop1<ORDER>::exec(
                                (ORDER % 2) ? int(x - kcenter_) : int(x + 0.5 - kcenter_),
                ix_);
        detail::SplineVolumeViewUnrollLoop1<ORDER>::exec(
                                (ORDER % 2) ? int(y - kcenter_) : int(y + 0.5 - kcenter_), 
                iy_);
        detail::SplineVolumeViewUnrollLoop1<ORDER>::exec(
                                (ORDER % 2) ? int(z - kcenter_) : int(z + 0.5 - kcenter_), 
                iz_);

        cx_ = x - ix_[kcenter_];
        cy_ = y - iy_[kcenter_];
        cz_ = z - iz_[kcenter_];
    }
    else
    {
      vigra_precondition(isValid(x,y,z),
                    "SplineVolumeView::calculateIndices(): coordinates out of range.");
        
        int xCenter = (ORDER % 2) ?
                      (int)VIGRA_CSTD::floor(x) :
                      (int)VIGRA_CSTD::floor(x + 0.5);
        int yCenter = (ORDER % 2) ?
                      (int)VIGRA_CSTD::floor(y) :
                      (int)VIGRA_CSTD::floor(y + 0.5);
        int zCenter = (ORDER % 2) ?
                      (int)VIGRA_CSTD::floor(z) :
                      (int)VIGRA_CSTD::floor(z + 0.5);

        
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
        if(z >= z1_)
        {
            for(int i = 0; i < ksize_; ++i)
                iz_[i] = d1_ - vigra::abs(d1_ - zCenter - (i - kcenter_));
        }
        else
        {
            for(int i = 0; i < ksize_; ++i)
                iz_[i] = vigra::abs(zCenter - (kcenter_ - i));
        }

        cx_ = x - xCenter;
        cy_ = y - yCenter;
    cz_ = z - zCenter;
    }
    x_ = x;
    y_ = y;
    z_ = z;
}

template <int ORDER, class VALUETYPE>
void SplineVolumeView<ORDER, VALUETYPE>::coefficients(double t, double * const & c) const
{
    t += kcenter_;
    for(int i = 0; i<ksize_; ++i)
        c[i] = k_(t-i);
}

template <int ORDER, class VALUETYPE>
void SplineVolumeView<ORDER, VALUETYPE>::derivCoefficients(double t, 
                                               unsigned int d, double * const & c) const
{
    t += kcenter_;
    for(int i = 0; i<ksize_; ++i)
        c[i] = k_(t-i, d);
}

template <int ORDER, class VALUETYPE>
VALUETYPE SplineVolumeView<ORDER, VALUETYPE>::convolve() const
{
  InternalValue sum = vigra::NumericTraits<InternalValue>::zero();
    
  for(int k=0; k<ksize_; ++k) {
    InternalValue sumy = vigra::NumericTraits<InternalValue>::zero();
    for(int j=0; j<ksize_; ++j) {
      InternalValue sumx = vigra::NumericTraits<InternalValue>::zero();
      for(int i=0; i<ksize_; ++i) {
    sumx += kx_[i] * volume_(ix_[i],iy_[j],iz_[k]);
      }
      sumy += ky_[j] * sumx;
    }
    sum += kz_[k] * sumy;
  }
      
  return NumericTraits<VALUETYPE>::fromRealPromote(sum);
}

// template <int ORDER, class VALUETYPE>
// template <class Array>
// void 
// SplineVolumeView<ORDER, VALUETYPE>::coefficientArray(double x, double y, double z,
// 						     Array & res) const
// {
//     typename Spline::WeightMatrix & weights = Spline::weights();
//     InternalValue tmp[ksize_][ksize_][ksize_]; 
    
//     calculateIndices(x, y, z);
//     for(int l=0; l<ksize_; ++l) 
//       {
// 	for(int j=0; j<ksize_; ++j)
// 	  {
// 	    for(int i=0; i<ksize_; ++i)
// 	      {
// 		tmp[i][j][l] = 0.0;
// 		for(int k=0; k<ksize_; ++k)
// 		  {
// 		    tmp[i][j][l] += weights[i][k]*image_(ix_[k], iy_[j]);
// 		  }
// 	      }       
// 	  }
//       }
//     res.resize(ksize_, ksize_);
//     for(int j=0; j<ksize_; ++j)
//     {
//         for(int i=0; i<ksize_; ++i)
//         {
//             res(i,j) = 0.0;
//             for(int k=0; k<ksize_; ++k)
//             {
//                 res(i,j) += weights[j][k]*tmp[i][k];
//             }
//         }       
//     }
// }
  
template <int ORDER, class VALUETYPE>
VALUETYPE SplineVolumeView<ORDER, VALUETYPE>::operator()(double x, double y, double z) const
{
  calculateIndices(x, y, z);
  coefficients(cx_, kx_);
  coefficients(cy_, ky_);
  coefficients(cz_, kz_);
  return convolve();
  //return volume_(int(x),int(y),int(z));
}
  
template <int ORDER, class VALUETYPE>
VALUETYPE SplineVolumeView<ORDER, VALUETYPE>::operator()
  (double x, double y, double z,
   unsigned int dx, unsigned int dy, unsigned int dz) const
{
  calculateIndices(x, y, z);
  derivCoefficients(cx_, dx, kx_);
  derivCoefficients(cy_, dy, ky_);
  derivCoefficients(cz_, dz, kz_);
  return convolve();
}

} // namespace vigra


#endif /* VIGRA_SPLINEIMAGEVIEW_HXX */
