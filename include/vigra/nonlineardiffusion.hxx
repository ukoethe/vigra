/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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

#ifndef VIGRA_NONLINEARDIFFUSION_HXX
#define VIGRA_NONLINEARDIFFUSION_HXX

#include <vector>
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/imageiteratoradapter.hxx"

namespace vigra {

template <class SrcIterator, class SrcAccessor,
          class CoeffIterator, class DestIterator>
void internalNonlinearDiffusionDiagonalSolver(
    SrcIterator sbegin, SrcIterator send, SrcAccessor sa,
    CoeffIterator diag, CoeffIterator upper, CoeffIterator lower,
    DestIterator dbegin)
{
    int w = send - sbegin - 1;
    
    int i;
    
    for(i=0; i<w; ++i)
    {
        lower[i] = lower[i] / diag[i];
        
        diag[i+1] = diag[i+1] - lower[i] * upper[i];
    }
    
    dbegin[0] = sa(sbegin);
    
    for(i=1; i<=w; ++i)
    {
        dbegin[i] = sa(sbegin, i) - lower[i-1] * dbegin[i-1];
    }
    
    dbegin[w] = dbegin[w] / diag[w];
    
    for(i=w-1; i>=0; --i)
    {
        dbegin[i] = (dbegin[i] - upper[i] * dbegin[i+1]) / diag[i];
    }
}


template <class SrcIterator, class SrcAccessor, 
          class WeightIterator, class WeightAccessor, 
          class DestIterator, class DestAccessor>
void internalNonlinearDiffusionAOSStep(
                   SrcIterator sul, SrcIterator slr, SrcAccessor as,
                   WeightIterator wul, WeightAccessor aw,
                   DestIterator dul, DestAccessor ad, double timestep)
{
    // use traits to determine SumType as to prevent possible overflow
    typedef typename 
        NumericTraits<typename DestAccessor::value_type>::RealPromote 
        DestType;
    
    typedef typename 
        NumericTraits<typename WeightAccessor::value_type>::RealPromote 
        WeightType;
        
    // calculate width and height of the image
    int w = slr.x - sul.x;  
    int h = slr.y - sul.y;  
    int d = (w < h) ? h : w;

    std::vector<WeightType> lower(d),
                            diag(d),
                            upper(d);

    std::vector<DestType> res(d);

    int x,y;
    
    WeightType one = NumericTraits<WeightType>::one();
    
     // create y iterators
    SrcIterator ys = sul;
    WeightIterator yw = wul;
    DestIterator yd = dul;
    
    // x-direction
    for(y=0; y<h; ++y, ++ys.y, ++yd.y, ++yw.y)
    {
        typename SrcIterator::row_iterator xs = ys.rowIterator();
        typename WeightIterator::row_iterator xw = yw.rowIterator();
        typename DestIterator::row_iterator xd = yd.rowIterator();

        // fill 3-diag matrix
        
        diag[0] = one + timestep * (aw(xw) + aw(xw, 1));
        for(x=1; x<w-1; ++x)
        {
            diag[x] = one + timestep * (2.0 * aw(xw, x) + aw(xw, x+1) + aw(xw, x-1));
        }
        diag[w-1] = one + timestep * (aw(xw, w-1) + aw(xw, w-2));

        for(x=0; x<w-1; ++x)
        {
            lower[x] = -timestep * (aw(xw, x) + aw(xw, x+1));
            upper[x] = lower[x];
        }
        
        internalNonlinearDiffusionDiagonalSolver(xs, xs+w, as, 
                            diag.begin(), upper.begin(), lower.begin(), res.begin());
                            
        for(x=0; x<w; ++x, ++xd)
        {
            ad.set(res[x], xd);
        }
    }
        
    // y-direction
    ys = sul;
    yw = wul;
    yd = dul;
    
    for(x=0; x<w; ++x, ++ys.x, ++yd.x, ++yw.x)
    {
        typename SrcIterator::column_iterator xs = ys.columnIterator();
        typename WeightIterator::column_iterator xw = yw.columnIterator();
        typename DestIterator::column_iterator xd = yd.columnIterator();

        // fill 3-diag matrix
        
        diag[0] = one + timestep * (aw(xw) + aw(xw, 1));
        for(y=1; y<h-1; ++y)
        {
            diag[y] = one + timestep * (2.0 * aw(xw, y) + aw(xw, y+1) + aw(xw, y-1));
        }
        diag[h-1] = one + timestep * (aw(xw, h-1) + aw(xw, h-2));

        for(y=0; y<h-1; ++y)
        {
            lower[y] = -timestep * (aw(xw, y) + aw(xw, y+1));
            upper[y] = lower[y];
        }
        
        internalNonlinearDiffusionDiagonalSolver(xs, xs+h, as, 
                            diag.begin(), upper.begin(), lower.begin(), res.begin());
                            
        for(y=0; y<h; ++y, ++xd)
        {
            ad.set(0.5 * (ad(xd) + res[y]), xd);
        }
    }
}

/** \addtogroup NonLinearDiffusion Non-linear Diffusion
    
    Perform edge-preserving smoothing.
*/
//@{

/********************************************************/
/*                                                      */
/*                  nonlinearDiffusion                  */
/*                                                      */
/********************************************************/

/** \brief Perform edge-preserving smoothing at the given scale.

    The algorithm solves the non-linear diffusion equation 
    
    \f[
        \frac{\partial}{\partial t} u =
        \frac{\partial}{\partial x}
          \left( g(|\nabla u|) 
                 \frac{\partial}{\partial x} u
          \right)
    \f]
    
    where <em> t</em> is the time, <b> x</b> is the location vector, 
    <em> u(</em><b> x</b><em> , t)</em> is the smoothed image at time <em> t</em>, and 
    <em> g(.)</em> is the location dependent diffusivity. At time zero, the image
    <em> u(</em><b> x</b><em> , 0)</em> is simply the original image. The time is
    propotional to the square of the scale parameter: <em> t = s^{2</em> / 2}.
    The diffusion
    equation is solved iteratively according 
    to the Additive Operator Splitting Scheme (AOS) from
    
    {\small
    J. Weickert: {\em "Recursive Separable Schemes for Nonlinear Diffusion 
    Filters"}, 
    in: B. ter Haar Romeny, L. Florack, J. Koenderingk, M. Viergever (eds.):
        1st Intl. Conf. on Scale-Space Theory in Computer Vision 1997,
        Springer LNCS 1252}

    <TT>DiffusivityFunctor</TT> implements the gradient dependent local diffusivity. 
    It is passed
    as an argument to \ref gradientBasedTransform(). The return value must be
    between 0 and 1 and determines the weight a pixel gets when
    its neighbors are smoothed. Weickert recommends the use of the diffusivity
    implemented by class \ref DiffusivityFunctor. It's also possible to use 
    other functors, for example one that always returns 1, in which case 
    we obtain the solution to the linear diffusion equation, i.e. 
    Gaussian convolution.
    
    The source value type must be a 
    linear space with internal addition, scalar multiplication, and
    NumericTraits defined. The value_type of the DiffusivityFunctor must be the 
    scalar field over wich the source value type's linear space is defined.
    
    In addition to <TT>nonlinearDiffusion()</TT>, there is an algorithm 
    <TT>nonlinearDiffusionExplicit()</TT> which implements the Explicit Scheme
    described in the above article. Both algorithms have the same interface,
    but the explicit scheme gives slightly more accurate approximations of 
    the diffusion process at the cost of much slower processing.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor,
                  class DiffusivityFunctor>
        void nonlinearDiffusion(SrcIterator sul, SrcIterator slr, SrcAccessor as,
                           DestIterator dul, DestAccessor ad,
                           DiffusivityFunctor const & weight, double scale);
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor,
                  class DiffusivityFunctor>
        void nonlinearDiffusion(
            triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            DiffusivityFunctor const & weight, double scale);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> "<a href="nonlineardiffusion_8hxx-source.html">vigra/nonlineardiffusion.hxx</a>"
    
    
    \code
    FImage src(w,h), dest(w,h); 
    float edge_threshold, scale;   
    ...
    
    nonlinearDiffusion(srcImageRange(src), destImage(dest),
                       DiffusivityFunctor<float>(edge_threshold), scale);
    \endcode

    <b> Required Interface:</b>
    
    <ul>
    
    <li> <TT>SrcIterator</TT> and <TT>DestIterator</TT> are models of ImageIterator
    <li> <TT>SrcAccessor</TT> and <TT>DestAccessor</TT> are models of StandardAccessor
    <li> <TT>SrcAccessor::value_type</TT> is a linear space
    <li> <TT>DiffusivityFunctor</TT> conforms to the requirements of
          \ref gradientBasedTransform(). Its range is between 0 and 1.
    <li> <TT>DiffusivityFunctor::value_type</TT> is an algebraic field
    
    </ul>
    
    <b> Precondition:</b>
    
    <TT>scale > 0</TT>
*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class DiffusivityFunc>
void nonlinearDiffusion(SrcIterator sul, SrcIterator slr, SrcAccessor as,
                   DestIterator dul, DestAccessor ad,
                   DiffusivityFunc const & weight, double scale)
{
    vigra_precondition(scale > 0.0, "nonlinearDiffusion(): scale must be > 0");
    
    double total_time = scale*scale/2.0;
    static const double time_step = 5.0;
    int number_of_steps = (int)(total_time / time_step);
    double rest_time = total_time - time_step * number_of_steps;
    
    int w = slr.x - sul.x;  
    int h = slr.y - sul.y;  
    Diff2D range(w,h);

    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote 
        TmpType;
    typedef typename DiffusivityFunc::value_type WeightType;
    
    BasicImage<TmpType> smooth1(w,h);
    BasicImage<TmpType> smooth2(w,h);
    
    BasicImage<WeightType> weights(w,h);
    
    typename BasicImage<TmpType>::Iterator s1 = smooth1.upperLeft(),
                                  s2 = smooth2.upperLeft();
    typename BasicImage<TmpType>::Accessor a = smooth1.accessor();
    
    typename BasicImage<WeightType>::Iterator wi = weights.upperLeft();
    typename BasicImage<WeightType>::Accessor wa = weights.accessor();

    gradientBasedTransform(sul, slr, as, wi, wa, weight);

    internalNonlinearDiffusionAOSStep(sul, slr, as, wi, wa, s1, a, rest_time);

    for(int i = 0; i < number_of_steps; ++i)
    {
        gradientBasedTransform(s1, s1+range, a, wi, wa, weight);
                      
        internalNonlinearDiffusionAOSStep(s1, s1+range, a, wi, wa, s2, a, time_step);
    
		std::swap(s1, s2);
    }
    
    copyImage(s1, s1+range, a, dul, ad);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class DiffusivityFunc>
inline
void nonlinearDiffusion(
    triple<SrcIterator, SrcIterator, SrcAccessor> src,
    pair<DestIterator, DestAccessor> dest,
    DiffusivityFunc const & weight, double scale)
{
    nonlinearDiffusion(src.first, src.second, src.third,
                           dest.first, dest.second,
                           weight, scale);
}

template <class SrcIterator, class SrcAccessor, 
          class WeightIterator, class WeightAccessor, 
          class DestIterator, class DestAccessor>
void internalNonlinearDiffusionExplicitStep(
                   SrcIterator sul, SrcIterator slr, SrcAccessor as,
                   WeightIterator wul, WeightAccessor aw,
                   DestIterator dul, DestAccessor ad,
                   double time_step)
{
    // use traits to determine SumType as to prevent possible overflow
    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote 
        SumType;
    
    typedef typename 
        NumericTraits<typename WeightAccessor::value_type>::RealPromote 
        WeightType;
        
    // calculate width and height of the image
    int w = slr.x - sul.x;  
    int h = slr.y - sul.y;  

    int x,y;
    
    static const Diff2D left(-1, 0);
    static const Diff2D right(1, 0);
    static const Diff2D top(0, -1);
    static const Diff2D bottom(0, 1);
    
    WeightType gt, gb, gl, gr, g0;
    WeightType one = NumericTraits<WeightType>::one();
    SumType sum;
    
    time_step /= 2.0;
    
    // create y iterators
    SrcIterator ys = sul;
    WeightIterator yw = wul;
    DestIterator yd = dul;
        
    SrcIterator xs = ys;
    WeightIterator xw = yw;
    DestIterator xd = yd;
    
    gt = (aw(xw) + aw(xw, bottom)) * time_step;
    gb = (aw(xw) + aw(xw, bottom)) * time_step;
    gl = (aw(xw) + aw(xw, right)) * time_step;
    gr = (aw(xw) + aw(xw, right)) * time_step;
    g0 = one - gt - gb - gr - gl;

    sum = g0 * as(xs);
    sum += gt * as(xs, bottom);
    sum += gb * as(xs, bottom);
    sum += gl * as(xs, right);
    sum += gr * as(xs, right);

    ad.set(sum, xd);

    for(x=2, ++xs.x, ++xd.x, ++xw.x; x<w; ++x, ++xs.x, ++xd.x, ++xw.x)
    {
        gt = (aw(xw) + aw(xw, bottom)) * time_step;
        gb = (aw(xw) + aw(xw, bottom)) * time_step;
        gl = (aw(xw) + aw(xw, left)) * time_step;
        gr = (aw(xw) + aw(xw, right)) * time_step;
        g0 = one - gt - gb - gr - gl;

        sum = g0 * as(xs);
        sum += gt * as(xs, bottom);
        sum += gb * as(xs, bottom);
        sum += gl * as(xs, left);
        sum += gr * as(xs, right);

        ad.set(sum, xd);
    }

    gt = (aw(xw) + aw(xw, bottom)) * time_step;
    gb = (aw(xw) + aw(xw, bottom)) * time_step;
    gl = (aw(xw) + aw(xw, left)) * time_step;
    gr = (aw(xw) + aw(xw, left)) * time_step;
    g0 = one - gt - gb - gr - gl;

    sum = g0 * as(xs);
    sum += gt * as(xs, bottom);
    sum += gb * as(xs, bottom);
    sum += gl * as(xs, left);
    sum += gr * as(xs, left);

    ad.set(sum, xd);
    
    for(y=2, ++ys.y, ++yd.y, ++yw.y; y<h; ++y, ++ys.y, ++yd.y, ++yw.y)
    {
        xs = ys;
        xd = yd;
        xw = yw;
        
        gt = (aw(xw) + aw(xw, top)) * time_step;
        gb = (aw(xw) + aw(xw, bottom)) * time_step;
        gl = (aw(xw) + aw(xw, right)) * time_step;
        gr = (aw(xw) + aw(xw, right)) * time_step;
        g0 = one - gt - gb - gr - gl;

        sum = g0 * as(xs);
        sum += gt * as(xs, top);
        sum += gb * as(xs, bottom);
        sum += gl * as(xs, right);
        sum += gr * as(xs, right);

        ad.set(sum, xd);
        
        for(x=2, ++xs.x, ++xd.x, ++xw.x; x<w; ++x, ++xs.x, ++xd.x, ++xw.x)
        {
            gt = (aw(xw) + aw(xw, top)) * time_step;
            gb = (aw(xw) + aw(xw, bottom)) * time_step;
            gl = (aw(xw) + aw(xw, left)) * time_step;
            gr = (aw(xw) + aw(xw, right)) * time_step;
            g0 = one - gt - gb - gr - gl;
            
            sum = g0 * as(xs);
            sum += gt * as(xs, top);
            sum += gb * as(xs, bottom);
            sum += gl * as(xs, left);
            sum += gr * as(xs, right);
            
            ad.set(sum, xd);
        }
        
        gt = (aw(xw) + aw(xw, top)) * time_step;
        gb = (aw(xw) + aw(xw, bottom)) * time_step;
        gl = (aw(xw) + aw(xw, left)) * time_step;
        gr = (aw(xw) + aw(xw, left)) * time_step;
        g0 = one - gt - gb - gr - gl;

        sum = g0 * as(xs);
        sum += gt * as(xs, top);
        sum += gb * as(xs, bottom);
        sum += gl * as(xs, left);
        sum += gr * as(xs, left);

        ad.set(sum, xd);
    }
    
    xs = ys;
    xd = yd;
    xw = yw;

    gt = (aw(xw) + aw(xw, top)) * time_step;
    gb = (aw(xw) + aw(xw, top)) * time_step;
    gl = (aw(xw) + aw(xw, right)) * time_step;
    gr = (aw(xw) + aw(xw, right)) * time_step;
    g0 = one - gt - gb - gr - gl;

    sum = g0 * as(xs);
    sum += gt * as(xs, top);
    sum += gb * as(xs, top);
    sum += gl * as(xs, right);
    sum += gr * as(xs, right);

    ad.set(sum, xd);

    for(x=2, ++xs.x, ++xd.x, ++xw.x; x<w; ++x, ++xs.x, ++xd.x, ++xw.x)
    {
        gt = (aw(xw) + aw(xw, top)) * time_step;
        gb = (aw(xw) + aw(xw, top)) * time_step;
        gl = (aw(xw) + aw(xw, left)) * time_step;
        gr = (aw(xw) + aw(xw, right)) * time_step;
        g0 = one - gt - gb - gr - gl;

        sum = g0 * as(xs);
        sum += gt * as(xs, top);
        sum += gb * as(xs, top);
        sum += gl * as(xs, left);
        sum += gr * as(xs, right);

        ad.set(sum, xd);
    }

    gt = (aw(xw) + aw(xw, top)) * time_step;
    gb = (aw(xw) + aw(xw, top)) * time_step;
    gl = (aw(xw) + aw(xw, left)) * time_step;
    gr = (aw(xw) + aw(xw, left)) * time_step;
    g0 = one - gt - gb - gr - gl;

    sum = g0 * as(xs);
    sum += gt * as(xs, top);
    sum += gb * as(xs, top);
    sum += gl * as(xs, left);
    sum += gr * as(xs, left);

    ad.set(sum, xd);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class DiffusivityFunc>
void nonlinearDiffusionExplicit(SrcIterator sul, SrcIterator slr, SrcAccessor as,
                   DestIterator dul, DestAccessor ad,
                   DiffusivityFunc const & weight, double scale)
{
    vigra_precondition(scale > 0.0, "nonlinearDiffusionExplicit(): scale must be > 0");
    
    double total_time = scale*scale/2.0;
    static const double time_step = 0.25;
    int number_of_steps = total_time / time_step;
    double rest_time = total_time - time_step * number_of_steps;
    
    int w = slr.x - sul.x;  
    int h = slr.y - sul.y;  
    Diff2D range(w,h);

    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote 
        TmpType;
    typedef typename DiffusivityFunc::value_type WeightType;
    
    BasicImage<TmpType> smooth1(w,h);
    BasicImage<TmpType> smooth2(w,h);
    
    BasicImage<WeightType> weights(w,h);
    
    typename BasicImage<TmpType>::Iterator s1 = smooth1.upperLeft(),
                                  s2 = smooth2.upperLeft();
    typename BasicImage<TmpType>::Accessor a = smooth1.accessor();
    
    typename BasicImage<WeightType>::Iterator wi = weights.upperLeft();
    typename BasicImage<WeightType>::Accessor wa = weights.accessor();

    gradientBasedTransform(sul, slr, as, wi, wa, weight);

    internalNonlinearDiffusionExplicitStep(sul, slr, as, wi, wa, s1, a, rest_time);

    for(int i = 0; i < number_of_steps; ++i)
    {
        gradientBasedTransform(s1, s1+range, a, wi, wa, weight);
                      
        internalNonlinearDiffusionExplicitStep(s1, s1+range, a, wi, wa, s2, a, time_step);
    
        swap(s1, s2);
    }
    
    copyImage(s1, s1+range, a, dul, ad);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class DiffusivityFunc>
inline
void nonlinearDiffusionExplicit(
    triple<SrcIterator, SrcIterator, SrcAccessor> src,
    pair<DestIterator, DestAccessor> dest,
    DiffusivityFunc const & weight, double scale)
{
    nonlinearDiffusionExplicit(src.first, src.second, src.third,
                           dest.first, dest.second,
                           weight, scale);
}

/********************************************************/
/*                                                      */
/*                   DiffusivityFunctor                 */
/*                                                      */
/********************************************************/

/** \brief Diffusivity functor for non-linear diffusion. 

    This functor implements the diffusivity recommended by Weickert:
    
    \f[
        g(|\nabla u|) = 1 - 
           \exp{\left(\frac{-3.315}{(|\nabla u| / thresh)^4}\right)}
    \f]
    
    
    where <TT>thresh</TT> is a threshold that determines whether a specific gradient 
    magnitude is interpreted as a significant edge (above the threshold) 
    or as noise. It is meant to be used with \ref nonlinearDiffusion().
*/
template <class Value>
class DiffusivityFunctor
{
  public:
         /** the functors value type (must be a field with <TT>sqrt()</TT> defined)
         */
    typedef Value value_type;
    
         /** init functor with given edge threshold
         */
    DiffusivityFunctor(Value const & thresh)
    : weight_(thresh*thresh), 
      one_(NumericTraits<Value>::one()),
      zero_(NumericTraits<Value>::zero())
    {}
    
         /** calculate diffusivity from scalar arguments
         */
    Value operator()(Value const & gx, Value const & gy) const
    {
        Value mag = (gx*gx + gy*gy) / weight_;
                     
        return (mag == zero_) ? one_ : one_ - exp(-3.315 / mag / mag);
    }
    
         /** calculate diffusivity from RGB arguments
         */
    Value operator()(RGBValue<Value> const & gx, RGBValue<Value> const & gy) const
    {
        Value mag = (gx.red()*gx.red() + 
                     gx.green()*gx.green() +
                     gx.blue()*gx.blue() +
                     gy.red()*gy.red() + 
                     gy.green()*gy.green() +
                     gy.blue()*gy.blue()) / weight_;

        return (mag == zero_) ? one_ : one_ - exp(-3.315 / mag / mag);
    }
    
    Value weight_;
    Value one_;
    Value zero_;
};

//@}

} // namespace vigra

#endif /* VIGRA_NONLINEARDIFFUSION_HXX */
