/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2014 by Benjamin Seppke                 */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_SHOCKFILTER_HXX
#define VIGRA_SHOCKFILTER_HXX

#include "basicimage.hxx"
#include "convolution.hxx"
#include "tensorutilities.hxx"

namespace vigra {
    

/********************************************************/
/*                                                      */
/*           Coherence enhancing shock filter           */
/*                                                      */
/********************************************************/

/**  
    This function calculates of the coherence enhancing shock filter proposed by 
    J. Weickert (2002): Coherence-Enhancing Show Filters. 
    It uses the structure tensor information of an image and an iterative discrete upwinding scheme
    instead of pure dilation and erosion to perform the necessary morphological operations
    on the image. 
*/
//@{

/** \brief This function calculates discrete upwinding scheme proposed by J. Weickert (2002) in Coherence-Enhancing Show Filters.

    An image is upwinded positively (dilation), if the given second src is positive.
    Otherwise it is upwinds negatively (eroded). The effect can be steered by an upwinding
    factor.
    
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2
                  class T3, class S3>
        void
        upwindImage(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> const & src2,
                    MultiArrayView<2, T3, S3> dest,
                    float upwind_factor_h);

    }
    \endcode

    \deprecatedAPI{upwindImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator,  class SrcAccessor, 
                  class Src2Iterator, class Src2Accessor,
                  class DestIterator, class DestAccessor>
        void upwindImage(SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,
                         Src2Iterator s2_ul, Src2Accessor s2_acc, 
                         DestIterator d_ul, DestAccessor d_acc,
                         float upwind_factor_h)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator,  class SrcAccessor, 
                  class Src2Iterator, class Src2Accessor,
                  class DestIterator, class DestAccessor>
        void
        upwindImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<Src2Iterator, Src2Accessor> src2,
                    pair<DestIterator, DestAccessor> dest,
                    float upwind_factor_h);
    }
    \endcode
    \deprecatedEnd
*/


doxygen_overloaded_function(template <...> void upwindImage)

template <class SrcIterator,  class SrcAccessor, 
          class Src2Iterator, class Src2Accessor,
          class DestIterator, class DestAccessor>
void upwindImage(SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,
                 Src2Iterator s2_ul, Src2Accessor s2_acc, 
                 DestIterator d_ul, DestAccessor d_acc,
                 float upwind_factor_h)
{
    using namespace std;
    
    typedef typename SrcIterator::difference_type  DiffType;
    
    DiffType shape = s_lr - s_ul;
    
    typedef typename SrcAccessor::value_type  SrcType;
    typedef typename DestAccessor::value_type ResultType;
    
    SrcType upper, lower, left, right, center;
    ResultType fx, fy;  
    
    
    for(int y=0; y<shape[1]; ++y)
    {
        for(int x=0; x<shape[0]; ++x)
        {
            upper  = s_acc(s_ul + Diff2D(x, max(0, y-1)));
            lower  = s_acc(s_ul + Diff2D(x, min(shape[1]-1, y+1)));
            left   = s_acc(s_ul + Diff2D(max(0, x-1), y));
            right  = s_acc(s_ul + Diff2D(min(shape[0]-1, x+1), y));
            center = s_acc(s_ul + Diff2D(x, y));
            
            if(s2_acc(s2_ul+Diff2D(x,y))<0)
            {
                fx = max(max(right - center, left  - center), 0.0f);
                fy = max(max(lower - center, upper - center), 0.0f);
                d_acc.set (center + upwind_factor_h*sqrt( fx*fx + fy*fy), d_ul+Diff2D(x,y));
            }
            else
            {
                fx = max(max(center - right, center - left), 0.0f);
                fy = max(max(center - lower, center - upper), 0.0f);
                d_acc.set (center - upwind_factor_h*sqrt( fx*fx + fy*fy), d_ul+Diff2D(x,y));
            }               
        }
    }
}

template <class SrcIterator,  class SrcAccessor, 
          class Src2Iterator, class Src2Accessor, 
          class DestIterator, class DestAccessor>
inline void upwindImage(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                        pair<Src2Iterator, Src2Accessor> s2, 
                        pair<DestIterator, DestAccessor> d,
                        float upwind_factor_h)
{
    upwindImage(s.first, s.second, s.third, s2.first, s2.second, d.first, d.second, upwind_factor_h);
}

template <class T1, class S1, 
          class T2, class S2,
          class T3, class S3>
inline void upwindImage(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, T2, S2> const & src2,
                        MultiArrayView<2, T3, S3> dest,
                        float upwind_factor_h)
{
    vigra_precondition(src.shape() == src2.shape() && src.shape() == dest.shape(),
                        "vigra::upwindImage(): shape mismatch between input and output.");
    upwindImage(srcImageRange(src),
                srcImage(src2),
                destImage(dest), 
                upwind_factor_h);
}


/** \brief This function calculates of the coherence enhancing shock filter proposed by J. Weickert (2002) in Coherence-Enhancing Show Filters.
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        shockFilter(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> dest,
                    float sigma, float rho, float upwind_factor_h, 
                    unsigned int iterations);

    }
    \endcode

    \deprecatedAPI{shockFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void shockFilter(SrcIterator supperleft,
                         SrcIterator slowerright, SrcAccessor sa,
                         DestIterator dupperleft, DestAccessor da,
                         float sigma, float rho, float upwind_factor_h, 
                         unsigned int iterations);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        shockFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    float sigma, float rho, float upwind_factor_h, 
                    unsigned int iterations);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/shockilter.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int w=1000, h=1000;
    MultiArray<2, float> src(w,h), dest(w,h);
    ...
    
 
    
    // apply a shock-filter:
    shockFilter(src, dest, 1.0, 5.0, 0.3, 5);
    \endcode

    <b> Preconditions:</b>

    The image must be larger than the window size of the filter.
*/

doxygen_overloaded_function(template <...> void upwindImage)

template <class SrcIterator,  class SrcAccessor, 
          class DestIterator, class DestAccessor>
void shockFilter(SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,
                 DestIterator d_ul, DestAccessor d_acc,
                 float sigma, float rho, float upwind_factor_h, 
                 unsigned int iterations)
{
    
    typedef typename SrcIterator::difference_type  DiffType;
    DiffType shape = s_lr - s_ul;
        
    unsigned int    w = shape[0],
                    h = shape[1];
    
    FVector3Image tensor(w,h);
    FVector3Image eigen(w,h);
    FImage hxx(w,h), hxy(w,h), hyy(w,h), temp(w,h) ,result(w,h);
    
    float c, s, v_xx, v_xy, v_yy;
    
    copyImage(srcIterRange(s_ul, s_lr, s_acc), destImage(result));
                     
    for(unsigned int i = 0; i<iterations; ++i)
    {   
        
        structureTensor(srcImageRange(result), destImage(tensor), sigma, rho);
        tensorEigenRepresentation(srcImageRange(tensor), destImage(eigen));
        hessianMatrixOfGaussian(srcImageRange(result),
                                destImage(hxx), destImage(hxy), destImage(hyy), sigma);
        
        for(int y=0; y<shape[1]; ++y)
        {
            for(int x=0; x<shape[0]; ++x)
            {
                c = cos(eigen(x,y)[2]);
                s = sin(eigen(x,y)[2]);
                v_xx = hxx(x,y);
                v_xy = hxy(x,y);
                v_yy = hyy(x,y);
                //store signum image in hxx (safe, because no other will ever access hxx(x,y)
                hxx(x,y) = c*c*v_xx + 2*c*s*v_xy + s*s*v_yy;
            }
        }
        upwindImage(srcImageRange(result),srcImage(hxx), destImage(temp), upwind_factor_h);
        result = temp;

    }
    copyImage(srcImageRange(result), destIter(d_ul, d_acc));
}

template <class SrcIterator,  class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void shockFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                                        pair<DestIterator, DestAccessor> d,
                                        float sigma, float rho, float upwind_factor_h, 
                                        unsigned int iterations)
{
    shockFilter(s.first, s.second, s.third, 
                d.first, d.second, 
                sigma, rho, upwind_factor_h, 
                iterations);
}       

template <class T1, class S1, 
          class T2, class S2>
inline void shockFilter(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, T2, S2> dest,
                        float sigma, float rho, float upwind_factor_h, 
                        unsigned int iterations)
{
    vigra_precondition(src.shape() == dest.shape(),
                        "vigra::shockFilter(): shape mismatch between input and output.");
    shockFilter(srcImageRange(src),
                destImage(dest), 
                sigma, rho, upwind_factor_h, 
                iterations);
}

} //end of namespace vigra

#endif //VIGRA_SHOCKFILTER_HXX
