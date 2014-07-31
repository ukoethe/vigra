/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2014 by Benjamin Seppke                 */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_CORRELATION_HXX
#define VIGRA_CORRELATION_HXX

#include "stdimage.hxx"
#include "inspectimage.hxx"
#include "functorexpression.hxx"
#include "fftw3.hxx"

// "slow" correlation algorithms are performed using the windowing filter env.
#include "applywindowfunction.hxx"


namespace vigra
{

namespace detail
{

struct FourierCorrelationFunctor
{
    FourierCorrelationFunctor(double normalization = 1.0)
    : normalization_factor(normalization)
    {
    }
    template<class T>
    FFTWComplex<T> operator()(FFTWComplex<T> const& u, FFTWComplex<T> const& v) const
    {
        return u*conj(v)/ normalization_factor;
    }
    
    double normalization_factor;
};

} //End of namespace detail
  


/********************************************************/
/*                                                      */
/*  Fast normalized cross correlation of mask to image  */
/*                                                      */
/********************************************************/
/**
    This function performes a fast normalized cross-correlation using the Fast Fourier Transform
    and sum-image look-up-tables. The implementation is based on the article by J.P.Lewis (1995):
    "Fast Normalized Cross-Correlation".
 */

//@{

/** \brief This function performes a fast normalized cross-correlation

    To compute the normalized cross correlation in a fast way, it is using the Fast Fourier Transform
    and sum-image look-up-tables as it is suggested by J.P.Lewis (1995):
    "Fast Normalized Cross-Correlation".

    The input pixel type <tt>T1</tt> must be a \ref LinearSpace "linear space" over 
    the window functions' value_type <tt>T</tt>, i.e. addition of source values, multiplication with functions' values,
    and NumericTraits must be defined. The mask's value_type must be an \ref AlgebraicField "algebraic field",
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must be defined. 
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class T3, class S3>
        void
        fastNormalizedCrossCorrelation(MultiArrayView<2, T1, S1> const & mask,
                                       MultiArrayView<2, T2, S2> const & src,
                                       MultiArrayView<2, T3, S3> dest);

    }
    \endcode

    \deprecatedAPI{fastNormalizedCrossCorrelation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class MaskIterator, class MaskAccessor,
                  class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void fastNormalizedCrossCorrelation(MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc,
                                            SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,  
                                            DestIterator d_ul, DestAccessor d_acc);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class MaskIterator, class MaskAccessor,
                  class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        fastNormalizedCrossCorrelation(triple<MaskIterator, MaskIterator, MaskAccessor> mask,
                                       triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                       pair<DestIterator, DestAccessor> dest);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/correlation.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int m_w=51, m_h=51;
    unsigned int w=1000, h=1000;
    MultiArray<2, float> mask(m_w,m_h), src(w,h), dest(w,h);
    ...
    
    //compute normalized cross correlation of mask and image -> dest
    fastNormalizedCrossCorrelation(mask, src, dest);
    \endcode
    
    <b> Preconditions:</b>

    The image must be larger than the size of the mask.
*/

doxygen_overloaded_function(template <...> void fastNormalizedCrossCorrelation)

template <class MaskIterator, class MaskAccessor,
          class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void fastNormalizedCrossCorrelation(MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc,
                                    SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,  
                                    DestIterator d_ul, DestAccessor d_acc)
{
    
    using namespace vigra;
    
    typedef typename MaskIterator::difference_type MaskDiffType;
    typedef typename SrcIterator::difference_type  SrcDiffType;
    typedef typename DestIterator::difference_type DestDiffType;
    
    MaskDiffType    mask_shape  = m_lr - m_ul;
    SrcDiffType     image_shape = s_lr - s_ul;
    
    const int   s_w = image_shape.x,    //Image dimensions
                s_h = image_shape.y,    //--
    
                m_w = mask_shape.x,     //Mask dimensions
                m_h = mask_shape.y,     //--
                m_total = m_w*m_h,
    
                padded_w = s_w+m_w,     //Padded dimensions
                padded_h = s_h+m_h,     //--
                padded_total = padded_w*padded_h;
    
    vigra_precondition( m_w % 2 == 1 , "vigra::fastNormalizedCrossCorrelation(): Mask width has to be of odd size!");
    vigra_precondition( m_h % 2 == 1 , "vigra::fastNormalizedCrossCorrelation(): Mask height has to be of odd size!");
    
    vigra_precondition( m_w <= s_w && m_h <= s_h , "vigra::fastNormalizedCrossCorrelation(): Mask is larger than image!");
    
    //find mask mean
    FindAverage<double> average;   
    inspectImage(srcIterRange(m_ul, m_lr, m_acc), average);
    
    //find mask sum and mask^2 sum
    double  mask_sum  = 0.0,
    mask_sum2 = 0.0;    
    
    MaskIterator ym = m_ul;
    MaskIterator xm = ym;
    
    for( ; ym.y != m_lr.y; ym.y++)
    {   
        for(xm = ym; xm.x != m_lr.x; xm.x++)
        {
            mask_sum  += m_acc(xm);
            mask_sum2 += m_acc(xm) * m_acc(xm);
        }
    }       
    
    //calculate the fix part of the denumerator
    double fix_denumerator = sqrt(m_total*mask_sum2 - mask_sum*mask_sum);
    
    if (fix_denumerator == 0)
    {
        //set correlation result d to zeros only...
        DestIterator yd = d_ul;
        DestIterator xd = yd;
        DestIterator d_lr = d_ul + image_shape;
        
        for( ; yd.y != d_lr.y; yd.y++)
        {   
            for(xd = yd; xd.x != d_lr.x; xd.x++)
            {
                d_acc.set(0.0, xd);
            }
        }   
    }
    else
    {
        //padding images to enable border cases
        FImage   padded_mask(padded_w, padded_h),
                 padded_image(padded_w, padded_h);
        
        
        //fill padded mask  
        transformImage(srcIterRange(m_ul, m_lr, m_acc), destImage(padded_mask), functor::Arg1() - functor::Param(average()));
        
        //fill padded image
        copyImage(srcIterRange(s_ul, s_lr, s_acc), destImage(padded_image));
        
        
        //calculate (unnormalized) numerator:
        
        //create fourier images
        FFTWComplexImage fourier_image(padded_w,padded_h),
        fourier_mask(padded_w,padded_h);
        //transform both
        fourierTransform(srcImageRange(padded_image), destImage(fourier_image));
        fourierTransform(srcImageRange(padded_mask), destImage(fourier_mask));
        
        detail::FourierCorrelationFunctor corr(padded_total);
        
        //do the correlation in frequency space
        combineTwoImages(srcImageRange(fourier_image), srcImage(fourier_mask), 
                         destImage(fourier_image), 
                         corr);
        
        //and go back to normal space
        fourierTransformInverse(srcImageRange(fourier_image), destImage(fourier_mask));
        
        
        //Create fast sum tables for the variable denumerator
        DImage sum_table(padded_w,padded_h),
               sum_table2(padded_w,padded_h);
        for(int u=0; u<padded_w-1; u++)
        {
            for(int v=0; v<padded_h-1; v++)
            {
                sum_table( u+1,v+1)  = padded_image(u,v)                   + sum_table(u,v+1)  + sum_table(u+1,v)   - sum_table(u,v);
                sum_table2(u+1,v+1)  = padded_image(u,v)*padded_image(u,v) + sum_table2(u,v+1) + sum_table2(u+1,v)  - sum_table2(u,v);
            }
        }
        
        //calculate the result, use the sum tables for the denumerators
        for(int v=0; v<=s_h-m_h; v++)
        {
            for(int u=0; u<=s_w-m_w; u++)
            {
                //calculate e(window) and e(window^2)
                double e_uv   = sum_table( u+m_w, v+m_h) - sum_table( u, v+m_h) - sum_table( u+m_w, v)  + sum_table( u,v),
                       e_uv_2 = sum_table2(u+m_w, v+m_h) - sum_table2(u, v+m_h) - sum_table2(u+m_w, v)  + sum_table2(u,v),
                       var_denumerator = sqrt(m_total*e_uv_2 - e_uv*e_uv);
                
                //calclate overall result
                if(var_denumerator == 0)
                {
                    d_acc.set(0.0, d_ul + DestDiffType(u+m_w/2,v+m_h/2));
                }
                else
                {           
                    d_acc.set( m_total*(fourier_mask(u,v)).re()/(var_denumerator*fix_denumerator),
                              d_ul + DestDiffType(u+m_w/2,v+m_h/2));
                }
            }
        }
    }
}

template <class MaskIterator, class MaskAccessor,
          class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void fastNormalizedCrossCorrelation(triple<MaskIterator, MaskIterator, MaskAccessor> mask,
                                           triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                           std::pair<DestIterator, DestAccessor> dest)
{
    fastNormalizedCrossCorrelation(mask.first, mask.second, mask.third, 
                                   src.first, src.second, src.third, 
                                   dest.first, dest.second);
}

template <class T1, class S1,
          class T2, class S2, 
          class T3, class S3>
inline void fastNormalizedCrossCorrelation(MultiArrayView<2, T1, S1> const & mask,
                                           MultiArrayView<2, T2, S2> const & src,
                                           MultiArrayView<2, T3, S3> dest)
{   
vigra_precondition(src.shape() == dest.shape(),
                        "vigra::fastNormalizedCrossCorrelation(): shape mismatch between input and output.");
    fastNormalizedCrossCorrelation(srcImageRange(mask), 
                                   srcImageRange(src), 
                                   destImage(dest));
}



/********************************************************/
/*                                                      */
/*        Fast cross correlation of mask to image       */
/*                                                      */
/********************************************************/
/**
    This function performes a fast cross-correlation using the Fast Fourier Transform and
    the dependency of the convolution and the correlation in Fourier space.
 */
 
/** \brief This function performes a fast cross-correlation

    This function performes a fast cross-correlation using the Fast Fourier Transform and
    the dependency of the convolution and the correlation in Fourier space.

    The input pixel type <tt>T1</tt> must be a \ref LinearSpace "linear space" over 
    the window functions' value_type <tt>T</tt>, i.e. addition of source values, multiplication with functions' values,
    and NumericTraits must be defined. The mask's value_type must be an \ref AlgebraicField "algebraic field",
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must be defined. 
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class T3, class S3>
        void
        fastCrossCorrelation(MultiArrayView<2, T1, S1> const & mask,
                             MultiArrayView<2, T2, S2> const & src,
                             MultiArrayView<2, T3, S3> dest);

    }
    \endcode

    \deprecatedAPI{fastCrossCorrelation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class MaskIterator, class MaskAccessor,
                  class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void fastCrossCorrelation(MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc,
                                  SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,  
                                  DestIterator d_ul, DestAccessor d_acc);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class MaskIterator, class MaskAccessor,
                  class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        fastCrossCorrelation(triple<MaskIterator, MaskIterator, MaskAccessor> mask,
                             triple<SrcIterator, SrcIterator, SrcAccessor> src,
                             pair<DestIterator, DestAccessor> dest);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/correlation.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int m_w=51, m_h=51;
    unsigned int w=1000, h=1000;
    MultiArray<2, float> mask(m_w,m_h), src(w,h), dest(w,h);
    ...
    
    //compute fast cross correlation of mask and image -> dest
    fastCrossCorrelation(mask, src, dest);
    \endcode
    
    <b> Preconditions:</b>

    The image must be larger than the size of the mask.
*/

doxygen_overloaded_function(template <...> void fastCrossCorrelation)
template <class MaskIterator, class MaskAccessor,
          class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void fastCrossCorrelation(MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc,
                          SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,  
                          DestIterator d_ul, DestAccessor d_acc)
{
    
    using namespace vigra;
    
    typedef typename MaskIterator::difference_type MaskDiffType;
    typedef typename SrcIterator::difference_type  SrcDiffType;
    typedef typename DestIterator::difference_type DestDiffType;
    
    MaskDiffType    mask_shape  = m_lr - m_ul;
    SrcDiffType     image_shape = s_lr - s_ul;
    
    const int   s_w = image_shape.x,    //Image dimensions
                s_h = image_shape.y,    //--
    
                m_w = mask_shape.x,     //Mask dimensions
                m_h = mask_shape.y,     //--
                m_total = m_w*m_h,
    
                padded_w = s_w+m_w,     //Padded dimensions
                padded_h = s_h+m_h,     //--
                padded_total = padded_w*padded_h;
    
    vigra_precondition( m_w % 2 == 1 , "vigra::fastCrossCorrelation(): Mask width has to be of odd size!");
    vigra_precondition( m_h % 2 == 1 , "vigra::fastCrossCorrelation(): Mask height has to be of odd size!");
    
    vigra_precondition( m_w <= s_w && m_h <= s_h , "vigra::fastCrossCorrelation(): Mask is larger than image!");
    
    //padding images to enable border cases
    FImage   padded_mask(padded_w, padded_h),
             padded_image(padded_w, padded_h);
  
    
    //fill padded mask  
    copyImage(srcIterRange(m_ul, m_lr, m_acc), destImage(padded_mask));
    
    //fill padded image
    copyImage(srcIterRange(s_ul, s_lr, s_acc), destImage(padded_image));
    
    
    //create fourier images
    FFTWComplexImage fourier_image(padded_w,padded_h),
    fourier_mask(padded_w,padded_h);
    
    //transform both
    fourierTransform(srcImageRange(padded_image), destImage(fourier_image));
    fourierTransform(srcImageRange(padded_mask), destImage(fourier_mask));
    
    detail::FourierCorrelationFunctor corr(padded_total);
    
    //do the correlation in frequency space
    combineTwoImages(srcImageRange(fourier_image), srcImage(fourier_mask), 
                     destImage(fourier_image), 
                     corr);
    
    //and go back to normal space
    fourierTransformInverse(srcImageRange(fourier_image), destImage(fourier_mask));
    
    //calculate the result, use the sum tables for the denumerators
    for(int v=0; v<=s_h-m_h; v++)
    {
        for(int u=0; u<=s_w-m_w; u++)
        {
            d_acc.set( (fourier_mask(u,v)).re(), d_ul + DestDiffType(u+m_w/2,v+m_h/2));
        }
    }
}

template <class MaskIterator, class MaskAccessor,
          class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void fastCrossCorrelation(triple<MaskIterator, MaskIterator, MaskAccessor> mask,
                                 triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                 std::pair<DestIterator, DestAccessor> dest)
{
    fastCrossCorrelation(mask.first, mask.second, mask.third, 
                         src.first, src.second, src.third, 
                         dest.first, dest.second);
}

template <class T1, class S1,
          class T2, class S2, 
          class T3, class S3>
inline void fastCrossCorrelation(MultiArrayView<2, T1, S1> const & mask,
                                 MultiArrayView<2, T2, S2> const & src,
                                 MultiArrayView<2, T3, S3> dest)
{   
    vigra_precondition(src.shape() == dest.shape(),
                        "vigra::fastCrossCorrelation(): shape mismatch between input and output.");
    fastCrossCorrelation(srcImageRange(mask), 
                         srcImageRange(src), 
                         destImage(dest));
}



/********************************************************/
/*                                                      */
/*          Cross correlation of mask to image          */
/*                                                      */
/********************************************************/
/**
    This function performes a (slow) cross-correlation using the window function environment
    and comparison of the mask with the underlying image part for each pixel. This may 
    however be faster for very few comparisons.
 */
 
/** \brief This function performes a (slow) cross-correlation

    This function performes a (slow) cross-correlation using the window function environment
    and comparison of the mask with the underlying image part for each pixel. This may 
    however be faster for very few comparisons.
    
    The input pixel type <tt>T1</tt> must be a \ref LinearSpace "linear space" over 
    the window functions' value_type <tt>T</tt>, i.e. addition of source values, multiplication with functions' values,
    and NumericTraits must be defined. The mask's value_type must be an \ref AlgebraicField "algebraic field",
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must be defined. 
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class T3, class S3>
        void
        crossCorrelation(MultiArrayView<2, T1, S1> const & mask,
                         MultiArrayView<2, T2, S2> const & src,
                         MultiArrayView<2, T3, S3> dest);

    }
    \endcode

    \deprecatedAPI{crossCorrelation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class MaskIterator, class MaskAccessor,
                  class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void crossCorrelation(MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc,
                              SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,  
                              DestIterator d_ul, DestAccessor d_acc);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class MaskIterator, class MaskAccessor,
                  class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        crossCorrelation(triple<MaskIterator, MaskIterator, MaskAccessor> mask,
                         triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/correlation.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int m_w=51, m_h=51;
    unsigned int w=1000, h=1000;
    MultiArray<2, float> mask(m_w,m_h), src(w,h), dest(w,h);
    ...
    
    //compute (slow) cross correlation of mask and image -> dest
    crossCorrelation(mask, src, dest);
    \endcode
    
    <b> Preconditions:</b>

    The image must be larger than the size of the mask.
*/
template<class MaskIterator, class MaskAccessor>
class CorrelationFunctor
{
public:
    CorrelationFunctor(MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc)
    : m_mask_ul(m_ul),
      m_mask_lr(m_lr),
      m_mask_acc(m_acc)
    {
    }
    
    CorrelationFunctor(triple<MaskIterator,MaskIterator,MaskAccessor> m)
    : m_mask_ul(m.first),
      m_mask_lr(m.second),
      m_mask_acc(m.third)
    {
    }
    
    template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
    void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc) const
    {   
        using namespace vigra;
        
        SrcIterator s_ul = s - windowShape()/2,
                    s_lr = s + windowShape()/2+Diff2D(1,1);
        
        //find img2 mean
        FindAverage<double> average;   
        inspectImage(srcIterRange(s_ul, s_lr, s_acc), average);
        
        SrcIterator ys = s_ul;
        SrcIterator xs = ys;
        
        MaskIterator ym = m_mask_ul;
        MaskIterator xm = ym;
        
        double  res=0;      
        
        for( ; ys.y != s_lr.y; ys.y++, ym.y++)
        {   
            for(xs = ys, xm = ym; xs.x != s_lr.x; xs.x++, xm.x++)
            {
                res += m_mask_acc(xm)*s_acc(xs);
            }       
        }
        d_acc.set(res,d);
    }
    
    Diff2D windowShape() const
    {
        return m_mask_lr - m_mask_ul;
    }
    
private:
    MaskIterator m_mask_ul;
    MaskIterator m_mask_lr;
    MaskAccessor m_mask_acc;    
};


template <class MaskIterator, class MaskAccessor, 
          class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void crossCorrelation(MaskIterator m_ul,  MaskIterator m_lr,   MaskAccessor m_acc,
                             SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                             DestIterator d_ul, DestAccessor d_acc,
                             BorderTreatmentMode border = BORDER_TREATMENT_AVOID)
{
    CorrelationFunctor<MaskIterator, MaskAccessor> func(m_ul, m_lr, m_acc);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class MaskIterator, class MaskAccessor, 
          class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void crossCorrelation(triple<MaskIterator, MaskIterator, MaskAccessor> m,
                             triple<SrcIterator, SrcIterator, SrcAccessor> s,
                             pair<DestIterator, DestAccessor> d, 
                             BorderTreatmentMode border = BORDER_TREATMENT_AVOID)
{
    crossCorrelation(m.first, m.second, m.third,
                     s.first, s.second, s.third,
                     d.first, d.second, 
                     border);
}

template <class T1, class S1,
          class T2, class S2, 
          class T3, class S3>
inline void crossCorrelation(MultiArrayView<2, T1, S1> const & mask,
                             MultiArrayView<2, T2, S2> const & src,
                             MultiArrayView<2, T3, S3> dest)
{   
    vigra_precondition(src.shape() == dest.shape(),
                        "vigra::crossCorrelation(): shape mismatch between input and output.");
    crossCorrelation(srcImageRange(mask), 
                     srcImageRange(src), 
                     destImage(dest));
}




/********************************************************/
/*                                                      */
/*     Normalized Cross correlation of mask to image    */
/*                                                      */
/********************************************************/
/**
    This function performes a (slow) normalized cross-correlation using the window function
    environment and comparison of the mask with the underlying image part for each pixel. 
    This may however be faster for very few comparisons.
 */
 
/** \brief This function performes a (slow) normalized cross-correlation

    This function performes a (slow) normalized cross-correlation using the window function
    environment and comparison of the mask with the underlying image part for each pixel. 
    This may however be faster for very few comparisons.
    
    The input pixel type <tt>T1</tt> must be a \ref LinearSpace "linear space" over 
    the window functions' value_type <tt>T</tt>, i.e. addition of source values, multiplication with functions' values,
    and NumericTraits must be defined. The mask's value_type must be an \ref AlgebraicField "algebraic field",
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must be defined. 
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class T3, class S3>
        void
        normalizedCrossCorrelation(MultiArrayView<2, T1, S1> const & mask,
                                   MultiArrayView<2, T2, S2> const & src,
                                   MultiArrayView<2, T3, S3> dest);

    }
    \endcode

    \deprecatedAPI{normalizedCrossCorrelation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class MaskIterator, class MaskAccessor,
                  class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void normalizedCrossCorrelation(MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc,
                                        SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,  
                                        DestIterator d_ul, DestAccessor d_acc);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class MaskIterator, class MaskAccessor,
                  class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        normalizedCrossCorrelation(triple<MaskIterator, MaskIterator, MaskAccessor> mask,
                                   triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                   pair<DestIterator, DestAccessor> dest);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/correlation.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int m_w=51, m_h=51;
    unsigned int w=1000, h=1000;
    MultiArray<2, float> mask(m_w,m_h), src(w,h), dest(w,h);
    ...
    
    //compute (slow) normalized cross correlation of mask and image -> dest
    normalizedCrossCorrelation(mask, src, dest);
    \endcode
    
    <b> Preconditions:</b>

    The image must be larger than the size of the mask.
*/
template<class MaskIterator, class MaskAccessor>
class NormalizedCorrelationFunctor
{
public:
    
    NormalizedCorrelationFunctor(MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc)
    : m_mask_ul(m_ul),
      m_mask_lr(m_lr),
      m_mask_acc(m_acc),
      m_s11(0.0),
      m_avg1(0.0)
    {
        init_s11();
    }
    
    NormalizedCorrelationFunctor(triple<MaskIterator,MaskIterator,MaskAccessor> m)
    : m_mask_ul(m.first),
      m_mask_lr(m.second),
      m_mask_acc(m.third),
      m_s11(0.0),
      m_avg1(0.0)
    {
        init_s11();
    }
    
    template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
    void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc) const
    {   
        using namespace vigra;
        
        SrcIterator s_ul = s - windowShape()/2,
        s_lr = s + windowShape()/2+Diff2D(1,1);
        
        if(m_s11 == 0.0)
        {
            d_acc.set(0,d);
        }
        else
        {
            //find img2 mean
            FindAverage<double> average;   
            inspectImage(srcIterRange(s_ul, s_lr, s_acc), average);
            
            SrcIterator ys = s_ul;
            SrcIterator xs = ys;
            
            MaskIterator ym = m_mask_ul;
            MaskIterator xm = ym;
            
            double  s1=0,s2=0, s12=0, s22=0;        
            
            for( ; ys.y != s_lr.y; ys.y++, ym.y++)
            {   
                for(xs = ys, xm = ym; xs.x != s_lr.x; xs.x++, xm.x++)
                {
                    s1 = m_mask_acc(xm);
                    s2 = s_acc(xs);
                    s12 += (s1-m_avg1)*(s2-average());
                    s22 += (s2-average())*(s2-average());
                }       
            }
            if(s22 == 0.0)
            {
                d_acc.set(0,d);
            }
            else
            {
                d_acc.set(s12/sqrt(m_s11*s22),d);
            }
        }
    }
    
    Diff2D windowShape() const
    {
        return m_mask_lr - m_mask_ul;
    }
    
private:
    void init_s11()
    {
        //find mask mean
        FindAverage<double> average;   
        inspectImage(srcIterRange(m_mask_ul, m_mask_lr, m_mask_acc), average);
        
        MaskIterator ym = m_mask_ul;
        MaskIterator xm = ym;
        
        m_avg1 = average();
        
        for( ; ym.y != m_mask_lr.y; ym.y++)
        {   
            for(xm = ym; xm.x != m_mask_lr.x; xm.x++)
            {
                m_s11 += (m_mask_acc(xm)-m_avg1)*(m_mask_acc(xm)-m_avg1);
            }
        }
    }
    
    MaskIterator m_mask_ul;
    MaskIterator m_mask_lr;
    MaskAccessor m_mask_acc;
    
    double m_s11;
    double m_avg1;
};


template <class MaskIterator, class MaskAccessor, 
          class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void normalizedCrossCorrelation(MaskIterator m_ul,  MaskIterator m_lr,   MaskAccessor m_acc,
                                       SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                                       DestIterator d_ul, DestAccessor d_acc,
                                       BorderTreatmentMode border = BORDER_TREATMENT_AVOID)
{
    NormalizedCorrelationFunctor<MaskIterator, MaskAccessor> func(m_ul, m_lr, m_acc);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class MaskIterator, class MaskAccessor, 
          class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void normalizedCrossCorrelation(triple<MaskIterator, MaskIterator, MaskAccessor> m,
                                       triple<SrcIterator, SrcIterator, SrcAccessor> s,
                                       pair<DestIterator, DestAccessor> d, 
                                       BorderTreatmentMode border = BORDER_TREATMENT_AVOID)
{
    normalizedCrossCorrelation(m.first, m.second, m.third,
                               s.first, s.second, s.third,
                               d.first, d.second, 
                               border);
}
    
template <class T1, class S1,
          class T2, class S2, 
          class T3, class S3>
inline void normalizedCrossCorrelation(MultiArrayView<2, T1, S1> const & mask,
                                       MultiArrayView<2, T2, S2> const & src,
                                       MultiArrayView<2, T3, S3> dest)
{   
    vigra_precondition(src.shape() == dest.shape(),
                        "vigra::normalizedCrossCorrelation(): shape mismatch between input and output.");
    normalizedCrossCorrelation(srcImageRange(mask), 
                               srcImageRange(src), 
                               destImage(dest));
}
    
//@}

} //end of namespace vigra

#endif //VIGRA_CORRELATION_HXX