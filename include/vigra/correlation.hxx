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
#include "multi_math.hxx"
#include "multi_fft.hxx"

// "slow" correlation algorithms are performed using the windowing filter env.
#include "applywindowfunction.hxx"


namespace vigra
{
    
    /********************************************************/
    /*                                                      */
    /*    Fast cross correlation of an image w.r.t a mask   */
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
     
     By default, the borders are filled with zeros. Use the clearBorders switch to change that behavior if you need to.
     
     <b> Declarations:</b>
     
     pass 2D array views:
     \code
     namespace vigra {
       template <class T1, class S1,
                class T2, class S2,
                class T3, class S3>
       void
       fastCrossCorrelation(MultiArrayView<2, T1, S1> const & in,
                            MultiArrayView<2, T2, S2> const & mask,
                            MultiArrayView<2, T3, S3> out,
                            bool clearBorders=true);
     
     }
     \endcode
     
     
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
    template <class T1, class S1,
    class T2, class S2,
    class T3, class S3>
    inline void fastCrossCorrelation(MultiArrayView<2, T1, S1> const & in,
                                     MultiArrayView<2, T2, S2> const & mask,
                                     MultiArrayView<2, T3, S3> out,
                                     bool clearBorders=true)
    {
        vigra_precondition(in.shape() == out.shape(),
                           "vigra::fastCrossCorrelation(): shape mismatch between input and output.");
        correlateFFT(in, mask, out);
        
        if(clearBorders)
        {
            Shape2 border(mask.width()/2, mask.height()/2);
            initMultiArrayBorder(out, border, border, T3());
        }
    }
    
    

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
    
    By default, the borders are filled with zeros. Use the clearBorders switch to change that behavior if you need to.
 
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class T3, class S3>
        void
        fastNormalizedCrossCorrelation(MultiArrayView<2, T1, S1> const & in,
                                       MultiArrayView<2, T2, S2> const & mask,
                                       MultiArrayView<2, T3, S3> out,
                                       bool clearBorders=true);

    }
    \endcode


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

template <class T1, class S1,
          class T2, class S2, 
          class T3, class S3>
inline void fastNormalizedCrossCorrelation(MultiArrayView<2, T1, S1> const & in,
                                           MultiArrayView<2, T2, S2> const & mask,
                                           MultiArrayView<2, T3, S3> out,
                                           bool clearBorders=true)
{
    
    using namespace vigra::multi_math;
    
    vigra_precondition(in.shape() == out.shape(),
                        "vigra::fastNormalizedCrossCorrelation(): shape mismatch between input and output.");

    unsigned int m_w = mask.width(),
                 m_h = mask.height(),
                 m_total = m_w*m_h,
                 i_w = in.width(),
                 i_h = in.height();
    
    vigra_precondition( m_w % 2 == 1 , "vigra::fastNormalizedCrossCorrelation(): Mask width has to be of odd size!");
    vigra_precondition( m_h % 2 == 1 , "vigra::fastNormalizedCrossCorrelation(): Mask height has to be of odd size!");
    
    vigra_precondition( m_w <= i_w && m_h <= i_h , "vigra::fastNormalizedCrossCorrelation(): Mask is larger than image!");
    
    //find mask sum and mask^2 sum and mask average
    double  mask_sum  = 0.0,
            mask_sum2 = 0.0,
            mask_avg  = 0.0;
    
    for(unsigned int y=0; y<m_h ;++y)
    {
        for(unsigned int x=0; x<m_w ;++x)
        {
            //cache current value
            mask_avg = mask(x,y);
            
            mask_sum  += mask_avg;
            mask_sum2 += mask_avg*mask_avg;
        }
    }
    
    mask_avg = mask_sum/(m_total);
    
    //calculate the fix part of the denumerator
    double fix_denumerator = sqrt(m_total*mask_sum2 - mask_sum*mask_sum);
    
    if(fix_denumerator == 0)
    {
        out = 0;
    }
    else
    {
        //pre-normalize the mask
        MultiArray<2, T2> norm_mask(m_w,m_h);
        norm_mask = mask - mask_avg;
        
        //calculate (semi-normalized) numerator:
        MultiArray<2, T3> corr_result(i_w,i_h);
        fastCrossCorrelation(in, norm_mask, corr_result,clearBorders);
        
        //Create fast sum tables for the variable denumerator
        MultiArray<2, double> sum_table(i_w+1,i_h+1),
                              sum_table2(i_w+1,i_h+1);
        
        for(unsigned int v=0; v<i_h; v++)
        {
            for(unsigned int u=0; u<i_w; u++)
            {
                sum_table( u+1,v+1)  = in(u,v)         +  sum_table(u,v+1)  + sum_table(u+1,v)  -  sum_table(u,v);
                sum_table2(u+1,v+1)  = in(u,v)*in(u,v) + sum_table2(u,v+1) + sum_table2(u+1,v)  - sum_table2(u,v);
            }
        }
        //calculate the result, use the sum tables for the denumerators
        for(unsigned int v=m_h/2; v<i_h-m_h/2; v++)
        {
            for(unsigned int u=m_w/2; u<i_w-m_w/2; u++)
            {
                //calculate e(window) and e(window^2)
                double e_uv   = sum_table( u+m_w/2+1, v+m_h/2+1) - sum_table( u-m_w/2, v+m_h/2+1) - sum_table( u+m_w/2+1, v-m_h/2)  + sum_table( u-m_w/2,v-m_h/2),
                       e_uv_2 = sum_table2(u+m_w/2+1, v+m_h/2+1) - sum_table2(u-m_w/2, v+m_h/2+1) - sum_table2(u+m_w/2+1, v-m_h/2)  + sum_table2(u-m_w/2,v-m_h/2),
                       var_denumerator = sqrt(m_total*e_uv_2 - e_uv*e_uv);
                
                //calclate overall result
                if(var_denumerator == 0)
                {
                    out(u,v) = 0;
                }
                else
                {
                    out(u,v) = m_total*corr_result(u,v)/(var_denumerator*fix_denumerator);
                }
            }
        }
    }
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
        crossCorrelation(MultiArrayView<2, T1, S1> const & in,
                         MultiArrayView<2, T2, S2> const & mask,
                         MultiArrayView<2, T3, S3> out);

    }
    \endcode

    \deprecatedAPI{crossCorrelation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void crossCorrelation(SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,
                              MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc,
                              DestIterator d_ul, DestAccessor d_acc);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void
        crossCorrelation(triple<MaskIterator, MaskIterator, MaskAccessor> src,
                         triple<SrcIterator, SrcIterator, SrcAccessor> mask,
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


template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void crossCorrelation(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                             MaskIterator m_ul,  MaskIterator m_lr,   MaskAccessor m_acc,
                             DestIterator d_ul, DestAccessor d_acc,
                             BorderTreatmentMode border = BORDER_TREATMENT_AVOID)
{
    CorrelationFunctor<MaskIterator, MaskAccessor> func(m_ul, m_lr, m_acc);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void crossCorrelation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                             triple<MaskIterator, MaskIterator, MaskAccessor> mask,
                             pair<DestIterator, DestAccessor> dest,
                             BorderTreatmentMode border = BORDER_TREATMENT_AVOID)
{
    crossCorrelation(src.first, src.second, src.third,
                     mask.first, mask.second, mask.third,
                     dest.first, dest.second,
                     border);
}

template <class T1, class S1,
          class T2, class S2, 
          class T3, class S3>
inline void crossCorrelation(MultiArrayView<2, T1, S1> const & in,
                             MultiArrayView<2, T2, S2> const & mask,
                             MultiArrayView<2, T3, S3> out)
{   
    vigra_precondition(in.shape() == out.shape(),
                        "vigra::crossCorrelation(): shape mismatch between input and output.");
    crossCorrelation(srcImageRange(in),
                     srcImageRange(mask),
                     destImage(out));
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
        normalizedCrossCorrelation(MultiArrayView<2, T1, S1> const & in,
                                   MultiArrayView<2, T2, S2> const & mask,
                                   MultiArrayView<2, T3, S3> out);

    }
    \endcode

    \deprecatedAPI{normalizedCrossCorrelation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void normalizedCrossCorrelation(SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,
                                        MaskIterator m_ul, MaskIterator m_lr, MaskAccessor m_acc,
                                        DestIterator d_ul, DestAccessor d_acc);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void
        normalizedCrossCorrelation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                   triple<MaskIterator, MaskIterator, MaskAccessor> mask,
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
    
    NormalizedCorrelationFunctor(triple<MaskIterator,MaskIterator,MaskAccessor> mask)
    : m_mask_ul(mask.first),
      m_mask_lr(mask.second),
      m_mask_acc(mask.third),
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


template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void normalizedCrossCorrelation(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                                       MaskIterator m_ul,  MaskIterator m_lr,   MaskAccessor m_acc,
                                       DestIterator d_ul, DestAccessor d_acc,
                                       BorderTreatmentMode border = BORDER_TREATMENT_AVOID)
{
    NormalizedCorrelationFunctor<MaskIterator, MaskAccessor> func(m_ul, m_lr, m_acc);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void normalizedCrossCorrelation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                       triple<MaskIterator, MaskIterator, MaskAccessor> mask,
                                       pair<DestIterator, DestAccessor> dest,
                                       BorderTreatmentMode border = BORDER_TREATMENT_AVOID)
{
    normalizedCrossCorrelation(src.first, src.second, src.third,
                               mask.first, mask.second, mask.third,
                               dest.first, dest.second,
                               border);
}
    
template <class T1, class S1,
          class T2, class S2, 
          class T3, class S3>
inline void normalizedCrossCorrelation(MultiArrayView<2, T1, S1> const & in,
                                       MultiArrayView<2, T2, S2> const & mask,
                                       MultiArrayView<2, T3, S3> out)
{   
    vigra_precondition(in.shape() == out.shape(),
                        "vigra::normalizedCrossCorrelation(): shape mismatch between input and output.");
    normalizedCrossCorrelation(srcImageRange(in),
                               srcImageRange(mask),
                               destImage(out));
}
    
//@}

} //end of namespace vigra

#endif //VIGRA_CORRELATION_HXX