/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2014 by Benjamin Seppke                 */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_MEDIANFILTER_HXX
#define VIGRA_MEDIANFILTER_HXX

#include <vector>
#include <algorithm>

#include "applywindowfunction.hxx"

namespace vigra
{

/********************************************************/
/*                                                      */
/*              Generic median filter                   */
/*                                                      */
/********************************************************/
/**  
    This function calculates the median of a window of given size for the complete image. 
    It also allows a correct border handling, since it uses the \ref applyWindowFunction 
    environment for computation! 
*/
//@{

/** \brief This function calculates the median of a window of given size for the complete image.

    All \ref BorderTreatmentMode "border treatment modes"  (except BORDER_TREATMENT_CLIP)  are supported.

    The input pixel type <tt>T1</tt> must be a \ref LinearSpace "linear space" over 
    the window functions' value_type <tt>T</tt>. Especially, the values must be sortable by
    std::sort, to derive the mean values aka the median.
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        medianFilter(MultiArrayView<2, T1, S1> const & src,
                     MultiArrayView<2, T2, S2> dest,
                     Diff2D window_shape, 
                     BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);

    }
    \endcode

    \deprecatedAPI{medianFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void medianFilter(SrcIterator supperleft,
                          SrcIterator slowerright, SrcAccessor sa,
                          DestIterator dupperleft, DestAccessor da,
                          Diff2D window_shape, 
                          BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        medianFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                     pair<DestIterator, DestAccessor> dest,
                     Diff2D window_shape, 
                     BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/medianfilter.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int w=1000, h=1000;
    MultiArray<2, float> src(w,h), dest(w,h);
    ...
    
    // apply a median filter with a window size of 5x5
    medianFilter(src, dest, Diff2D(5,5));
    \endcode
    
    <b> Preconditions:</b>

    The image must be larger than the window size of the filter.
*/

doxygen_overloaded_function(template <...> void medianFilter)

template<class VALUETYPE>
class MedianFunctor
{
public:
    MedianFunctor(Diff2D window_shape)
    : m_window_shape(window_shape),
      m_buffer(window_shape.x*window_shape.y)
    {
    }

    template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
    void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc)
    {
        SrcIterator s_ul = s - m_window_shape/2,
                    s_lr = s_ul + m_window_shape;
        
        std::fill(m_buffer.begin(), m_buffer.end(), VALUETYPE());
        
        SrcIterator ys = s_ul;
        SrcIterator xs = ys;
        
        typename std::vector<VALUETYPE>::iterator iter = m_buffer.begin();
        
        for( ; ys.y != s_lr.y; ys.y++)
        {   
            for(xs = ys; xs.x != s_lr.x; xs.x++, iter++)
            {
                *iter = s_acc(xs);
            }       
        }
        
        std::sort(m_buffer.begin(), m_buffer.end());
        d_acc.set(m_buffer[m_buffer.size()/2],d);
    }
    
    Diff2D windowShape() const
    {
        return m_window_shape;
    }
    
private:
    Diff2D m_window_shape;
    std::vector<VALUETYPE> m_buffer;    
};


template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void medianFilter(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                         DestIterator d_ul, DestAccessor d_acc, 
                         Diff2D window_shape,
                         BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    MedianFunctor<typename SrcIterator::value_type> func(window_shape);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void medianFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                         pair<DestIterator, DestAccessor> d, 
                         Diff2D window_shape,
                         BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    medianFilter(s.first, s.second, s.third,
                 d.first, d.second, 
                 window_shape,
                 border);
}

template <class T1, class S1, 
          class T2, class S2>
inline void medianFilter(MultiArrayView<2, T1, S1> const & src,
                         MultiArrayView<2, T2, S2> dest, 
                         Diff2D window_shape,
                         BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    vigra_precondition(src.shape() == dest.shape(),
                        "vigra::medianFilter(): shape mismatch between input and output.");
    medianFilter(srcImageRange(src),
                 destImage(dest), 
                 window_shape, 
                 border);
}

//@}

} //end of namespace vigra

#endif //VIGRA_MEDIANFILTER_HXX
