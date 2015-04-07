/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2014 by Benjamin Seppke                 */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_APPLYWINDOWFUCTION_HXX
#define VIGRA_APPLYWINDOWFUCTION_HXX

#include "basicimage.hxx"
#include "copyimage.hxx"
#include "basicgeometry.hxx"
#include "initimage.hxx"
#include "bordertreatment.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*             Apply window filters to images           */
/*                                                      */
/********************************************************/

/**  
    This function calculates the results for a window function (given as a functor) when
    applied to the complete image. Also allows a correct border handling! 
    See \ref medianFilter() for an example of a quite basic window function and its application.
*/
//@{

/** \brief Apply a window function to each pixels of a given image.

    If you pass a functor to this function, which implements the two functions:
    <ol>
        <li>Diff2D windowShape() const<br/>
        to return the filter window size, which has to be odd in each dimension and</li>
        <li>void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc)<br/>
        to compute the results of the current window,</li>
    </ol>
    this function calculates the results for the complete image.

    All \ref BorderTreatmentMode "border treatment modes"  (except BORDER_TREATMENT_CLIP)  are supported.

    The input pixel type <tt>T1</tt> must be a \ref LinearSpace "linear space" over 
    the window functions' value_type <tt>T</tt>, i.e. addition of source values, multiplication with functions' values,
    and NumericTraits must be defined. The filters' value_type must be an \ref AlgebraicField "algebraic field",
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must be defined. 
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2,
                  class ProcessingFunctor>
        void
        applyWindowFunction(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            ProcessingFunctor func, 
                            BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);

    }
    \endcode

    \deprecatedAPI{applyWindowFunction}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class ProcessingFunctor>
        void applyWindowFunction(SrcIterator supperleft,
                                 SrcIterator slowerright, SrcAccessor sa,
                                 DestIterator dupperleft, DestAccessor da,
                                 ProcessingFunctor func, 
                                 BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class ProcessingFunctor>
        void
        applyWindowFunction(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            ProcessingFunctor func, 
                            BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/applywindowfunction.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int w=1000, h=1000;
    MultiArray<2, float> src(w,h), dest(w,h);
    ...
    
    template<class VALUETYPE>
    class AvgFunctor
    {
        public:
            MedianFunctor(Diff2D window_shape)
            : m_window_shape(window_shape)
            {
            }
            
            template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
            void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc)
            {
                SrcIterator s_ul = s - m_window_shape/2,
                            s_lr = s_ul + m_window_shape;
        
                VALUETYPE result = NumericTraits<float>::zero();
                                
                SrcIterator ys = s_ul;
                SrcIterator xs = ys;
            
                for( ; ys.y != s_lr.y; ys.y++)
                {   
                    for(xs = ys; xs.x != s_lr.x; xs.x++, iter++)
                    {
                        res += s_acc(xs);
                    }       
                }
                
                d_acc.set(res/(m_window_shape.x*m_window_shape.y),d);
            }
    
            Diff2D windowShape() const
            {
                return m_window_shape;
            }
        private:
            Diff2D m_window_shape;  
    };
 
    
    // create an AverageFilter function for a 5x5 filter
    AvgFunctor func(Diff2D(5,5));
    
    
    // apply the filter function to the input image
    applyWindowFunction(src, dest, func);
    \endcode

    <b> Preconditions:</b>

    The image must be larger than the window size.
*/

doxygen_overloaded_function(template <...> void applyWindowFunction)

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, 
          class ProcessingFunctor>
void applyWindowFunction(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                         DestIterator d_ul, DestAccessor d_acc, 
                         ProcessingFunctor func, 
                         BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    vigra_precondition((border == BORDER_TREATMENT_AVOID   ||
                      //border == BORDER_TREATMENT_CLIP    ||
                        border == BORDER_TREATMENT_REPEAT  ||
                        border == BORDER_TREATMENT_REFLECT ||
                        border == BORDER_TREATMENT_WRAP    ||
                        border == BORDER_TREATMENT_ZEROPAD   
                        ),  
                       "vigra::applyWindowFunction():\n"
                       "  Border treatment must be one of follow treatments:\n"
                       "  - BORDER_TREATMENT_AVOID\n"
                       //"  - BORDER_TREATMENT_CLIP\n" 
                       "  - BORDER_TREATMENT_REPEAT\n"
                       "  - BORDER_TREATMENT_REFLECT\n"
                       "  - BORDER_TREATMENT_WRAP\n"
                       "  - BORDER_TREATMENT_ZEROPAD\n"
                       );
    
    typename SrcIterator::difference_type img_shape = s_lr - s_ul;
    Diff2D win_shape = func.windowShape();
    
    vigra_precondition( win_shape.x % 2 == 1 , "vigra::applyWindowFunction(): Filter window width has to be of odd size!");
    vigra_precondition( win_shape.y % 2 == 1 , "vigra::applyWindowFunction(): Filter window height has to be of odd size!");
    
    vigra_precondition( win_shape.x <= img_shape.x && win_shape.y <= img_shape.y , "vigra::applyWindowFunction(): Filter window is larger than image!");
    
    /**********************************************************************************
     *                                                                                *
     *  COMPUTE ALL "SAFE" PIXELS, WHERE THE MASK FITS COMPLETELY INTO THE IMAGE      *
     *                                                                                *
     **********************************************************************************/
    
    
    SrcIterator     ys  = s_ul, 
                    xs  = ys;
    
    DestIterator    yd  = d_ul, 
                    xd  = yd;
    
    SrcIterator     end = s_ul + img_shape - win_shape/2;
    
    ys.y += win_shape.y/2;
    yd.y += win_shape.y/2;
    
    unsigned int y=0;
    
    for( ; ys.y != end.y; ys.y++, yd.y++, y++)
    {   
        xs = ys;
        xs.x += win_shape.x/2;
        
        xd = yd;
        xd.x += win_shape.x/2;
                
        for( ; xs.x != end.x; xs.x++, xd.x++)
        {
            func(xs, s_acc, xd, d_acc);
        }
    }
    
    
    
    /**********************************************************************************
     *                                                                                *
     *                      HANDLE THE EIGHT BORDER CASES SEPARATELY                  *
     *                                                                                *
     *                                                                                *
     *  and do this diffently according to the used border treatment type, of course  *
     *                                                                                *
     **********************************************************************************/
    
    
    if(border == BORDER_TREATMENT_AVOID)
        return; // skip processing near the border
    
    
    //Do some preparation for the special cases:
    //Create window of width = image width and height = win_shape
    BasicImage<typename SrcIterator::PixelType> temp(img_shape.x, win_shape.y + win_shape.y/2);
    typedef typename BasicImage<typename SrcIterator::PixelType>::Iterator TempIterator;
    
    TempIterator    t_ul  = temp.upperLeft(), 
                    t_lr  = temp.lowerRight();
    
    DestAccessor    t_acc = temp.accessor();
    
    TempIterator    yt = t_ul,
                    xt = yt;
    
    
    /**********************************************************************************
     *                                                                                * 
     *                             SIDE CASE 1 (1/8):                                 * 
     *                                                                                *                                                                               *
     *       HANDLE UPPER PIXELS WHERE THE MASK IS MISSING MINIMUM y-COORDINATES      *
     *                                                                                *
     **********************************************************************************/
    if(border == BORDER_TREATMENT_REPEAT)
    {
        ys = s_ul;
        Diff2D lineDiff(img_shape.x,1);
        
        for( ; yt.y != t_lr.y - win_shape.y; ++yt.y)
        {
            copyImage(ys, ys+lineDiff, s_acc, yt, t_acc);
        }
        copyImage(ys, ys+Diff2D(img_shape.x,win_shape.y), s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        reflectImage(s_ul, s_ul+Diff2D(img_shape.x,win_shape.y/2), s_acc, t_ul, t_acc, horizontal);
        copyImage(s_ul, s_ul+Diff2D(img_shape.x,win_shape.y), s_acc, t_ul+Diff2D(0,win_shape.y/2), t_acc);
        
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
        copyImage(s_ul+Diff2D(0, img_shape.y-win_shape.y/2), s_lr, s_acc, t_ul, t_acc);
        copyImage(s_ul, s_ul+Diff2D(img_shape.x,win_shape.y), s_acc, t_ul+Diff2D(0,win_shape.y/2), t_acc);
        
    }
    else if(border == BORDER_TREATMENT_ZEROPAD)
    {
        initImage(t_ul, t_lr, t_acc, 0);
        copyImage(s_ul, s_ul+Diff2D(img_shape.x,win_shape.y), s_acc, t_ul+Diff2D(0,win_shape.y/2), t_acc);
        
    } 
    
    yt = t_ul;
    yt.y += win_shape.y/2;
    yd = d_ul;
    
    for( ; yt.y != t_lr.y-win_shape.y/2; ++yd.y,  ++yt.y)
    {
        xt = yt;
        xt.x += win_shape.x/2;
        
        xd = yd;
        xd.x += win_shape.x/2;
        
        for( ; xt.x != t_lr.x-win_shape.x/2; xd.x++, xt.x++)
        {
            func(xt, t_acc, xd, d_acc);     
        }
    }
    
    /**********************************************************************************
     *                                                                                * 
     *                             SIDE CASE 2 (2/8):                                 * 
     *                                                                                *                                                                               *
     *       HANDLE LOWER PIXELS WHERE THE MASK IS MISSING MAXIMUM y-COORDINATES      *
     *                                                                                *
     **********************************************************************************/
    if(border == BORDER_TREATMENT_REPEAT)
    {
        ys = s_ul + Diff2D(0, img_shape.y-1);
        yt = t_ul + Diff2D(0, win_shape.x);
        
        Diff2D lineDiff(img_shape.x,1);
        
        for( ; yt.y != t_lr.y ; ++yt.y)
        {
            copyImage(ys, ys+lineDiff, s_acc, yt, t_acc);
        }
        ys = s_ul + Diff2D(0, img_shape.y-win_shape.y);
        yt = t_ul;
        copyImage(ys, ys+Diff2D(img_shape.x,win_shape.y), s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        reflectImage(s_ul+Diff2D(0, img_shape.y-win_shape.y/2), s_lr, s_acc, t_ul+Diff2D(0, win_shape.y), t_acc, horizontal);
        copyImage(s_ul+Diff2D(0, img_shape.y-win_shape.y), s_lr, s_acc, t_ul, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
        copyImage(s_ul, s_ul + Diff2D(img_shape.x, win_shape.y/2), s_acc, t_ul+Diff2D(0, win_shape.y), t_acc);
        copyImage(s_ul+Diff2D(0, img_shape.y-win_shape.y), s_lr, s_acc, t_ul, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_ZEROPAD)
    {
        initImage(t_ul, t_lr, t_acc, 0);        
        copyImage(s_ul+Diff2D(0, img_shape.y-win_shape.y), s_lr, s_acc, t_ul, t_acc);
        
    } 
    
    
    yt = t_ul;
    yt.y += win_shape.y/2;
    yd = d_ul;
    yd.y += img_shape.y-win_shape.y/2-1;
    
    for( ; yt.y != t_lr.y - win_shape.y/2 ; ++yd.y, ++yt.y)
    {
        xt = yt;
        xt.x += win_shape.x/2;
        
        xd = yd;
        xd.x += win_shape.x/2;
        
        for( ; xt.x != t_lr.x-win_shape.x/2; xd.x++, xt.x++)
        {
            func(xt, t_acc, xd, d_acc);                 
        }
    }


    //Preparation needed for left and right processing
    temp.resize(win_shape.x+win_shape.x/2,img_shape.y);
    t_ul = temp.upperLeft(); t_lr = temp.lowerRight();
    t_acc = temp.accessor();
    
    
    
    /**********************************************************************************
     *                                                                                * 
     *                             SIDE CASE 3 (3/8):                                 * 
     *                                                                                *
     *       HANDLE LEFT PIXELS WHERE THE MASK IS MISSING MINIMUM x-COORDINATES       *
     *                                                                                *
     **********************************************************************************/
    if(border == BORDER_TREATMENT_REPEAT)
    {
        xs = s_ul;
        xt = t_ul;
        xd = d_ul;
        
        Diff2D colDiff(1,img_shape.y);
        
        for( ; xt.x != t_lr.x - win_shape.x; ++xt.x)
        {
            copyImage(xs, xs+colDiff, s_acc, xt, t_acc);
        }
        copyImage(xs, xs+Diff2D(win_shape.x, img_shape.y), s_acc, xt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        reflectImage(s_ul, s_ul+Diff2D(win_shape.x/2, img_shape.y), s_acc, t_ul, t_acc, vertical);
        copyImage(s_ul, s_ul+Diff2D(win_shape.x, img_shape.y), s_acc, t_ul + Diff2D(win_shape.x/2, 0), t_acc);
        
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
        copyImage(s_ul+Diff2D(img_shape.x-win_shape.x/2, 0), s_lr, s_acc, t_ul, t_acc);
        copyImage(s_ul, s_ul+Diff2D(win_shape.x, img_shape.y), s_acc, t_ul + Diff2D(win_shape.x/2, 0), t_acc);
        
    }
    else if(border == BORDER_TREATMENT_ZEROPAD)
    {
        initImage(t_ul, t_lr, t_acc, 0);        
        copyImage(s_ul, s_ul+Diff2D(win_shape.x, img_shape.y), s_acc, t_ul + Diff2D(win_shape.x/2, 0), t_acc);
        
    } 
        
    
    yt = t_ul;
    yt.y += win_shape.y/2;
    yd = d_ul;
    yd.y += win_shape.y/2;
    
    for( ; yt.y != t_lr.y-win_shape.y/2; ++yd.y,  ++yt.y)
    {
        xt = yt;
        xt.x += win_shape.x/2;
        
        xd = yd;
        
        for( ; xt.x != t_lr.x-win_shape.x/2; xd.x++, xt.x++)
        {
            func(xt, t_acc, xd, d_acc);                 
        }
    }
    
    
    /**********************************************************************************
     *                                                                                * 
     *                             SIDE CASE 4 (4/8):                                 * 
     *                                                                                *
     *       HANDLE RIGHT PIXELS WHERE THE MASK IS MISSING MAXIMUM x-COORDINATES      *
     *                                                                                *
     **********************************************************************************/
    if(border == BORDER_TREATMENT_REPEAT)
    {
        xs = s_ul + Diff2D(img_shape.x-1,0);
        xt = t_ul + Diff2D(win_shape.x,0);
        xd = d_ul;
        
        Diff2D colDiff(1,img_shape.y);
        
        for( ; xt.x != t_lr.x ; ++xt.x)
        {
            copyImage(xs, xs+colDiff, s_acc, xt, t_acc);
        }
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,0);
        yt = t_ul;
        copyImage(ys, ys+Diff2D(win_shape.x,img_shape.y), s_acc, yt, t_acc);
        
        
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        reflectImage(s_ul+Diff2D(img_shape.x-win_shape.x/2,0), s_lr, s_acc, t_ul+Diff2D(win_shape.x,0), t_acc, vertical);
        copyImage(s_ul+Diff2D(img_shape.x-win_shape.x,0), s_lr, s_acc, t_ul, t_acc);
        
        
    }   
    else if(border == BORDER_TREATMENT_WRAP)
    {
        copyImage(s_ul, s_ul+Diff2D(win_shape.x/2,img_shape.y), s_acc, t_ul+Diff2D(win_shape.x,0), t_acc);
        copyImage(s_ul+Diff2D(img_shape.x-win_shape.x,0), s_lr, s_acc, t_ul, t_acc);
        
        
    }   
    else if(border == BORDER_TREATMENT_ZEROPAD)
    {
        initImage(t_ul, t_lr, t_acc, 0);        
        copyImage(s_ul+Diff2D(img_shape.x-win_shape.x,0), s_lr, s_acc, t_ul, t_acc);
        
    } 
    
    
    yt = t_ul;
    yt.y += win_shape.x/2;
    yd = d_ul;
    yd.x += img_shape.x-win_shape.x/2-1;
    yd.y += win_shape.y/2;
    
    for( ; yt.y != t_lr.y-win_shape.y/2; ++yd.y,  ++yt.y)
    {
        xt = yt;
        xt.x += win_shape.x/2;
        
        xd = yd;
        
        for( ; xt.x != t_lr.x-win_shape.x/2; xd.x++, xt.x++)
        {
            func(xt, t_acc, xd, d_acc);             
        }
    }
    
    //Do some preaparations for the corner cases
    temp.resize(win_shape+win_shape/2);
    t_ul = temp.upperLeft(); t_lr = temp.lowerRight();
    t_acc = temp.accessor();
    
    
    
    /**********************************************************************************
     *                                                                                * 
     *                             CORNER CASE 1 (5/8):                               * 
     *                                                                                *
     *            HANDLE UPPERLEFT PIXELS WHERE THE MASK IS MISSING MINIMUM           *
     *                 x-COORDINATES AND MINIMUM y-COORDINATES                        *
     *                                                                                *
     **********************************************************************************/
    if(border == BORDER_TREATMENT_REPEAT)
    {
        //init upperleft rect with single value
        ys = s_ul;
        yt = t_ul;
        initImage(yt,yt+win_shape/2, t_acc, s_acc(ys));
        
        //init upperright rect with vertical stripes
        ys = s_ul;
        yt = t_ul + Diff2D(win_shape.x/2,0);
        Diff2D lineDiff(win_shape.x,1);
        for( ; yt.y != t_lr.y-win_shape.y ; ++yt.y)
        {
            copyImage(ys, ys+lineDiff, s_acc, yt, t_acc);
        }
        
        //init lowerleft rect with horizontal stripes
        xs = s_ul;
        xt = t_ul + Diff2D(0,win_shape.y/2);
        Diff2D rowDiff(1, win_shape.y);
        for( ; xt.x != t_lr.x-win_shape.x ; ++xt.x)
        {
            copyImage(xs, xs+rowDiff, s_acc, xt, t_acc);
        }
        
        //copy image patch in lower right patch
        ys = s_ul;
        yt = t_ul + win_shape/2;
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        //init upperleft rect with double reflect image
        ys = s_ul;
        yt = t_ul;
        rotateImage(ys,ys+win_shape/2, s_acc, yt, t_acc, 180);
        
        //init upperright rect with horizontal reflected image
        ys = s_ul;
        yt = t_ul + Diff2D(win_shape.x/2,0);
        reflectImage(ys, ys+Diff2D(win_shape.x, win_shape.y/2), s_acc, yt, t_acc, horizontal);
        
        //init lowerleft rect with vertical reflected image
        xs = s_ul;
        xt = t_ul + Diff2D(0,win_shape.y/2);
        reflectImage(xs, xs+Diff2D(win_shape.x/2, win_shape.y), s_acc, xt, t_acc,vertical);
        
        //copy image patch in lower right patch
        ys = s_ul;
        yt = t_ul + win_shape/2;
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {   
        //init upperleft rect with lower right image part
        ys = s_ul+ img_shape - win_shape/2;
        yt = t_ul;
        copyImage(ys, s_lr, s_acc, yt, t_acc);
        
        //init upperright rect with images lower left part
        ys = s_ul + Diff2D(0, img_shape.y-win_shape.y/2);
        yt = t_ul + Diff2D(win_shape.x/2,0);
        copyImage(ys, ys+Diff2D(win_shape.x, win_shape.y/2), s_acc, yt, t_acc);
        
        //init lowerleft rect with with images upper right part
        xs = s_ul + Diff2D(img_shape.x-win_shape.x/2, 0);
        xt = t_ul + Diff2D(0,win_shape.y/2);
        copyImage(xs, xs+Diff2D(win_shape.x/2, win_shape.y), s_acc, xt, t_acc);
        
        //copy image patch in lower right patch
        ys = s_ul;
        yt = t_ul + win_shape/2;
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_ZEROPAD)
    {
        initImage(t_ul, t_lr, t_acc, 0);
        
        //copy image patch in lower right patch
        ys = s_ul;
        yt = t_ul + win_shape/2;
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    } 
    
    
    yt = t_ul;
    yt.y += win_shape.y/2;
    yd = d_ul;
    
    for( ; yt.y != t_lr.y-win_shape.y/2; ++yd.y,  ++yt.y)
    {
        xt = yt;
        xt.x += win_shape.x/2;
        
        xd = yd;
        
        for( ; xt.x != t_lr.x-win_shape.x/2; xd.x++, xt.x++)
        {
            func(xt, t_acc, xd, d_acc);             
        }
    }
    
    
    /**********************************************************************************
     *                                                                                * 
     *                             CORNER CASE 2 (6/8):                               * 
     *                                                                                *
     *            HANDLE UPPERRIGHT PIXELS WHERE THE MASK IS MISSING MAXIMUM          *
     *                 x-COORDINATES AND MINIMUM y-COORDINATES                        *
     *                                                                                *
     **********************************************************************************/
    if(border == BORDER_TREATMENT_REPEAT)
    {
        //init upperright rect with single value
        ys = s_ul + Diff2D(img_shape.x-1,0);
        yt = t_ul + Diff2D(win_shape.x,0);
        initImage(yt, yt+win_shape/2, t_acc, s_acc(ys));
        
        //init upperleft rect with vertical stripes
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,0);;
        yt = t_ul;      
        Diff2D lineDiff(win_shape.x,1);
        for( ; yt.y != t_lr.y-win_shape.y ; ++yt.y)
        {
            copyImage(ys, ys+lineDiff, s_acc, yt, t_acc);
        }
        
        //init lowerright rect with horizontal stripes
        xs = s_ul + Diff2D(img_shape.x-1,0);;
        xt = t_ul + Diff2D(win_shape.x,win_shape.y/2);
        Diff2D rowDiff(1, win_shape.y);
        for( ; xt.x != t_lr.x; ++xt.x)
        {
            copyImage(xs, xs+rowDiff, s_acc, xt, t_acc);
        }
        
        //copy image patch in lower left patch
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,0);
        yt = t_ul + Diff2D(0, win_shape.y/2);
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        //init upperright rect with double flipped image
        ys = s_ul + Diff2D(img_shape.x-win_shape.x/2,0);
        yt = t_ul + Diff2D(win_shape.x,0);      
        rotateImage(ys, ys+win_shape/2, s_acc, yt, t_acc, 180);
        
        //init upperleft rect with horizontal reflected image
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,0);
        yt = t_ul;
        reflectImage(ys, ys+Diff2D(win_shape.x, win_shape.y/2), s_acc, yt, t_acc, horizontal);
        
        //init lowerright rect with  vertical reflected image
        xs = s_ul + Diff2D(img_shape.x-win_shape.x/2,0);
        xt = t_ul + Diff2D(win_shape.x,win_shape.y/2);      
        reflectImage(xs, xs+Diff2D(win_shape.x/2, win_shape.y), s_acc, xt, t_acc,vertical);
        
        //copy image patch in lower left patch
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,0);
        yt = t_ul + Diff2D(0, win_shape.y/2);
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
        //init upperright rect with lower left image part
        ys = s_ul + Diff2D(0, img_shape.y-win_shape.y/2);
        yt = t_ul + Diff2D(win_shape.x,0);      
        copyImage(ys, ys+win_shape/2, s_acc, yt, t_acc);
        
        //init upperleft rect with lower right image part
        ys = s_ul + Diff2D(img_shape.x-win_shape.x, img_shape.y-win_shape.y/2);
        yt = t_ul;
        copyImage(ys, ys+Diff2D(win_shape.x, win_shape.y/2), s_acc, yt, t_acc);
        
        //init lowerright rect with upperleft image part
        xs = s_ul;
        xt = t_ul + Diff2D(win_shape.x,win_shape.y/2);      
        copyImage(xs, xs+Diff2D(win_shape.x/2, win_shape.y), s_acc, xt, t_acc);
        
        //copy image patch in lower left patch
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,0);
        yt = t_ul + Diff2D(0, win_shape.y/2);
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_ZEROPAD)
    {
        initImage(t_ul, t_lr, t_acc, 0);
        
        //copy image patch in lower left patch
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,0);
        yt = t_ul + Diff2D(0, win_shape.y/2);
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    } 

    
    yt = t_ul;
    yt.y += win_shape.y/2;
    yd = d_ul + Diff2D(img_shape.x-win_shape.x/2-1, 0);     
    
    for( ; yt.y != t_lr.y-win_shape.y/2; ++yd.y,  ++yt.y)
    {
        xt = yt;
        xt.x += win_shape.x/2;
        
        xd = yd;
        
        for( ; xt.x != t_lr.x-win_shape.x/2; xd.x++, xt.x++)
        {
            func(xt, t_acc, xd, d_acc);                         
        }
    }
    
    /**********************************************************************************
     *                                                                                * 
     *                             CORNER CASE 3 (7/8):                               * 
     *                                                                                *
     *            HANDLE LOWERLEFT PIXELS WHERE THE MASK IS MISSING MINIMUM           *
     *                 x-COORDINATES AND MAXIMUM y-COORDINATES                        *
     *                                                                                *
     **********************************************************************************/
    if(border == BORDER_TREATMENT_REPEAT)
    {
        //init lowerleft rect with single value
        ys = s_ul + Diff2D(0,img_shape.y-1);
        yt = t_ul + Diff2D(0,win_shape.y);
        initImage(yt, yt+win_shape/2, t_acc, s_acc(ys));
        
        //init lowerright rect with vertical stripes
        ys = s_ul + Diff2D(0,img_shape.y-1);
        yt = t_ul + Diff2D(win_shape.x/2,win_shape.y);
        Diff2D lineDiff(win_shape.x,1);
        for( ; yt.y != t_lr.y ; ++yt.y)
        {
            copyImage(ys, ys+lineDiff, s_acc, yt, t_acc);
        }
        
        //init upperleft rect with horizontal stripes
        xs = s_ul + Diff2D(0,img_shape.y-win_shape.y);
        xt = t_ul;
        Diff2D rowDiff(1, win_shape.y);
        for( ; xt.x != t_lr.x-win_shape.x; ++xt.x)
        {
            copyImage(xs, xs+rowDiff, s_acc, xt, t_acc);
        }
        
        //copy image patch in upper right patch
        yt = t_ul + Diff2D(win_shape.x/2,0);
        ys = s_ul + Diff2D(0,img_shape.y-win_shape.y);
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        //init lowerleft rect with double reflected image
        ys = s_ul + Diff2D(0,img_shape.y-win_shape.y/2);
        yt = t_ul + Diff2D(0,win_shape.y);
        rotateImage(ys, ys+win_shape/2, s_acc, yt, t_acc, 180);
        
        //init lowerright rect with horizontal reflected image
        ys = s_ul + Diff2D(0,img_shape.y-win_shape.y/2);
        yt = t_ul + Diff2D(win_shape.x/2,win_shape.y);
        reflectImage(ys, ys+Diff2D(win_shape.x, win_shape.y/2), s_acc, yt, t_acc, horizontal);
        
        //init upperleft rect with vertical reflected image
        xs = s_ul + Diff2D(0,img_shape.y-win_shape.y);
        xt = t_ul;
        reflectImage(xs, xs+Diff2D(win_shape.x/2, win_shape.y), s_acc, xt, t_acc,vertical);
        
        //copy image patch in upper right patch
        yt = t_ul + Diff2D(win_shape.x/2,0);
        ys = s_ul + Diff2D(0,img_shape.y-win_shape.y);
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
        //init lowerleft rect with upper right image part
        ys = s_ul + Diff2D(img_shape.x-win_shape.x/2,0);
        yt = t_ul + Diff2D(0,win_shape.y);
        copyImage(ys, ys+win_shape/2, s_acc, yt, t_acc);
        
        //init lowerright rect with upper left image part
        ys = s_ul;
        yt = t_ul + Diff2D(win_shape.x/2,win_shape.y);
        copyImage(ys, ys+Diff2D(win_shape.x, win_shape.y/2), s_acc, yt, t_acc);
        
        //init upperleft rect with lower right image part
        xs = s_ul + Diff2D(img_shape.x-win_shape.x/2,img_shape.y-win_shape.y);
        xt = t_ul;
        copyImage(xs, xs+Diff2D(win_shape.x/2, win_shape.y), s_acc, xt, t_acc);
        
        //copy image patch in upper right patch
        yt = t_ul + Diff2D(win_shape.x/2,0);
        ys = s_ul + Diff2D(0,img_shape.y-win_shape.y);
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_ZEROPAD)
    {
        initImage(t_ul, t_lr, t_acc, 0);
        
        //copy image patch in upper right patch
        yt = t_ul + Diff2D(win_shape.x/2,0);
        ys = s_ul + Diff2D(0,img_shape.y-win_shape.y);
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    } 
    
    
    yt = t_ul;
    yt.y += win_shape.y/2;
    yd = d_ul + Diff2D(0, img_shape.y-win_shape.y/2-1);
    
    for( ; yt.y != t_lr.y-win_shape.y/2; ++yd.y,  ++yt.y)
    {
        xt = yt;
        xt.x += win_shape.x/2;
        
        xd = yd;
        
        for( ; xt.x != t_lr.x-win_shape.x/2; xd.x++, xt.x++)
        {
            func(xt, t_acc, xd, d_acc);                             
        }
    }
    
    /**********************************************************************************
     *                                                                                * 
     *                             CORNER CASE 4 (8/8):                               * 
     *                                                                                *
     *            HANDLE LOWERRIGHT PIXELS WHERE THE MASK IS MISSING MAXIMUM          *
     *                 x-COORDINATES AND MAXIMUM y-COORDINATES                        *
     *                                                                                *
     **********************************************************************************/
    if(border == BORDER_TREATMENT_REPEAT)
    {
        //init lowerright rect with single value
        ys = s_ul + Diff2D(img_shape.x-1,img_shape.y-1);
        yt = t_ul + win_shape;
        initImage(yt,yt+win_shape/2, t_acc, s_acc(ys));
        
        //init lowerleft rect with vertical stripes
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,img_shape.y-1);
        yt = t_ul + Diff2D(0,win_shape.y);
        Diff2D lineDiff(win_shape.x,1);
        for( ; yt.y != t_lr.y ; ++yt.y)
        {
            copyImage(ys, ys+lineDiff, s_acc, yt, t_acc);
        }
        
        //init upperright rect with horizontal stripes
        xs = s_ul + Diff2D(img_shape.x-1,img_shape.y-win_shape.y);
        xt = t_ul + Diff2D(win_shape.x,0);
        Diff2D rowDiff(1, win_shape.y);
        for( ; xt.x != t_lr.x; ++xt.x)
        {
            copyImage(xs, xs+rowDiff, s_acc, xt, t_acc);
        }
        
        //copy image patch in upperleft patch
        ys = s_ul  + img_shape - win_shape;
        yt = t_ul;
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_REFLECT)
    {
        //init lowerright rect  with double reflected image
        ys = s_ul + img_shape-win_shape/2;
        yt = t_ul + win_shape;
        rotateImage(ys, ys+win_shape/2, s_acc, yt, t_acc, 180);
        
        //init lowerleft rect with horizontal reflected image
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,img_shape.y-win_shape.y/2);
        yt = t_ul + Diff2D(0,win_shape.y);
        reflectImage(ys, ys+Diff2D(win_shape.x, win_shape.y/2), s_acc, yt, t_acc, horizontal);
        
        //init upperright rect with vertical reflected image
        xs = s_ul + Diff2D(img_shape.x-win_shape.x/2,img_shape.y-win_shape.y);
        xt = t_ul + Diff2D(win_shape.x,0);
        reflectImage(xs, xs+Diff2D(win_shape.x/2, win_shape.y), s_acc, xt, t_acc,vertical);
        
        //copy image patch in upperleft patch
        ys = s_ul  + img_shape - win_shape;
        yt = t_ul;
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_WRAP)
    {
        //init lowerright with upperleft image part
        ys = s_ul;
        yt = t_ul + win_shape;
        copyImage(ys, ys+win_shape/2, s_acc, yt, t_acc);
        
        //init lowerleft rect with upper right image part
        ys = s_ul + Diff2D(img_shape.x-win_shape.x,0);
        yt = t_ul + Diff2D(0,win_shape.y);
        copyImage(ys, ys+Diff2D(win_shape.x, win_shape.y/2), s_acc, yt, t_acc);
        
        //init upperright rect with lower left image part
        xs = s_ul + Diff2D(0,img_shape.y-win_shape.y);
        xt = t_ul + Diff2D(win_shape.x,0);
        copyImage(xs, xs+Diff2D(win_shape.x/2, win_shape.y), s_acc, xt, t_acc);
        
        //copy image patch in upperleft patch
        ys = s_ul  + img_shape-win_shape;
        yt = t_ul;
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    }
    else if(border == BORDER_TREATMENT_ZEROPAD)
    {
        initImage(t_ul, t_lr, t_acc, 0);        
        
        //copy image patch in upperleft patch
        ys = s_ul  + img_shape - win_shape;
        yt = t_ul;
        copyImage(ys, ys+win_shape, s_acc, yt, t_acc);
        
    } 
    
    
    yt = t_ul;
    yt.y += win_shape.y/2;
    yd = d_ul + img_shape-win_shape/2-Diff2D(1,1);
    
    for( ; yt.y != t_lr.y-win_shape.y/2; ++yd.y,  ++yt.y)
    {
        xt = yt;
        xt.x += win_shape.x/2;
        
        xd = yd;
        
        for( ; xt.x != t_lr.x-win_shape.x/2; xd.x++, xt.x++)
        {
            func(xt, t_acc, xd, d_acc);                             
        }
    }
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class ProcessingFunctor>
inline void applyWindowFunction(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                                pair<DestIterator, DestAccessor> d, 
                                ProcessingFunctor func, 
                                BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    applyWindowFunction(s.first, s.second, s.third,
                         d.first, d.second, 
                         func, 
                         border);
}

template <class T1, class S1, 
          class T2, class S2,
          class ProcessingFunctor>
inline void applyWindowFunction(MultiArrayView<2, T1, S1> const & src,
                                MultiArrayView<2, T2, S2> dest, 
                                ProcessingFunctor func, 
                                BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    vigra_precondition(src.shape() == dest.shape(),
                        "vigra::applyWindowFunction(): shape mismatch between input and output.");
    applyWindowFunction(srcImageRange(src),
                            destImage(dest), 
                            func, 
                            border);
}

//@}

} //end of namespace vigra

#endif //VIGRA_APPLYWINDOWFUNCTION_HXX
