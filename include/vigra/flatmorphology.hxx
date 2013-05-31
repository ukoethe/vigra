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
 
 
#ifndef VIGRA_FLATMORPHOLOGY_HXX
#define VIGRA_FLATMORPHOLOGY_HXX

#include <cmath>
#include <vector>
#include "utilities.hxx"
#include "multi_shape.hxx"

namespace vigra {

/** \addtogroup Morphology Basic Morphological Operations
    Perform erosion, dilation, and median with disc structuring functions
    
    See also: \ref MultiArrayMorphology Separable morphology with parabola structuring functions in arbitrary dimensions
*/
//@{

/********************************************************/
/*                                                      */
/*                  discRankOrderFilter                 */
/*                                                      */
/********************************************************/

/** \brief Apply rank order filter with disc structuring function to the image.

    The pixel values of the source image <b> must</b> be in the range
    0...255. Radius must be >= 0. Rank must be in the range 0.0 <= rank 
    <= 1.0. The filter acts as a minimum filter if rank = 0.0, 
    as a median if rank = 0.5, and as a maximum filter if rank = 1.0.
    Accessor are used to access the pixel data.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        discRankOrderFilter(SrcIterator upperleft1, 
                            SrcIterator lowerright1, SrcAccessor sa,
                            DestIterator upperleft2, DestAccessor da,
                            int radius, float rank)
    }
    \endcode
    
    pass \ref ImageIterators and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        discRankOrderFilter(SrcIterator upperleft1, 
                            SrcIterator lowerright1, SrcAccessor sa,
                            DestIterator upperleft2, DestAccessor da,
                            int radius, float rank)
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        discRankOrderFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            int radius, float rank)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/flatmorphology.hxx\><br>
    Namespace: vigra
    
    \code
    vigra::CImage src, dest;
    
    // do median filtering
    vigra::discRankOrderFilter(srcImageRange(src), destImage(dest), 10, 0.5);
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcIterator src_upperleft;
    DestIterator dest_upperleft;
    int x, y;
    unsigned char value;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    // value_type of accessor must be convertible to unsigned char
    value = src_accessor(src_upperleft, x, y); 
    
    dest_accessor.set(value, dest_upperleft, x, y);
    \endcode
    
    <b> Preconditions:</b>
    
    \code
    for all source pixels: 0 <= value <= 255
    
    (rank >= 0.0) && (rank <= 1.0)
    radius >= 0
    \endcode
    
*/
doxygen_overloaded_function(template <...> void discRankOrderFilter)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
discRankOrderFilter(SrcIterator upperleft1, 
                    SrcIterator lowerright1, SrcAccessor sa,
                    DestIterator upperleft2, DestAccessor da,
                    int radius, float rank)
{
    vigra_precondition((rank >= 0.0) && (rank <= 1.0),
            "discRankOrderFilter(): Rank must be between 0 and 1"
            " (inclusive).");
    
    vigra_precondition(radius >= 0,
            "discRankOrderFilter(): Radius must be >= 0.");
    
    int i, x, y, xmax, ymax, xx, yy;
    int rankpos, winsize, leftsum;
    
    long hist[256];
    
    // prepare structuring function
    std::vector<int> struct_function(radius+1);
    struct_function[0] = radius;
    
    double r2 = (double)radius*radius;
    for(i=1; i<=radius; ++i)
    {
        double r = (double) i - 0.5;
        struct_function[i] = (int)(VIGRA_CSTD::sqrt(r2 - r*r) + 0.5);
    }

    int w = lowerright1.x - upperleft1.x;
    int h = lowerright1.y - upperleft1.y;
    
    SrcIterator ys(upperleft1);
    DestIterator yd(upperleft2);

    for(y=0; y<h; ++y, ++ys.y, ++yd.y)
    {
        SrcIterator xs(ys);
        DestIterator xd(yd);
        
        // first column
        int x0 = 0;
        int y0 = y;
        int x1 = w - 1;
        int y1 = h - y - 1;

        // clear histogram
        for(i=0; i<256; ++i) hist[i] = 0;
        winsize = 0;
        
        // init histogram
        ymax = (y1 < radius) ? y1 : radius;
        for(yy=0; yy<=ymax; ++yy)
        {
            xmax = (x1 < struct_function[yy]) ? x1 : struct_function[yy];
            for(xx=0; xx<=xmax; ++xx)
            {
                hist[sa(xs, Diff2D(xx, yy))]++;
                winsize++;
            }
        }
        
        ymax = (y0 < radius) ? y0 : radius;
        for(yy=1; yy<=ymax; ++yy)
        {
            xmax = (x1 < struct_function[yy]) ? x1 : struct_function[yy];
            for(xx=0; xx<=xmax; ++xx)
            {
                hist[sa(xs, Diff2D(xx, -yy))]++;
                winsize++;
            }
        }
    
        // find the desired histogram bin 
        leftsum = 0;
        if(rank == 0.0)
        {
            for(i=0; i<256; i++)
            {
                if(hist[i]) break;
            }
            rankpos = i;
        }
        else
        {
            for(i=0; i<256; i++)
            {
                if((float)(hist[i]+leftsum) / winsize >= rank) break;
                leftsum += hist[i];
            }
            rankpos = i;
        }
        
        da.set(rankpos, xd);
        
        ++xs.x;
        ++xd.x;
        
        // inner columns
        for(x=1; x<w; ++x, ++xs.x, ++xd.x)
        {
            x0 = x;
            y0 = y;
            x1 = w - x - 1;
            y1 = h - y - 1;
            
            // update histogram 
            // remove pixels at left border 
            yy = (y1 < radius) ? y1 : radius;
            for(; yy>=0; yy--)
            {
                unsigned char cur;
                xx = struct_function[yy]+1;
                if(xx > x0) break;
                
                cur = sa(xs, Diff2D(-xx, yy));
                
                hist[cur]--;
                if(cur < rankpos) leftsum--;
                winsize--;
            }
            yy = (y0 < radius) ? y0 : radius;
            for(; yy>=1; yy--)
            {
                unsigned char cur;
                xx = struct_function[yy]+1;
                if(xx > x0) break;
                
                cur = sa(xs, Diff2D(-xx, -yy));
                
                hist[cur]--;
                if(cur < rankpos) leftsum--;
                winsize--;
            }
            
            // add pixels at right border 
            yy = (y1 < radius) ? y1 : radius;
            for(; yy>=0; yy--)
            {
                unsigned char cur;
                xx = struct_function[yy];
                if(xx > x1) break;
                
                cur = sa(xs, Diff2D(xx, yy));
                
                hist[cur]++;
                if(cur < rankpos) leftsum++;
                winsize++;
            }
            yy = (y0 < radius) ? y0 : radius;
            for(; yy>=1; yy--)
            {
                unsigned char cur;
                xx = struct_function[yy];
                if(xx > x1) break;
                
                cur = sa(xs, Diff2D(xx, -yy));
                
                hist[cur]++;
                if(cur < rankpos) leftsum++;
                winsize++;
            }
        
            // find the desired histogram bin 
            if(rank == 0.0)
            {
                if(leftsum == 0)
                {
                    // search to the right 
                    for(i=rankpos; i<256; i++)
                    {
                        if(hist[i]) break;
                    }
                    rankpos = i;
                }
                else
                {
                    // search to the left 
                    for(i=rankpos-1; i>=0; i--)
                    {
                        leftsum -= hist[i];
                        if(leftsum == 0) break;
                    }
                    rankpos = i;
                }
            }
            else  // rank > 0.0 
            {
                if((float)leftsum / winsize < rank)
                {
                    // search to the right 
                    for(i=rankpos; i<256; i++)
                    {
                        if((float)(hist[i]+leftsum) / winsize >= rank) break;
                        leftsum+=hist[i];
                    }
                    rankpos = i;
                }
                else
                {
                    // search to the left 
                    for(i=rankpos-1; i>=0; i--)
                    {
                        leftsum-=hist[i];
                        if((float)leftsum / winsize < rank) break;
                    }
                    rankpos = i;
                }
            }
                    
            da.set(rankpos, xd);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
discRankOrderFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    int radius, float rank)
{
    discRankOrderFilter(src.first, src.second, src.third,
                        dest.first, dest.second,
                        radius, rank);
}

template <class T1, class S1,
          class T2, class S2>
void
discRankOrderFilter(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> dest,
                    int radius, float rank)
{
    discRankOrderFilter(srcImageRange(src),
                        destImage(dest),
                        radius, rank);
}

/********************************************************/
/*                                                      */
/*                      discErosion                     */
/*                                                      */
/********************************************************/

/** \brief Apply erosion (minimum) filter with disc of given radius to image.

    This is an abbreviation vor the rank order filter with rank = 0.0.
    See \ref discRankOrderFilter() for more information.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discErosion(SrcIterator upperleft1, 
                    SrcIterator lowerright1, SrcAccessor sa,
                    DestIterator upperleft2, DestAccessor da,
                    int radius)
    }
    \endcode
    
    pass \ref ImageIterators and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discErosion(SrcIterator upperleft1, 
                    SrcIterator lowerright1, SrcAccessor sa,
                    DestIterator upperleft2, DestAccessor da,
                    int radius)
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        discErosion(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    int radius)
    }
    \endcode

*/
doxygen_overloaded_function(template <...> void discErosion)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
discErosion(SrcIterator upperleft1, 
            SrcIterator lowerright1, SrcAccessor sa,
            DestIterator upperleft2, DestAccessor da,
            int radius)
{
    vigra_precondition(radius >= 0, "discErosion(): Radius must be >= 0.");
    
    discRankOrderFilter(upperleft1, lowerright1, sa, 
                        upperleft2, da, radius, 0.0);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
discErosion(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
            int radius)
{
    vigra_precondition(radius >= 0, "discErosion(): Radius must be >= 0.");
    
    discRankOrderFilter(src.first, src.second, src.third,
                        dest.first, dest.second,
                        radius, 0.0);
}

template <class T1, class S1,
          class T2, class S2>
inline void
discErosion(MultiArrayView<2, T1, S1> const & src,
            MultiArrayView<2, T2, S2> dest,
            int radius)
{
    discErosion(srcImageRange(src), destImage(dest), radius);
}

/********************************************************/
/*                                                      */
/*                     discDilation                     */
/*                                                      */
/********************************************************/

/** \brief Apply dilation (maximum) filter with disc of given radius to image.

    This is an abbreviation vor the rank order filter with rank = 1.0.
    See \ref discRankOrderFilter() for more information.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discDilation(SrcIterator upperleft1, 
                    SrcIterator lowerright1, SrcAccessor sa,
                    DestIterator upperleft2, DestAccessor da,
                    int radius)
    }
    \endcode
    
    pass \ref ImageIterators and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discDilation(SrcIterator upperleft1, 
                    SrcIterator lowerright1, SrcAccessor sa,
                    DestIterator upperleft2, DestAccessor da,
                    int radius)
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        discDilation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    int radius)
    }
    \endcode

*/
doxygen_overloaded_function(template <...> void discDilation)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
discDilation(SrcIterator upperleft1, 
            SrcIterator lowerright1, SrcAccessor sa,
            DestIterator upperleft2, DestAccessor da,
            int radius)
{
    vigra_precondition(radius >= 0, "discDilation(): Radius must be >= 0.");
    
    discRankOrderFilter(upperleft1, lowerright1, sa, 
                        upperleft2, da, radius, 1.0);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
discDilation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
             pair<DestIterator, DestAccessor> dest,
             int radius)
{
    vigra_precondition(radius >= 0, "discDilation(): Radius must be >= 0.");
    
    discRankOrderFilter(src.first, src.second, src.third,
                        dest.first, dest.second,
                        radius, 1.0);
}

template <class T1, class S1,
          class T2, class S2>
inline void
discDilation(MultiArrayView<2, T1, S1> const & src,
             MultiArrayView<2, T2, S2> dest,
             int radius)
{
    discDilation(srcImageRange(src), destImage(dest), radius);
}

/********************************************************/
/*                                                      */
/*                      discMedian                      */
/*                                                      */
/********************************************************/

/** \brief Apply median filter with disc of given radius to image.

    This is an abbreviation vor the rank order filter with rank = 0.5.
    See \ref discRankOrderFilter() for more information.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discMedian(SrcIterator upperleft1, 
                    SrcIterator lowerright1, SrcAccessor sa,
                    DestIterator upperleft2, DestAccessor da,
                    int radius)
    }
    \endcode
    
    pass \ref ImageIterators and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discMedian(SrcIterator upperleft1, 
                    SrcIterator lowerright1, SrcAccessor sa,
                    DestIterator upperleft2, DestAccessor da,
                    int radius)
    }
    \endcode
    
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        discMedian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    int radius)
    }
    \endcode

*/
doxygen_overloaded_function(template <...> void discMedian)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
discMedian(SrcIterator upperleft1, 
            SrcIterator lowerright1, SrcAccessor sa,
            DestIterator upperleft2, DestAccessor da,
            int radius)
{
    vigra_precondition(radius >= 0, "discMedian(): Radius must be >= 0.");
    
    discRankOrderFilter(upperleft1, lowerright1, sa, 
                        upperleft2, da, radius, 0.5);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
discMedian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
           pair<DestIterator, DestAccessor> dest,
           int radius)
{
    vigra_precondition(radius >= 0, "discMedian(): Radius must be >= 0.");
    
    discRankOrderFilter(src.first, src.second, src.third,
                        dest.first, dest.second,
                        radius, 0.5);
}

template <class T1, class S1,
          class T2, class S2>
inline void
discMedian(MultiArrayView<2, T1, S1> const & src,
           MultiArrayView<2, T2, S2> dest,
           int radius)
{
    discMedian(srcImageRange(src), destImage(dest), radius);
}

/********************************************************/
/*                                                      */
/*            discRankOrderFilterWithMask               */
/*                                                      */
/********************************************************/

/** \brief Apply rank order filter with disc structuring function to the image
    using a mask.
    
    The pixel values of the source image <b> must</b> be in the range
    0...255. Radius must be >= 0. Rank must be in the range 0.0 <= rank 
    <= 1.0. The filter acts as a minimum filter if rank = 0.0, 
    as a median if rank = 0.5, and as a maximum filter if rank = 1.0.
    Accessor are used to access the pixel data.
    
    The mask is only applied to th input image, i.e. the function
    generates an output wherever the current disc contains at least 
    one pixel with mask value 'true'. Source pixels with mask value
    'false' are ignored during the calculation of the rank order.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void
        discRankOrderFilterWithMask(SrcIterator upperleft1, 
                                    SrcIterator lowerright1, SrcAccessor sa,
                                    MaskIterator upperleftm, MaskAccessor mask,
                                    DestIterator upperleft2, DestAccessor da,
                                    int radius, float rank)
    }
    \endcode
    
    pass \ref ImageIterators and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void
        discRankOrderFilterWithMask(SrcIterator upperleft1, 
                                    SrcIterator lowerright1, SrcAccessor sa,
                                    MaskIterator upperleftm, MaskAccessor mask,
                                    DestIterator upperleft2, DestAccessor da,
                                    int radius, float rank)
    }
    \endcode
    
    
    group arguments (use in conjunction with \ref ArgumentObjectFactories):
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void
        discRankOrderFilterWithMask(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                    pair<MaskIterator, MaskAccessor> mask,
                                    pair<DestIterator, DestAccessor> dest,
                                    int radius, float rank)
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/flatmorphology.hxx\><br>
    Namespace: vigra
    
    \code
    vigra::CImage src, dest, mask;
    
    // do median filtering
    vigra::discRankOrderFilterWithMask(srcImageRange(src), 
                                maskImage(mask), destImage(dest), 10, 0.5);
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcIterator src_upperleft;
    DestIterator dest_upperleft;
    MaskIterator mask_upperleft;
    int x, y;
    unsigned char value;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    MaskAccessor mask_accessor;
                     
    mask_accessor(mask_upperleft, x, y) // convertible to bool
    
    // value_type of accessor must be convertible to unsigned char
    value = src_accessor(src_upperleft, x, y); 
    
    dest_accessor.set(value, dest_upperleft, x, y);
    \endcode
    
    <b> Preconditions:</b>
    
    \code
    for all source pixels: 0 <= value <= 255
    
    (rank >= 0.0) && (rank <= 1.0)
    radius >= 0
    \endcode
    
*/
doxygen_overloaded_function(template <...> void discRankOrderFilterWithMask)

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
void
discRankOrderFilterWithMask(SrcIterator upperleft1, 
                            SrcIterator lowerright1, SrcAccessor sa,
                            MaskIterator upperleftm, MaskAccessor mask,
                            DestIterator upperleft2, DestAccessor da,
                            int radius, float rank)
{
    vigra_precondition((rank >= 0.0) && (rank <= 1.0),
                 "discRankOrderFilter(): Rank must be between 0 and 1"
                 " (inclusive).");
    
    vigra_precondition(radius >= 0, "discRankOrderFilter(): Radius must be >= 0.");
    
    int i, x, y, xmax, ymax, xx, yy;
    int rankpos, winsize, leftsum;
    
    long hist[256];
    
    // prepare structuring function
    std::vector<int> struct_function(radius+1);
    struct_function[0] = radius;
    
    double r2 = (double)radius*radius;
    for(i=1; i<=radius; ++i)
    {
        double r = (double) i - 0.5;
        struct_function[i] = (int)(VIGRA_CSTD::sqrt(r2 - r*r) + 0.5);
    }

    int w = lowerright1.x - upperleft1.x;
    int h = lowerright1.y - upperleft1.y;
        
    SrcIterator ys(upperleft1);
    MaskIterator ym(upperleftm);
    DestIterator yd(upperleft2);

    for(y=0; y<h; ++y, ++ys.y, ++yd.y, ++ym.y)
    {
        SrcIterator xs(ys);
        MaskIterator xm(ym);
        DestIterator xd(yd);
        
        // first column
        int x0 = 0;
        int y0 = y;
        int x1 = w - 1;
        int y1 = h - y - 1;

        // clear histogram
        for(i=0; i<256; ++i) hist[i] = 0;
        winsize = 0;
        leftsum = 0;
        rankpos = 0;
        
        // init histogram
        ymax = (y1 < radius) ? y1 : radius;
        for(yy=0; yy<=ymax; ++yy)
        {
            xmax = (x1 < struct_function[yy]) ? x1 : struct_function[yy];
            for(xx=0; xx<=xmax; ++xx)
            {
                Diff2D pos(xx, yy);
                if(mask(xm, pos))
                {
                    hist[sa(xs, pos)]++;
                    winsize++;
                }
            }
        }
        
        ymax = (y0 < radius) ? y0 : radius;
        for(yy=1; yy<=ymax; ++yy)
        {
            xmax = (x1 < struct_function[yy]) ? x1 : struct_function[yy];
            for(xx=0; xx<=xmax; ++xx)
            {
                Diff2D pos(xx, -yy);
                if(mask(xm, pos))
                {
                    hist[sa(xs, pos)]++;
                    winsize++;
                }
            }
        }
    
        // find the desired histogram bin 
        if(winsize) 
        {
            if(rank == 0.0)
            {
                for(i=0; i<256; i++)
                {
                    if(hist[i]) break;
                }
                rankpos = i;
            }
            else
            {
                for(i=0; i<256; i++)
                {
                    if((float)(hist[i]+leftsum) / winsize >= rank) break;
                    leftsum += hist[i];
                }
                rankpos = i;
            }
            
            da.set(rankpos, xd);
        }
            
        ++xs.x;
        ++xd.x;
        ++xm.x;
        
        // inner columns
        for(x=1; x<w; ++x, ++xs.x, ++xd.x, ++xm.x)
        {
            x0 = x;
            y0 = y;
            x1 = w - x - 1;
            y1 = h - y - 1;
            
            // update histogram 
            // remove pixels at left border 
            yy = (y1 < radius) ? y1 : radius;
            for(; yy>=0; yy--)
            {
                unsigned char cur;
                xx = struct_function[yy]+1;
                if(xx > x0) break;
                
                Diff2D pos(-xx, yy);
                if(mask(xm, pos))
                {
                    cur = sa(xs, pos);
                    
                    hist[cur]--;
                    if(cur < rankpos) leftsum--;
                    winsize--;
                }
            }
            yy = (y0 < radius) ? y0 : radius;
            for(; yy>=1; yy--)
            {
                unsigned char cur;
                xx = struct_function[yy]+1;
                if(xx > x0) break;
                
                Diff2D pos(-xx, -yy);
                if(mask(xm, pos))
                {
                    cur = sa(xs, pos);
                    
                    hist[cur]--;
                    if(cur < rankpos) leftsum--;
                    winsize--;
                }
            }
            
            // add pixels at right border 
            yy = (y1 < radius) ? y1 : radius;
            for(; yy>=0; yy--)
            {
                unsigned char cur;
                xx = struct_function[yy];
                if(xx > x1) break;
                
                Diff2D pos(xx, yy);
                if(mask(xm, pos))
                {
                    cur = sa(xs, pos);
                    
                    hist[cur]++;
                    if(cur < rankpos) leftsum++;
                    winsize++;
                }
            }
            yy = (y0 < radius) ? y0 : radius;
            for(; yy>=1; yy--)
            {
                unsigned char cur;
                xx = struct_function[yy];
                if(xx > x1) break;
                
                Diff2D pos(xx, -yy);
                if(mask(xm, pos))
                {
                    cur = sa(xs, pos);
                    
                    hist[cur]++;
                    if(cur < rankpos) leftsum++;
                    winsize++;
                }
            }
        
            // find the desired histogram bin 
            if(winsize) 
            {
                if(rank == 0.0)
                {
                    if(leftsum == 0)
                    {
                        // search to the right 
                        for(i=rankpos; i<256; i++)
                        {
                            if(hist[i]) break;
                        }
                        rankpos = i;
                    }
                    else
                    {
                        // search to the left 
                        for(i=rankpos-1; i>=0; i--)
                        {
                            leftsum -= hist[i];
                            if(leftsum == 0) break;
                        }
                        rankpos = i;
                    }
                }
                else  // rank > 0.0 
                {
                    if((float)leftsum / winsize < rank)
                    {
                        // search to the right 
                        for(i=rankpos; i<256; i++)
                        {
                            if((float)(hist[i]+leftsum) / winsize >= rank) break;
                            leftsum+=hist[i];
                        }
                        rankpos = i;
                    }
                    else
                    {
                        // search to the left 
                        for(i=rankpos-1; i>=0; i--)
                        {
                            leftsum-=hist[i];
                            if((float)leftsum / winsize < rank) break;
                        }
                        rankpos = i;
                    }
                }
                        
                da.set(rankpos, xd);
            }
            else
            {
                leftsum = 0;
                rankpos = 0;
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
void
discRankOrderFilterWithMask(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<MaskIterator, MaskAccessor> mask,
                            pair<DestIterator, DestAccessor> dest,
                            int radius, float rank)
{
    discRankOrderFilterWithMask(src.first, src.second, src.third,
                                mask.first, mask.second,
                                dest.first, dest.second,
                                radius, rank);
}

template <class T1, class S1,
          class TM, class SM,
          class T2, class S2>
void
discRankOrderFilterWithMask(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, TM, SM> const & mask,
                            MultiArrayView<2, T2, S2> dest,
                            int radius, float rank)
{
    discRankOrderFilterWithMask(srcImageRange(src),
                                maskImage(mask),
                                destImage(dest),
                                radius, rank);
}

/********************************************************/
/*                                                      */
/*                 discErosionWithMask                  */
/*                                                      */
/********************************************************/

/** \brief Apply erosion (minimum) filter with disc of given radius to image
    using a mask.
    
    This is an abbreviation vor the masked rank order filter with 
    rank = 0.0. See \ref discRankOrderFilterWithMask() for more information.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discErosionWithMask(SrcIterator upperleft1, 
                            SrcIterator lowerright1, SrcAccessor sa,
                            MaskIterator upperleftm, MaskAccessor mask,
                            DestIterator upperleft2, DestAccessor da,
                            int radius)
    }
    \endcode
    
    pass \ref ImageIterators and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discErosionWithMask(SrcIterator upperleft1, 
                            SrcIterator lowerright1, SrcAccessor sa,
                            MaskIterator upperleftm, MaskAccessor mask,
                            DestIterator upperleft2, DestAccessor da,
                            int radius)
    }
    \endcode
    
    
    group arguments (use in conjunction with \ref ArgumentObjectFactories):
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discErosionWithMask(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<MaskIterator, MaskAccessor> mask,
                            pair<DestIterator, DestAccessor> dest,
                            int radius)
    }
    \endcode

*/
doxygen_overloaded_function(template <...> void discErosionWithMask)

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void 
discErosionWithMask(SrcIterator upperleft1, 
                    SrcIterator lowerright1, SrcAccessor sa,
                    MaskIterator upperleftm, MaskAccessor mask,
                    DestIterator upperleft2, DestAccessor da,
                    int radius)
{
    vigra_precondition(radius >= 0, "discErosionWithMask(): Radius must be >= 0.");
    
    discRankOrderFilterWithMask(upperleft1, lowerright1, sa, 
                                upperleftm, mask,
                                upperleft2, da,
                                radius, 0.0);
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void 
discErosionWithMask(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<MaskIterator, MaskAccessor> mask,
                    pair<DestIterator, DestAccessor> dest,
                    int radius)
{
    vigra_precondition(radius >= 0, "discErosionWithMask(): Radius must be >= 0.");
    
    discRankOrderFilterWithMask(src.first, src.second, src.third,
                                mask.first, mask.second,
                                dest.first, dest.second,
                                radius, 0.0);
}

template <class T1, class S1,
          class TM, class SM,
          class T2, class S2>
inline void 
discErosionWithMask(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, TM, SM> const & mask,
                    MultiArrayView<2, T2, S2> dest,
                    int radius)
{
    discErosionWithMask(srcImageRange(src), maskImage(mask), destImage(dest), radius);
}

/********************************************************/
/*                                                      */
/*                discDilationWithMask                  */
/*                                                      */
/********************************************************/

/** \brief Apply dilation (maximum) filter with disc of given radius to image
    using a mask.
    
    This is an abbreviation vor the masked rank order filter with 
    rank = 1.0. See \ref discRankOrderFilterWithMask() for more information.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discDilationWithMask(SrcIterator upperleft1, 
                            SrcIterator lowerright1, SrcAccessor sa,
                            MaskIterator upperleftm, MaskAccessor mask,
                            DestIterator upperleft2, DestAccessor da,
                            int radius)
    }
    \endcode
    
    pass \ref ImageIterators and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discDilationWithMask(SrcIterator upperleft1, 
                            SrcIterator lowerright1, SrcAccessor sa,
                            MaskIterator upperleftm, MaskAccessor mask,
                            DestIterator upperleft2, DestAccessor da,
                            int radius)
    }
    \endcode
    
    
    group arguments (use in conjunction with \ref ArgumentObjectFactories):
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discDilationWithMask(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<MaskIterator, MaskAccessor> mask,
                            pair<DestIterator, DestAccessor> dest,
                            int radius)
    }
    \endcode

*/
doxygen_overloaded_function(template <...> void discDilationWithMask)

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void 
discDilationWithMask(SrcIterator upperleft1, 
                     SrcIterator lowerright1, SrcAccessor sa,
                     MaskIterator upperleftm, MaskAccessor mask,
                     DestIterator upperleft2, DestAccessor da,
                     int radius)
{
    vigra_precondition(radius >= 0, "discDilationWithMask(): Radius must be >= 0.");
    
    discRankOrderFilterWithMask(upperleft1, lowerright1, sa, 
                                upperleftm, mask,
                                upperleft2, da,
                                radius, 1.0);
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void 
discDilationWithMask(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                     pair<MaskIterator, MaskAccessor> mask,
                     pair<DestIterator, DestAccessor> dest,
                     int radius)
{
    vigra_precondition(radius >= 0, "discDilationWithMask(): Radius must be >= 0.");
    
    discRankOrderFilterWithMask(src.first, src.second, src.third,
                                mask.first, mask.second,
                                dest.first, dest.second,
                                radius, 1.0);
}

template <class T1, class S1,
          class TM, class SM,
          class T2, class S2>
inline void 
discDilationWithMask(MultiArrayView<2, T1, S1> const & src,
                     MultiArrayView<2, TM, SM> const & mask,
                     MultiArrayView<2, T2, S2> dest,
                     int radius)
{
    discDilationWithMask(srcImageRange(src), maskImage(mask), destImage(dest), radius);
}

/********************************************************/
/*                                                      */
/*                 discMedianWithMask                   */
/*                                                      */
/********************************************************/

/** \brief Apply median filter with disc of given radius to image
    using a mask.
    
    This is an abbreviation vor the masked rank order filter with 
    rank = 0.5. See \ref discRankOrderFilterWithMask() for more information.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discMedianWithMask(SrcIterator upperleft1, 
                            SrcIterator lowerright1, SrcAccessor sa,
                            MaskIterator upperleftm, MaskAccessor mask,
                            DestIterator upperleft2, DestAccessor da,
                            int radius)
    }
    \endcode
    
    pass \ref ImageIterators and \ref DataAccessors:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discMedianWithMask(SrcIterator upperleft1, 
                            SrcIterator lowerright1, SrcAccessor sa,
                            MaskIterator upperleftm, MaskAccessor mask,
                            DestIterator upperleft2, DestAccessor da,
                            int radius)
    }
    \endcode
    
    
    group arguments (use in conjunction with \ref ArgumentObjectFactories):
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class MaskIterator, class MaskAccessor,
                  class DestIterator, class DestAccessor>
        void 
        discMedianWithMask(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<MaskIterator, MaskAccessor> mask,
                            pair<DestIterator, DestAccessor> dest,
                            int radius)
    }
    \endcode

*/
doxygen_overloaded_function(template <...> void discMedianWithMask)

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void 
discMedianWithMask(SrcIterator upperleft1, 
                    SrcIterator lowerright1, SrcAccessor sa,
                    MaskIterator upperleftm, MaskAccessor mask,
                    DestIterator upperleft2, DestAccessor da,
                    int radius)
{
    vigra_precondition(radius >= 0, "discMedianWithMask(): Radius must be >= 0.");
    
    discRankOrderFilterWithMask(upperleft1, lowerright1, sa, 
                                upperleftm, mask,
                                upperleft2, da,
                                radius, 0.5);
}

template <class SrcIterator, class SrcAccessor,
          class MaskIterator, class MaskAccessor,
          class DestIterator, class DestAccessor>
inline void 
discMedianWithMask(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<MaskIterator, MaskAccessor> mask,
                   pair<DestIterator, DestAccessor> dest,
                   int radius)
{
    vigra_precondition(radius >= 0, "discMedianWithMask(): Radius must be >= 0.");
    
    discRankOrderFilterWithMask(src.first, src.second, src.third,
                        mask.first, mask.second,
                        dest.first, dest.second,
                        radius, 0.5);
}

template <class T1, class S1,
          class TM, class SM,
          class T2, class S2>
inline void 
discMedianWithMask(MultiArrayView<2, T1, S1> const & src,
                   MultiArrayView<2, TM, SM> const & mask,
                   MultiArrayView<2, T2, S2> dest,
                   int radius)
{
    discMedianWithMask(srcImageRange(src), maskImage(mask), destImage(dest), radius);
}

//@}

} // namespace vigra

#endif // VIGRA_FLATMORPHOLOGY_HXX
