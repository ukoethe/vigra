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
 
 
#ifndef VIGRA_LOCALMINMAX_HXX
#define VIGRA_LOCALMINMAX_HXX

#include <vector>
#include "vigra/utilities.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/labelimage.hxx"

namespace vigra {

/** \addtogroup LocalMinMax Local Minima and Maxima

    Detect local minima and maxima of the gray level,
    including extremal plateaus larger than 1 pixel
*/
//@{

/********************************************************/
/*                                                      */
/*                       localMinima                    */
/*                                                      */
/********************************************************/

/** \brief Find local minima in an image.

    The minima are found only when the have a size of one pixel.
    Use \ref extendedLocalMinima() to find minimal plateaus. Minima are
    marked in the destination image with the given marker value
    (default is 1), all other destination pixels remain unchanged.
    <TT>SrcAccessor::value_type</TT> must be less-comparable.
    A pixel at the image border will never be marked as minimum. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
	          class DestIterator, class DestAccessor, 
	          class DestValue = DestAccessor::value_type>
        void 
        localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
		    DestIterator dul, DestAccessor da, 
		    DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
	          class DestIterator, class DestAccessor, 
	          class DestValue = DestAccessor::value_type>
        void 
        localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
		    pair<DestIterator, DestAccessor> dest,
		    DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BImage src(w,h), minima(w,h);
    
    // init destiniation image
    minima = 0;
    
    vigra::localMinima(srcImageRange(src), destImage(minima));
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);
    
    u < u
    
    DestValue marker;
    dest_accessor.set(marker, dest_upperleft);
    \endcode

*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
void 
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    int w = slr.x - sul.x - 2;
    int h = slr.y - sul.y - 2;

    static const Diff2D dist[] = {
        Diff2D( 1, 0), Diff2D( 1, -1), Diff2D( 0, -1), Diff2D( -1, -1), 
        Diff2D( -1, 0), Diff2D( -1, 1), Diff2D( 0, 1), Diff2D( 1, 1)};
    
    int i,x,y;
    
    sul += Diff2D(1,1);
    dul += Diff2D(1,1);
    
    for(y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
	SrcIterator  sx = sul;
	DestIterator dx = dul;
	
	for(x=0; x<w; ++x, ++sx.x, ++dx.x)
	{
	    for(i=0; i<8; ++i)
	    {
		if(!(sa(sx) < sa(sx, dist[i]))) break;
	    }
	
	    if(i == 8) da.set(marker, dx);
	}
    }
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
localMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    localMinima(sul, slr, sa, dul, da, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
inline void 
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
	    DestValue marker)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
localMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    localMinima(src.first, src.second, src.third,
                dest.first, dest.second,
                NumericTraits<typename DestAccessor::value_type>::one());
}

/********************************************************/
/*                                                      */
/*                       localMaxima                    */
/*                                                      */
/********************************************************/

/** \brief Find local maxima in an image.

    The maxima are found only when the have a size of one pixel.
    Use \ref extendedLocalMaxima() to find maximal plateaus. Maxima are
    marked in the destination image with the given marker value
    (default is 1), all other destination pixels remain unchanged.
    <TT>SrcAccessor::value_type</TT> must be less-comparable.
    A pixel at the image border will never be marked as maximum. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
	          class DestIterator, class DestAccessor, 
	          class DestValue = DestAccessor::value_type>
        void 
        localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
		    DestIterator dul, DestAccessor da, 
		    DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
	          class DestIterator, class DestAccessor, 
	          class DestValue = DestAccessor::value_type>
        void 
        localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
		    pair<DestIterator, DestAccessor> dest,
		    DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BImage src(w,h), maxima(w,h);
    
    // init destiniation image
    maxima = 0;
    
    vigra::localMaxima(srcImageRange(src), destImage(maxima));
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);
    
    u < u
    
    DestValue marker;
    dest_accessor.set(marker, dest_upperleft);
    \endcode

*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
void 
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    int w = slr.x - sul.x - 2;
    int h = slr.y - sul.y - 2;

    static const Diff2D dist[] = {
        Diff2D( 1, 0), Diff2D( 1, -1), Diff2D( 0, -1), Diff2D( -1, -1), 
        Diff2D( -1, 0), Diff2D( -1, 1), Diff2D( 0, 1), Diff2D( 1, 1)};
    
    int i,x,y;
    
    sul += Diff2D(1,1);
    dul += Diff2D(1,1);
    
    for(y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
	SrcIterator  sx = sul;
	DestIterator dx = dul;
	
	for(x=0; x<w; ++x, ++sx.x, ++dx.x)
	{
	    for(i=0; i<8; ++i)
	    {
		if(!(sa(sx, dist[i]) < sa(sx))) break;
	    }
	
	    if(i == 8) da.set(marker, dx);
	}
    }
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
localMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    localMaxima(sul, slr, sa, dul, da, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
inline void 
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
	    DestValue marker)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
localMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    localMaxima(src.first, src.second, src.third,
                dest.first, dest.second, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

/********************************************************/
/*                                                      */
/*                 extendedLocalMinima                  */
/*                                                      */
/********************************************************/

/** \brief Find local minimal regions in an image.

    This function finds regions of uniform pixel value
    whose neighboring regions are all have larger values
    (minimal plateaus of arbitrary size). Minimal regions are
    marked in the destination image with the given marker value
    (default is 1), all other destination pixels remain unchanged.
    <TT>SrcAccessor::value_type</TT> must be equality-comparable and
    less-comparable.
    A pixel at the image border will never be marked as minimum or 
    minimal plateau. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
	          class DestIterator, class DestAccessor, 
	          class DestValue = DestAccessor::value_type>
        void 
        extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
		            DestIterator dul, DestAccessor da, 
		            DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
	          class DestIterator, class DestAccessor, 
	          class DestValue = DestAccessor::value_type>
        void 
        extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
		            pair<DestIterator, DestAccessor> dest,
		            DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BImage src(w,h), minima(w,h);
    
    // init destiniation image
    minima = 0;
    
    vigra::extendedLocalMinima(srcImageRange(src), destImage(minima));
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);
    
    u == u
    u < u
    
    DestValue marker;
    dest_accessor.set(marker, dest_upperleft);
    \endcode

*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
void 
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    static const Diff2D dist[] = {
        Diff2D( 1, 0), Diff2D( 1, -1), Diff2D( 0, -1), Diff2D( -1, -1), 
        Diff2D( -1, 0), Diff2D( -1, 1), Diff2D( 0, 1), Diff2D( 1, 1),
        Diff2D( 1, 0), Diff2D( 1, -1), Diff2D( 0, -1), Diff2D( -1, -1)};
    
    int i,x,y;
    
    BasicImage<int> labels(w,h);
    
    labels = 0;
    
    int number_of_regions = 
        labelImage(sul, slr, sa, labels.upperLeft(), labels.accessor(), true);

    std::vector<unsigned char> processed(number_of_regions+1, (unsigned char)0);
    
    SrcIterator ys = sul + Diff2D(1,1);
    BasicImage<int>::Iterator lul = labels.upperLeft();
    BasicImage<int>::Iterator ly = lul + Diff2D(1,1);
    
    initImage(lul, lul+Diff2D(1,h), labels.accessor(), 0);

    for(y=1; y<h-1; ++y, ++ys.y, ++ly.y)
    {
	SrcIterator  sx = ys;
	BasicImage<int>::Iterator lx(ly);
	
	for(x=1; x<w-1; ++x, ++sx.x, ++lx.x)
	{
	    int lab = *lx;
	    typename SrcAccessor::value_type v = sa(sx);
	    
	    if(processed[lab]) continue;
	    
	    processed[lab] = 1; // assume minimum until opposite is proved
	    
	    int is_plateau = 0;
	    
	    for(i=11; i>3; --i)
	    {
		if(lx[dist[i]] == lab)
		{
		    is_plateau = i;
		}
		else if(sa(sx, dist[i]) < v) 
		{
		    break;
		}
	    }

	    if(i > 3) 
	    {
	        processed[lab] = 2;  // not a minimum
		continue;
	    }
	    
	    if(!is_plateau) continue;  // is a minimum
	    
	    if((x == 1) && (is_plateau == 4) &&
	       (lx[dist[3]] == lab)) is_plateau = 3;
	    
	    // is a plateau - do contour-following
	    int xa = x;
	    int ya = y;
	    int first_dir = is_plateau & 7;
	    int dir = first_dir;
	    int xc = xa;
	    int yc = ya;
	    
	    do
	    {
		xc = xc + dist[dir].x;
		yc = yc + dist[dir].y;
		
		dir = (dir + 6) & 7;
		
		for (; true; dir = (dir + 1) & 7)
		{
		    int xn = xc + dist[dir].x;
		    int yn = yc + dist[dir].y;
                    Diff2D dn(xn, yn);
		    
		    if((xn >= 0) && (xn < w) && (yn >= 0) && (yn < h))
		    {
		        if(lul[dn] == lab) break;
			
			if(dir & 1) continue;
			
			if(sa(sul, dn) < v) 
			{
			    processed[lab] = 2; // current region not a minimum
			}
			else
			{
			    processed[lul[dn]] = 2; // other region not 
			                                // a minimum
			}
		    }
		}
            }
	    while((xc != xa) || (yc != ya) || (dir != first_dir));
	}
    }

    for(y=0; y<h; ++y, ++dul.y, ++lul.y)
    {
	DestIterator  xd = dul;
	BasicImage<int>::Iterator lx(lul);
	
	for(x=0; x<w; ++x, ++xd.x, ++lx.x)
	{
	    if(processed[*lx] == 1) da.set(marker, xd);
	}
    }
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
extendedLocalMinima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    extendedLocalMinima(sul, slr, sa, dul, da, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
inline void 
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
	    DestValue marker)
{
    extendedLocalMinima(src.first, src.second, src.third,
                dest.first, dest.second, marker);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
extendedLocalMinima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    extendedLocalMinima(src.first, src.second, src.third,
                dest.first, dest.second, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

/********************************************************/
/*                                                      */
/*                 extendedLocalMaxima                  */
/*                                                      */
/********************************************************/

/** \brief Find local maximal regions in an image.

    This function finds regions of uniform pixel value
    whose neighboring regions are all have smaller values
    (maximal plateaus of arbitrary size). Maximal regions are
    marked in the destination image with the given marker value
    (default is 1), all other destination pixels remain unchanged.
    <TT>SrcAccessor::value_type</TT> must be equality-comparable and
    less-comparable.
    A pixel at the image border will never be marked as maximum or 
    maximal plateau. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
	          class DestIterator, class DestAccessor, 
	          class DestValue = DestAccessor::value_type>
        void 
        extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
		            DestIterator dul, DestAccessor da, 
		            DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
	          class DestIterator, class DestAccessor, 
	          class DestValue = DestAccessor::value_type>
        void 
        extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
		            pair<DestIterator, DestAccessor> dest,
		            DestValue marker = NumericTraits<DestValue>::one())
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="localminmax_8hxx-source.html">vigra/localminmax.hxx</a>"<br>
    Namespace: vigra
    
    \code
    vigra::BImage src(w,h), maxima(w,h);
    
    // init destiniation image
    maxima = 0;
    
    vigra::extendedLocalMaxima(srcImageRange(src), destImage(maxima));
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    
    SrcAccessor::value_type u = src_accessor(src_upperleft);
    
    u == u
    u < u
    
    DestValue marker;
    dest_accessor.set(marker, dest_upperleft);
    \endcode

*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
void 
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da, DestValue marker)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    static const Diff2D dist[] = {
        Diff2D( 1, 0), Diff2D( 1, -1), Diff2D( 0, -1), Diff2D( -1, -1), 
        Diff2D( -1, 0), Diff2D( -1, 1), Diff2D( 0, 1), Diff2D( 1, 1),
        Diff2D( 1, 0), Diff2D( 1, -1), Diff2D( 0, -1), Diff2D( -1, -1)};
    
    int i,x,y;
    
    BasicImage<int> labels(w,h);
    
    int number_of_regions = 
        labelImage(sul, slr, sa, labels.upperLeft(), labels.accessor(), true);

    std::vector<unsigned char> processed(number_of_regions+1, (unsigned char)0);
    
    SrcIterator ys = sul + Diff2D(1,1);
    BasicImage<int>::Iterator lul = labels.upperLeft();
    BasicImage<int>::Iterator ly = lul + Diff2D(1,1);

    lul(0,1) = 0;

    for(y=1; y<h-1; ++y, ++ys.y, ++ly.y)
    {
	SrcIterator  sx = ys;
	BasicImage<int>::Iterator lx(ly);
	
	for(x=1; x<w-1; ++x, ++sx.x, ++lx.x)
	{
	    int lab = *lx;
	    typename SrcAccessor::value_type v = sa(sx);
	    
	    if(processed[lab]) continue;
	    
	    processed[lab] = 1; // assume maximum until opposite is proved
	    
	    int is_plateau = 0;
	    
	    for(i=11; i>3; --i)
	    {
		if(lx[dist[i]] == lab)
		{
		    is_plateau = i;
		}
		else if(v < sa(sx, dist[i])) 
		{
		    break;
		}
	    }

	    if(i > 3) 
	    {
	        processed[lab] = 2;  // not a maximum
		continue;
	    }
	    
	    if(!is_plateau) continue;  // is a maximum
	    
	    if((x == 1) && (is_plateau == 4) &&
	       (lx[dist[3]] == lab)) is_plateau = 3;
	    
	    // is a plateau - do contour-following
	    int xa = x;
	    int ya = y;
	    int first_dir = is_plateau & 7;
	    int dir = first_dir;
	    int xc = xa;
	    int yc = ya;
	    
	    do
	    {
		xc = xc + dist[dir].x;
		yc = yc + dist[dir].y;
		
		dir = (dir + 6) & 7;
		
		for (; true; dir = (dir + 1) & 7)
		{
		    int xn = xc + dist[dir].x;
		    int yn = yc + dist[dir].y;
                    Diff2D dn(xn, yn);
		    
		    if((xn >= 0) && (xn < w) && (yn >= 0) && (yn < h))
		    {
		        if(lul[dn] == lab) break;
			
			if(dir & 1) continue;
			
			if(v < sa(sul, dn)) 
			{
			    processed[lab] = 2; // current region not a maximum
			}
			else
			{
			    processed[lul[dn]] = 2; // other region not 
			                                // a maximum
			}
		    }
		}
            }
	    while((xc != xa) || (yc != ya) || (dir != first_dir));
	}
    }

    for(y=0; y<h; ++y, ++dul.y, ++lul.y)
    {
	DestIterator  xd = dul;
	BasicImage<int>::Iterator lx(lul);
	
	for(x=0; x<w; ++x, ++xd.x, ++lx.x)
	{
	    if(processed[*lx] == 1) da.set(marker, xd);
	}
    }
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
extendedLocalMaxima(SrcIterator sul, SrcIterator slr, SrcAccessor sa,
            DestIterator dul, DestAccessor da)
{
    extendedLocalMaxima(sul, slr, sa, dul, da, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor, class DestValue>
inline void 
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest,
	    DestValue marker)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                dest.first, dest.second, marker);
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline void 
extendedLocalMaxima(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    extendedLocalMaxima(src.first, src.second, src.third,
                dest.first, dest.second, 
                NumericTraits<typename DestAccessor::value_type>::one());
}

//@}

} // namespace vigra

#endif // VIGRA_LOCALMINMAX_HXX
