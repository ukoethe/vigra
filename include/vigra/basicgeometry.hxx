/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
 
#ifndef VIGRA_BASICGEOMETRY_HXX
#define VIGRA_BASICGEOMETRY_HXX

#include "vigra/error.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/copyimage.hxx"
#include <cmath>

namespace vigra {

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
void rotateImage(SrcIterator is, SrcIterator end, SrcAccessor as,
                           DestIterator id, DestAccessor ad, int rotation)
{
    int x, y;
    int ws = end.x - is.x;
    int hs = end.y - is.y;

    vigra_precondition(rotation % 90 == 0, 
                "rotateImage(): "
                "This function rotates images only about multiples of 90 degree");

    rotation = rotation%360; 
    if (rotation < 0)
        rotation += 360;
    
    switch(rotation)
    {
        case 0:
            copyImage(is, end, as, id, ad);
            break;
        case 90: 
            is.x += (ws-1);
            for(x=0; x != ws; x++, is.x--, id.y++)
            {
                typename SrcIterator::column_iterator cs = is.columnIterator();
                typename DestIterator::row_iterator rd = id.rowIterator();
                for(y=0; y != hs; y++, cs++, rd++)
                {
                    ad.set(as(cs), rd);
                }
        
            }
            break;

        case 180:
            end.x--;
            end.y--;
            for(x=0; x != ws; x++, end.x--, id.x++)
            {
                typename SrcIterator::column_iterator cs = end.columnIterator();
                typename DestIterator::column_iterator cd = id.columnIterator();
                for(y=0; y != hs; y++, cs--, cd++)
                {
                    ad.set(as(cs), cd);
                }
        
            }
            break;

        case 270:  
            is.y += (hs-1);
            for(x=0; x != ws; x++, is.x++, id.y++)
            {
                typename SrcIterator::column_iterator cs = is.columnIterator();
                typename DestIterator::row_iterator rd = id.rowIterator();
                for(y=0; y != hs; y++, cs--, rd++)
                {
                    ad.set(as(cs), rd);
                }
        
            }
            break;
        default: //not needful, because of the exception handig in if-statement 
            vigra_fail("internal error"); 
    }
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void 
rotateImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
              pair<DestImageIterator, DestAccessor> dest, int rotation)
{
    rotateImage(src.first, src.second, src.third, dest.first, dest.second, rotation);
}

/*
*  Flip turns the image upside down.
*  Flop reverses the image.
*/

enum Reflect{horizontal = 1, vertical = 2};

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
void reflectImage(SrcIterator is, SrcIterator end, SrcAccessor as,
                  DestIterator id, DestAccessor ad, Reflect reflect)
{
    
    int ws = end.x - is.x;
    int hs = end.y - is.y;

    int x, y;

    if(reflect == horizontal)
    {//flipImage
        is.y += (hs-1);
        for(x=0; x<ws; ++x, ++is.x, ++id.x) 
        {
            typename SrcIterator::column_iterator  cs = is.columnIterator();
            typename DestIterator::column_iterator cd = id.columnIterator();
            for(y=0; y!=hs;y++, cs--, cd++)
            {
                ad.set(as(cs), cd);
            }
        }
    }
    else if(reflect == vertical)
    {//flopImage
        is.x += (ws-1);
        for(x=0; x < ws; ++x, --is.x, ++id.x) 
        {

            typename SrcIterator::column_iterator cs = is.columnIterator();
            typename DestIterator::column_iterator cd = id.columnIterator();
            for(y=0; y!=hs;y++, cs++, cd++)
            {
                ad.set(as(cs), cd);
            }
        }
    }
    else if(reflect == (horizontal | vertical))
    {//flipFlopImage   //???
        end.x--;
        end.y--;
        for(x=0; x != ws; x++, end.x--, id.x++)
        {
            typename SrcIterator::column_iterator cs = end.columnIterator();
            typename DestIterator::column_iterator cd = id.columnIterator();
            for(y=0; y != hs; y++, cs--, cd++)
            {
                ad.set(as(cs), cd);
            }
        }
    }
    else 
        vigra_fail("reflectImage(): "
                   "This function reflects horizontal or vertical,"
                   "   'and' is included");
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void 
reflectImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
              pair<DestImageIterator, DestAccessor> dest, Reflect reflect)
{
    reflectImage(src.first, src.second, src.third, dest.first, dest.second, reflect);
}

enum Transpose{major = 1, minor = 2};

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
void transposeImage(SrcIterator is, SrcIterator end, SrcAccessor as,
                    DestIterator id, DestAccessor ad, Transpose transpose)
{
    int ws = end.x - is.x;
    int hs = end.y - is.y;

    int x, y;

    if(transpose == major)
    {//Die Funktion spiegelt das Bild um (0,0) (1,1) Diagonale
        for(x=0; x != ws; x++, is.x++, id.y++)
        {

            typename SrcIterator::column_iterator cs = is.columnIterator();
            typename DestIterator::row_iterator rd = id.rowIterator();
            for(y=0; y != hs; y++, cs++, rd++)
            {
                ad.set(as(cs), rd);
            }
        }
    }
    else if(transpose == minor)
    {//Die Funktion spiegelt das Bild (1,0) (0,1) Diagonale
        end.x--;
        end.y--;
        for(x=0; x != ws; x++, --end.x, ++id.y)
        {

            typename SrcIterator::column_iterator cs = end.columnIterator();
            typename DestIterator::row_iterator rd = id.rowIterator();
            for(y=0; y != hs; y++, --cs, ++rd)
            {
                ad.set(as(cs), rd);
            }
        }
    }
    else if(transpose == (major | minor))
    {//flipFlopImage  //???
        end.x--;
        end.y--;
        for(x=0; x != ws; x++, end.x--, id.x++)
        {
            typename SrcIterator::column_iterator cs = end.columnIterator();
            typename DestIterator::column_iterator cd = id.columnIterator();
            for(y=0; y != hs; y++, cs--, cd++)
            {
                ad.set(as(cs), cd);
            }
        }
    
    }
    else 
        vigra_fail("transposeImage(): "
                   "This function transposes major or minor,"
                   "   'and' is included");

}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void 
transposeImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
              pair<DestImageIterator, DestAccessor> dest, Transpose transpose)
{
    transposeImage(src.first, src.second, src.third, dest.first, dest.second, transpose);
}

/********************************************************/
/*                                                      */
/*                        resampleLine                  */
/*                                                      */
/********************************************************/


/**
* Vergroessert eine Linie um einen Faktor. 
* Ist z.B. der Faktor = 4 so werden in der
* neuen Linie(Destination) jedes Pixel genau 4 mal 
* vorkommen, also es findet auch keine Glaetung 
* statt (NoInterpolation). Als Parameter sollen
* der Anfangs-, der Enditerator und der Accessor
* der Ausgangslinie (Source line), der Anfangsiterator
* und Accessor der Ziellinie (destination line) und
* anschliessend der Faktor um den die Linie (Zeile)
* vergroessert bzw. verkleinert werden soll. 
*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
void resampleLine(SrcIterator src_iter, SrcIterator src_iter_end, SrcAccessor src_acc,
                  DestIterator dest_iter, DestAccessor dest_acc, double factor)
{
    // The width of the src line.      
    int src_width = src_iter_end - src_iter;
 
    vigra_precondition(src_width > 0,
                       "resampleLine(): input image too small.");
    vigra_precondition(factor > 0.0,
                       "resampleLine(): factor must be positive.");
    
    if (factor >= 1.0)
    {
        int int_factor = (int)factor;
        double dx = factor - int_factor;
        double saver = dx;
        for ( ; src_iter != src_iter_end ; ++src_iter, saver += dx)
        {
            if (saver >= 1.0)
            {
                saver = saver - (int)saver;
                dest_acc.set(src_acc(src_iter), dest_iter);
                ++dest_iter;
            }
            for(int i = 0 ; i < int_factor ; i++, ++dest_iter)
            {
                dest_acc.set(src_acc(src_iter), dest_iter);
            }
        }
    }
    else
    {
        DestIterator dest_end = dest_iter + (int)VIGRA_CSTD::ceil(src_width*factor);  
        factor = 1.0/factor;
        int int_factor = (int)factor;
        double dx = factor - int_factor;
        double saver = dx;
        src_iter_end -= 1;
        for ( ; src_iter != src_iter_end && dest_iter != dest_end ; 
              ++dest_iter, src_iter += int_factor, saver += dx)
        {
            if (saver >= 1.0)
            {
                saver = saver - (int)saver;
                ++src_iter;
            }
            dest_acc.set(src_acc(src_iter), dest_iter);
        }
        if (dest_iter != dest_end)
        {
            dest_acc.set(src_acc(src_iter_end), dest_iter);
        }
    }
}


/********************************************************/
/*                                                      */
/*                     resampleImage                    */
/*                                                      */
/********************************************************/

inline int sizeForResamplingFactor(int oldsize, double factor)
{
    return (factor < 1.0)
        ? (int)VIGRA_CSTD::ceil(oldsize * factor) 
        : (int)(oldsize * factor);
}

/** \brief Resample image by a given factor.

    This algorithm is very fast and does not require any arithmetic on the pixel types.
    
    The range must of the input image (resp. regions)
    must be given and it must have a size of at
    least 2x2. Destiniation pixels are directly copied from the appropriate
    source pixels. 
    The function uses accessors. 
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void 
        resampleImage(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                      DestIterator id, DestAccessor ad, double factor);
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
              class DestImageIterator, class DestAccessor>
        inline void 
        resampleImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                      pair<DestImageIterator, DestAccessor> dest, double factor);
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="basicgeometry_8hxx-source.html">vigra/basicgeometry.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::resampleImage(srcImageRange(src), destImage(dest), factor);
    
    \endcode

    <b> Required Interface:</b>
    
    \code
    SrcImageIterator src_upperleft, src_lowerright;
    DestImageIterator dest_upperleft;
    
    SrcAccessor src_accessor;
    
    dest_accessor.set(src_accessor(src_upperleft), dest_upperleft);

    \endcode
    
    <b> Preconditions:</b>
    
    \code
    src_lowerright.x - src_upperleft.x > 1
    src_lowerright.y - src_upperleft.y > 1
    \endcode
    
*/
template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
resampleImage(SrcIterator is, SrcIterator iend, SrcAccessor sa,
              DestIterator id, DestAccessor ad, double factor)
{
    int width_old = iend.x - is.x;
    int height_old = iend.y - is.y;
    
    //Bei Verkleinerung muss das dest-Bild ceiling(src*factor), da z.B.
    //aus 6x6 grossem Bild wird eins 18x18 grosses gemacht bei Vergroesserungsfaktor 3.1
    //umgekehrt damit wir vom 18x18 zu 6x6 (und nicht 5x5) bei Vergroesserung von 1/3.1
    //muss das kleinste Integer das groesser als 18/3.1 ist genommen werden.
    int height_new = sizeForResamplingFactor(height_old, factor);
    int width_new = sizeForResamplingFactor(width_old, factor);
    
    vigra_precondition((width_old > 1) && (height_old > 1),
                 "resampleImage(): "
                 "Source image to small.\n");
    vigra_precondition((width_new > 1) && (height_new > 1),
                 "resampleImage(): "
                 "Destination image to small.\n");
        
    typedef typename SrcAccessor::value_type SRCVT;
    typedef BasicImage<SRCVT> TmpImage;
    typedef typename TmpImage::traverser TmpImageIterator;

    BasicImage<SRCVT> tmp(width_old, height_new);
    
    int x,y;
    
    typename BasicImage<SRCVT>::Iterator yt = tmp.upperLeft();

    for(x=0; x<width_old; ++x, ++is.x, ++yt.x) 
    {
        typename SrcIterator::column_iterator c1 = is.columnIterator();
        typename TmpImageIterator::column_iterator ct = yt.columnIterator();
        resampleLine(c1, c1 + height_old, sa, ct, tmp.accessor(), factor);
    }

    yt = tmp.upperLeft();

    for(y=0; y < height_new; ++y, ++yt.y, ++id.y) 
    {
        typename DestIterator::row_iterator rd = id.rowIterator();
        typename TmpImageIterator::row_iterator rt = yt.rowIterator();
        resampleLine(rt, rt + width_old, tmp.accessor(), rd, ad, factor);
    }

}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void 
resampleImage(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
              pair<DestImageIterator, DestAccessor> dest, double factor)
{
    resampleImage(src.first, src.second, src.third, dest.first, dest.second, factor);
}


} // namespace vigra


#endif /* VIGRA_BASICGEOMETRY_HXX */
