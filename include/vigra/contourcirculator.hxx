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


#ifndef VIGRA_CONTOURCIRCULATOR_HXX
#define VIGRA_CONTOURCIRCULATOR_HXX

#include <vigra/pixelneighborhood.hxx>

namespace vigra
{

/** \addtogroup ImageIteratorAdapters
 */
//@{

/********************************************************/
/*                                                      */
/*                CrackContourCirculator                */
/*                                                      */
/********************************************************/

/** \brief Circulator that walks around a given region.

    The circulator follows the <em>crack contour</em> of a given region.
    Here, a region is an 8-connected component of pixels with the same
    value, such as the regions in a label image. 
    The crack contour is located between the inside and outside
    pixels, that is "on the crack" between the region and the background. 
    Thus, the circulator moves from pixel corner to pixel corner. By definition,
    the first corner (were the circulator was initialized) gets the
    coordinate (0,0), and calls to <tt>*circulator</tt> return the distance
    of the current corner to the initial one. 
    
    The circulator can be used to calculate the area of the region (in pixels):
    
    \code
    ImageIterator region_anchor = ...;
    int area = 0;
    
    // calculate area from following the crack contour of the region
    CrackContourCirculator<ImageIterator> crack(region_anchor);
    CrackContourCirculator<ImageIterator> crackend(crack);
    do
    {
        area += crack.diff().x * crack.pos().y -
                crack.diff().y * crack.pos().x;
    }
    while(++crack != crackend);

    area /= 2;
    std::cout << "Area of region " << *region_anchor << ": " << area << std::endl;
    \endcode

    <b>\#include</b> "<a href="contourcirculator_8hxx-source.html">vigra/contourcirculator.hxx</a>"<br>
    Namespace: vigra
*/
template <class IMAGEITERATOR>
class CrackContourCirculator
{
    typedef NeighborhoodCirculator<IMAGEITERATOR, EightNeighborCode> 
            NEIGHBORHOODCIRCULATOR;
    typedef typename IMAGEITERATOR::value_type label_type;

protected:
    NEIGHBORHOODCIRCULATOR neighborCirc_;
    label_type label_;
    Diff2D pos_;

public:
        /** the circulator's value type
        */
    typedef Diff2D value_type;

        /** the circulator's reference type (return type of <TT>*circ</TT>)
        */
    typedef Diff2D const & reference;

        /** the circulator's pointer type (return type of <TT>operator-></TT>)
        */
    typedef Diff2D const * pointer;

        /** the circulator tag
        */
    typedef forward_circulator_tag iterator_category;

        /** Initialize the circulator for a given region.
        
            The image iterator <tt>in_the_region</tt> must refer
            to a boundary pixel of the region to be analysed. The
            direction code <tt>dir</tt> must point to a pixel outside the 
            region (the default assumes that the pixel left of the 
            given region pixel belongs to the background).    
            The first corner of the crack contour is the corner to the 
            right of this direction (i.e. the north west corner of 
            the region pixel, if the direction was West).
        */
    CrackContourCirculator(IMAGEITERATOR const & in_the_region, 
                           FourNeighborCode::Direction dir = FourNeighborCode::West)
        : neighborCirc_(in_the_region, EightNeighborCode::code(dir)),
          label_(*(neighborCirc_.center())),
          pos_(0, 0)
    {
        neighborCirc_.turnLeft();
    }

        /** Move to the next crack corner of the contour (pre-increment).
        */
    CrackContourCirculator & operator++()
    {
        pos_ += neighborCirc_.diff();

        neighborCirc_--;

        if(*neighborCirc_ == label_)
        {
            neighborCirc_.moveCenterToNeighbor(); // TODO: simplify moveCenterToNeighbor()s
            --neighborCirc_;
        }
        else
        {
            neighborCirc_.moveCenterToNeighbor(); // jump out
            neighborCirc_ += 3;
            if(*neighborCirc_ == label_)
            {
                neighborCirc_.moveCenterToNeighbor();
                neighborCirc_.turnRight();
            }
            else
            {
                neighborCirc_.moveCenterToNeighbor();
                neighborCirc_.turnLeft();
                neighborCirc_.moveCenterToNeighbor();
                neighborCirc_.turnRight();
            }
        }

        return *this;
    }

        /** Move to the next crack corner of the contour (post-increment).
        */
    CrackContourCirculator operator++(int)
    {
        CrackContourCirculator ret(*this);
        ++(*this);
        return ret;
    }
    
        /** equality
        */
    bool operator==(CrackContourCirculator const & o) const
    {
        return neighborCirc_ == o.neighborCirc_;
    }

        /** inequality
        */
    bool operator!=(CrackContourCirculator const & o) const
    {
        return neighborCirc_ != o.neighborCirc_;
    }

        /** Get the coordinate of the current corner 
            (relative to the first corner).
        */
    reference pos() const 
        { return pos_; }

        /** Equivalent to pos()
        */
    reference operator*() const 
        { return pos_; }

        /** Access member of the current coordinate.
        */
    pointer operator->() const 
        { return &pos_; }

        /** Get the offset from the current corner of the contour 
            to the next one.
        */
    reference diff() const 
        { return neighborCirc_.diff(); }
};

//@}

} // namespace vigra

#endif /* VIGRA_CONTOURCIRCULATOR_HXX */
