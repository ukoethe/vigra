/************************************************************************/
/*                                                                      */
/*               Copyright 2002-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_TENSORUTILITIES_HXX
#define VIGRA_TENSORUTILITIES_HXX

#include <cmath>
#include "vigra/utilities.hxx"

namespace vigra {

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void vectorToTensor(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                    DestIterator dul, DestAccessor dest)
{
    vigra_precondition(src.size(sul) == 2,
                       "vectorToTensor(): input image must have 2 bands.");
    vigra_precondition(dest.size(dul) == 3,
                       "vectorToTensor(): output image must have 3 bands.");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    for(int y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        typename SrcIterator::row_iterator s = sul.rowIterator();
        typename SrcIterator::row_iterator send = s + w;
        typename DestIterator::row_iterator d = dul.rowIterator();
        for(; s < send; ++s, ++d)
        {
            dest.setComponent(sq(src.getComponent(s, 0)), d, 0);
            dest.setComponent(-src.getComponent(s, 0)*src.getComponent(s, 1), d, 1);
                           // ^ negative sign to turn left-handed into right-handed coordinates
            dest.setComponent(sq(src.getComponent(s, 1)), d, 2);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void vectorToTensor(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                     pair<DestIterator, DestAccessor> d)
{
    vectorToTensor(s.first, s.second, s.third, d.first, d.second);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void tensorEigenRepresentation(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                               DestIterator dul, DestAccessor dest)
{
    vigra_precondition(src.size(sul) == 3,
                       "tensorEigenRepresentation(): input image must have 3 bands.");
    vigra_precondition(dest.size(dul) == 3,
                       "tensorEigenRepresentation(): output image must have 3 bands.");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    for(int y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        typename SrcIterator::row_iterator s = sul.rowIterator();
        typename SrcIterator::row_iterator send = s + w;
        typename DestIterator::row_iterator d = dul.rowIterator();
        for(; s < send; ++s, ++d)
        {
            typedef typename 
               NumericTraits<typename SrcAccessor::component_type>::RealPromote TmpType;
            TmpType d1 = src.getComponent(s,0) + src.getComponent(s,2);
            TmpType d2 = src.getComponent(s,0) - src.getComponent(s,2);
            TmpType d3 = 2.0 * src.getComponent(s,1);
            TmpType d4 = VIGRA_CSTD::sqrt(sq(d2) + sq(d3));
            
            dest.setComponent(0.5 * (d1 + d4), d, 0); // large EV
            dest.setComponent(0.5 * (d1 - d4), d, 1); // small EV
            dest.setComponent(0.5 * VIGRA_CSTD::atan2(d3, d2), d, 2); // orientation
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void tensorEigenRepresentation(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                               pair<DestIterator, DestAccessor> d)
{
    tensorEigenRepresentation(s.first, s.second, s.third, d.first, d.second);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void tensorTrace(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                 DestIterator dul, DestAccessor dest)
{
    vigra_precondition(src.size(sul) == 3,
                       "tensorTrace(): input image must have 3 bands.");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    for(int y=0; y<h; ++y, ++sul.y, ++dul.y)
    {
        typename SrcIterator::row_iterator s = sul.rowIterator();
        typename SrcIterator::row_iterator send = s + w;
        typename DestIterator::row_iterator d = dul.rowIterator();
        for(; s < send; ++s, ++d)
        {
            dest.set(src.getComponent(s,0) + src.getComponent(s,2), d);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void tensorTrace(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                 pair<DestIterator, DestAccessor> d)
{
    tensorTrace(s.first, s.second, s.third, d.first, d.second);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator1, class DestAccessor1,
          class DestIterator2, class DestAccessor2>
void tensorToEdgeCorner(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                        DestIterator1 edgeul, DestAccessor1 edge,
                        DestIterator2 cornerul, DestAccessor2 corner)
{
    vigra_precondition(src.size(sul) == 3,
                       "tensorToEdgeCorner(): input image must have 3 bands.");
    vigra_precondition(edge.size(edgeul) == 2,
                       "tensorToEdgeCorner(): edge image must have 2 bands.");

    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    for(int y=0; y<h; ++y, ++sul.y, ++edgeul.y, ++cornerul.y)
    {
        typename SrcIterator::row_iterator s = sul.rowIterator();
        typename SrcIterator::row_iterator send = s + w;
        typename DestIterator1::row_iterator e = edgeul.rowIterator();
        typename DestIterator2::row_iterator c = cornerul.rowIterator();
        for(; s < send; ++s, ++e, ++c)
        {
            typedef typename 
               NumericTraits<typename SrcAccessor::component_type>::RealPromote TmpType;
            TmpType d1 = src.getComponent(s,0) + src.getComponent(s,2);
            TmpType d2 = src.getComponent(s,0) - src.getComponent(s,2);
            TmpType d3 = 2.0 * src.getComponent(s,1);
            TmpType d4 = VIGRA_CSTD::sqrt(sq(d2) + sq(d3));
            
            edge.setComponent(d4, e, 0); // edgeness = difference of EVs
            edge.setComponent(0.5 * VIGRA_CSTD::atan2(d3, d2), e, 1); // orientation
            corner.set(d1 - d4, c); // cornerness = 2 * small EV
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator1, class DestAccessor1,
          class DestIterator2, class DestAccessor2>
inline
void tensorToEdgeCorner(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                        pair<DestIterator1, DestAccessor1> edge,
                        pair<DestIterator2, DestAccessor2> corner)
{
    tensorToEdgeCorner(s.first, s.second, s.third, 
                       edge.first, edge.second, corner.first, corner.second);
}

} // namespace vigra

#endif /* VIGRA_TENSORUTILITIES_HXX */
