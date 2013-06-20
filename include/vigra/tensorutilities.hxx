/************************************************************************/
/*                                                                      */
/*               Copyright 2002-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_TENSORUTILITIES_HXX
#define VIGRA_TENSORUTILITIES_HXX

#include <cmath>
#include "utilities.hxx"
#include "mathutil.hxx"
#include "multi_shape.hxx"

namespace vigra {

/** \addtogroup TensorImaging Tensor Image Processing
*/
//@{

/********************************************************/
/*                                                      */
/*                      vectorToTensor                  */
/*                                                      */
/********************************************************/

/** \brief Calculate the tensor (outer) product of a 2D vector with itself.

    This function is useful to transform vector images into a tensor representation 
    that can be used as input to tensor based processing and analysis functions
    (e.g. tensor smoothing). The input pixel type must be vectors of length 2, whereas
    the output must contain vectors of length 3 which will represent the tensor components
    in the order t11, t12 (== t21 due to symmetry), t22.
    
    <b>Note:</b> In order to account for the left-handedness of the image coordinate system,
    the second tensor component (t12) can be negated by setting <tt>negateComponent2 = false</tt>.
    Angles will then be interpreted counter-clockwise rather than clockwise. By default,
    this behavior is switched off.
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        vectorToTensor(MultiArrayView<2, T1, S1> const & src,
                       MultiArrayView<2, T2, S2> dest,
                       bool negateComponent2 = false);
    }
    \endcode

    \deprecatedAPI{vectorToTensor}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void vectorToTensor(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                            DestIterator dul, DestAccessor dest,
                            bool negateComponent2 = false);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void vectorToTensor(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                            pair<DestIterator, DestAccessor> d,
                            bool negateComponent2 = false);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/tensorutilities.hxx\>

    \code
    FImage img(w,h);
    FVector2Image gradient(w,h);
    FVector3Image tensor(w,h);
    
    gaussianGradient(srcImageRange(img), destImage(gradient), 2.0);
    vectorToTensor(srcImageRange(gradient), destImage(tensor));
    \endcode

*/
doxygen_overloaded_function(template <...> void vectorToTensor)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void vectorToTensor(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                    DestIterator dul, DestAccessor dest,
                    bool negateComponent2)
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
        if(negateComponent2)
        {
            for(; s < send; ++s, ++d)
            {
                dest.setComponent(sq(src.getComponent(s, 0)), d, 0);
                dest.setComponent(-src.getComponent(s, 0)*src.getComponent(s, 1), d, 1);
                               // ^ negative sign to turn left-handed into right-handed coordinates
                dest.setComponent(sq(src.getComponent(s, 1)), d, 2);
            }
        }
        else
        {
            for(; s < send; ++s, ++d)
            {
                dest.setComponent(sq(src.getComponent(s, 0)), d, 0);
                dest.setComponent(src.getComponent(s, 0)*src.getComponent(s, 1), d, 1);
                dest.setComponent(sq(src.getComponent(s, 1)), d, 2);
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void vectorToTensor(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                    DestIterator dul, DestAccessor dest)
{
    vectorToTensor(sul, slr, src, dul, dest, false);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
vectorToTensor(triple<SrcIterator, SrcIterator, SrcAccessor> s,
               pair<DestIterator, DestAccessor> d,
               bool negateComponent2)
{
    vectorToTensor(s.first, s.second, s.third, d.first, d.second, negateComponent2);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
vectorToTensor(triple<SrcIterator, SrcIterator, SrcAccessor> s,
               pair<DestIterator, DestAccessor> d)
{
    vectorToTensor(s.first, s.second, s.third, d.first, d.second, false);
}

template <class T1, class S1,
          class T2, class S2>
inline void
vectorToTensor(MultiArrayView<2, T1, S1> const & src,
               MultiArrayView<2, T2, S2> dest,
               bool negateComponent2 = false)
{
    vigra_precondition(src.shape() == dest.shape(),
        "vectorToTensor(): shape mismatch between input and output.");
    vectorToTensor(srcImageRange(src), destImage(dest), negateComponent2);
}

/********************************************************/
/*                                                      */
/*               tensorEigenRepresentation              */
/*                                                      */
/********************************************************/

/** \brief Calculate eigen representation of a symmetric 2x2 tensor.

    This function turns a 3-band image representing the tensor components
    t11, t12 (== t21 due to symmetry), t22 into the a 3-band image holding the eigen
    representation e1, e2, and angle, where e1 \> e2. When the tensor is
    defined in a left-handed coordinate system (the default on images), the angle will
    then be given in clockwise orientation, starting at the x-axis. Otherwise, it
    will be given in counter-clockwise orientation.
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        tensorEigenRepresentation(MultiArrayView<2, T1, S1> const & src,
                                  MultiArrayView<2, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{tensorEigenRepresentation}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void tensorEigenRepresentation(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                       DestIterator dul, DestAccessor dest);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void tensorEigenRepresentation(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                                       pair<DestIterator, DestAccessor> d);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/tensorutilities.hxx\>

    \code
    FVector3Image tensor(w,h);
    FVector3Image eigen(w,h);
    
    tensorEigenRepresentation(srcImageRange(tensor), destImage(eigen));
    \endcode

*/
doxygen_overloaded_function(template <...> void tensorEigenRepresentation)

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
            TmpType d3 = TmpType(2.0) * src.getComponent(s,1);
            TmpType d4 = (TmpType)hypot(d2, d3);
            
            dest.setComponent(0.5 * (d1 + d4), d, 0); // large EV
            dest.setComponent(0.5 * (d1 - d4), d, 1); // small EV
            if(d2==0.0 && d3==0.0)
            {
                dest.setComponent(0, d, 2); // orientation
            }
            else
            {
                dest.setComponent(0.5 * VIGRA_CSTD::atan2(d3, d2), d, 2); // orientation
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
tensorEigenRepresentation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest)
{
    tensorEigenRepresentation(src.first, src.second, src.third, dest.first, dest.second);
}

template <class T1, class S1,
          class T2, class S2>
inline void
tensorEigenRepresentation(MultiArrayView<2, T1, S1> const & src,
                          MultiArrayView<2, T2, S2> dest)
{
    vigra_precondition(src.shape() == dest.shape(),
        "tensorEigenRepresentation(): shape mismatch between input and output.");
    tensorEigenRepresentation(srcImageRange(src), destImage(dest));
}

/********************************************************/
/*                                                      */
/*                      tensorTrace                     */
/*                                                      */
/********************************************************/

/** \brief Calculate the trace of a 2x2 tensor.

    This function turns a 3-band image representing the tensor components
    t11, t12 (== t21 due to symmetry), t22 into the a 1-band image holding the 
    tensor trace t11 + t22.
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        tensorTrace(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> dest);
    }
    \endcode

    \deprecatedAPI{tensorTrace}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void tensorTrace(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                         DestIterator dul, DestAccessor dest);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void tensorTrace(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                         pair<DestIterator, DestAccessor> d);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/tensorutilities.hxx\>

    \code
    FVector3Image tensor(w,h);
    FImage trace(w,h);
    
    tensorTrace(srcImageRange(tensor), destImage(trace));
    \endcode

*/
doxygen_overloaded_function(template <...> void tensorTrace)

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
inline void
tensorTrace(triple<SrcIterator, SrcIterator, SrcAccessor> src,
            pair<DestIterator, DestAccessor> dest)
{
    tensorTrace(src.first, src.second, src.third, dest.first, dest.second);
}

template <class T1, class S1,
          class T2, class S2>
inline void
tensorTrace(MultiArrayView<2, T1, S1> const & src,
            MultiArrayView<2, T2, S2> dest)
{
    vigra_precondition(src.shape() == dest.shape(),
        "tensorTrace(): shape mismatch between input and output.");
    tensorTrace(srcImageRange(src), destImage(dest));
}

/********************************************************/
/*                                                      */
/*                  tensorToEdgeCorner                  */
/*                                                      */
/********************************************************/

/** \brief Decompose a symmetric 2x2 tensor into its edge and corner parts.

    This function turns a 3-band image representing the tensor components
    t11, t12 (== t21 due to symmetry), t22 into the a 2-band image holding 
    the tensor's edgeness (difference of the tensor's 
    eigenvalues) and orientation, and a 1-band image representing its corner part 
    (equal to the twice the small eigen value). The original tensor must be 
    positive definite and defined in a right-handed coordinate system (e.g.
    the tensor resulting from \ref boundaryTensor()).
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T21, class S21,
                  class T22, class S22>
        void
        tensorToEdgeCorner(MultiArrayView<2, T1, S1> const & src,
                           MultiArrayView<2, T21, S21> edge,
                           MultiArrayView<2, T22, S22> corner);
    }
    \endcode

    \deprecatedAPI{tensorToEdgeCorner}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator1, class DestAccessor1,
                  class DestIterator2, class DestAccessor2>
        void tensorToEdgeCorner(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                DestIterator1 edgeul, DestAccessor1 edge,
                                DestIterator2 cornerul, DestAccessor2 corner);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator1, class DestAccessor1,
                  class DestIterator2, class DestAccessor2>
        void tensorToEdgeCorner(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                                pair<DestIterator1, DestAccessor1> edge,
                                pair<DestIterator2, DestAccessor2> corner);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/tensorutilities.hxx\>

    \code
    FVector3Image tensor(w,h);
    FVector2Image edgePart(w,h);
    FImage cornerPart(w,h);
    
    tensorTrace(srcImageRange(tensor), destImage(edgePart), destImage(cornerPart));
    \endcode

*/
doxygen_overloaded_function(template <...> void tensorToEdgeCorner)

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
            TmpType d4 = (TmpType)hypot(d2, d3);
            
            edge.setComponent(d4, e, 0); // edgeness = difference of EVs
            if(d2 == 0.0 && d3 == 0.0)
            {
                edge.setComponent(0.0, e, 1); // orientation
            }
            else
            {
                edge.setComponent(0.5 * VIGRA_CSTD::atan2(d3, d2), e, 1); // orientation
            }
            corner.set(d1 - d4, c); // cornerness = 2 * small EV
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator1, class DestAccessor1,
          class DestIterator2, class DestAccessor2>
inline void
tensorToEdgeCorner(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator1, DestAccessor1> edge,
                   pair<DestIterator2, DestAccessor2> corner)
{
    tensorToEdgeCorner(src.first, src.second, src.third, 
                       edge.first, edge.second, corner.first, corner.second);
}

template <class T1, class S1,
          class T21, class S21,
          class T22, class S22>
inline void
tensorToEdgeCorner(MultiArrayView<2, T1, S1> const & src,
                   MultiArrayView<2, T21, S21> edge,
                   MultiArrayView<2, T22, S22> corner)
{
    vigra_precondition(src.shape() == edge.shape(),
        "tensorToEdgeCorner(): shape mismatch between input and output.");
    tensorToEdgeCorner(srcImageRange(src), 
                       destImage(edge), destImage(corner));
}

//@}

} // namespace vigra

#endif /* VIGRA_TENSORUTILITIES_HXX */
