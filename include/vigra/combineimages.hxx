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
 
 
#ifndef VIGRA_COMBINEIMAGES_HXX
#define VIGRA_COMBINEIMAGES_HXX

#include "utilities.hxx"
#include "numerictraits.hxx"
#include "functortraits.hxx"
#include "multi_shape.hxx"

#include <cmath>

namespace vigra {

/** \addtogroup CombineAlgo Algorithms to Combine Images

    Apply functor to calculate a pixelwise transformation depending on multiple images.
    Note that the binary functors of the STL may be used with these functions.
*/
//@{

/********************************************************/
/*                                                      */
/*                    combine...Lines                   */
/*                                                      */
/********************************************************/

template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class DestIterator, class DestAccessor, class Functor>
void
combineTwoLines(SrcIterator1 s1, 
                SrcIterator1 s1end, SrcAccessor1 src1,
                SrcIterator2 s2, SrcAccessor2 src2,
                DestIterator d, DestAccessor dest,
                Functor const & f)
{
    for(; s1 != s1end; ++s1, ++s2, ++d)
        dest.set(f(src1(s1), src2(s2)), d);
}

template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class MaskIterator, class MaskAccessor, 
          class DestIterator, class DestAccessor, class Functor>
void
combineTwoLinesIf(SrcIterator1 s1, 
                  SrcIterator1 s1end, SrcAccessor1 src1,
                  SrcIterator2 s2, SrcAccessor2 src2,
                  MaskIterator m, MaskAccessor mask,
                  DestIterator d, DestAccessor dest,
                  Functor const & f)
{
    for(; s1 != s1end; ++s1, ++s2, ++m, ++d)
        if(mask(m))
            dest.set(f(src1(s1), src2(s2)), d);
}

template <class SrcIterator1, class SrcAccessor1,
          class SrcIterator2, class SrcAccessor2,
          class SrcIterator3, class SrcAccessor3,
          class DestIterator, class DestAccessor, class Functor>
void
combineThreeLines(SrcIterator1 s1, 
                  SrcIterator1 s1end, SrcAccessor1 src1,
                  SrcIterator2 s2, SrcAccessor2 src2,
                  SrcIterator3 s3, SrcAccessor3 src3,
                  DestIterator d, DestAccessor dest,
                  Functor const & f)
{
    for(; s1 != s1end; ++s1, ++s2, ++s3, ++d)
        dest.set(f(src1(s1), src2(s2), src3(s3)), d);
}

/********************************************************/
/*                                                      */
/*                    combineTwoImages                  */
/*                                                      */
/********************************************************/

/** \brief Combine two source images into destination image.

    After the introduction of arithmetic and algebraic \ref MultiMathModule "array expressions",
    this function is rarely needed. Moreover, \ref combineTwoMultiArrays() provides the 
    same functionality for arbitrary dimensional arrays.

    The transformation given by the functor is applied to the source 
    pixels and the result written into the corresponding destination pixel.
    This is typically used for operations like add and subtract.
    Note that the binary functors of the STL can be used in addition to
    the functors specifically defined in \ref CombineFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T11, class S11,
                  class T12, class S12,
                  class T2, class S2,
                  class Functor>
        void
        combineTwoImages(MultiArrayView<2, T11, S11> const & src1,
                         MultiArrayView<2, T12, S12> const & src2,
                         MultiArrayView<2, T2, S2> dest,
                         Functor const & f);
    }
    \endcode
    
    \deprecatedAPI{combineTwoImages}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class DestImageIterator, class DestAccessor,
              class Functor>
        void
        combineTwoImages(SrcImageIterator1 src1_upperleft, 
                 SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                 SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                 DestImageIterator dest_upperleft, DestAccessor da,
                 Functor const & f)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class DestImageIterator, class DestAccessor,
              class Functor>
        void
        combineTwoImages(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                 pair<SrcImageIterator2, SrcAccessor2> src2,
                 pair<DestImageIterator, DestAccessor> dest,
                 Functor const & f)
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/combineimages.hxx\><br>
    Namespace: vigra

    \code
    #include <functional>     // for plus
    MultiArray<2, float> src1(width, height), src2(width, height),
                         dest(width, height);
    ... // fill source images 
    
    combineTwoImages(src1, src2, dest,
                     std::plus<float>());
    \endcode

    \deprecatedUsage{combineTwoImages}
    \code
    #include <functional>     // for plus
    FImage src1(width, height), src2(width, height),
           dest(width, height);
    ... // fill source images 
    
    vigra::combineTwoImages(
                srcIterRange(src1.upperLeft(), src1.lowerRight()), 
                srcIter(src2.upperLeft()), 
                destIter(dest.upperLeft()),  
                std::plus<float>());
    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator1 src1_upperleft, src1_lowerright;
    SrcImageIterator2 src2_upperleft;
    DestImageIterator dest_upperleft;
    SrcImageIterator1::row_iterator sx1 = src1_upperleft.rowIterator();
    SrcImageIterator2::row_iterator sx2 = src2_upperleft.rowIterator();
    DestImageIterator::row_iterator dx = dest_upperleft.rowIterator();
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;

    dest_accessor.set(
          functor(src1_accessor(sx1), src2_accessor(sx2)), 
          dx);

    \endcode
    \deprecatedEnd
    
    \see TransformFunctor, MultiMathModule, \ref FunctorExpressions
*/
doxygen_overloaded_function(template <...> void combineTwoImages)

template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class DestImageIterator, class DestAccessor,
          class Functor>
void
combineTwoImages(SrcImageIterator1 src1_upperleft, 
                 SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                 SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                 DestImageIterator dest_upperleft, DestAccessor da,
                 Functor const & f)
{
    int w = src1_lowerright.x - src1_upperleft.x;
    
    for(; src1_upperleft.y < src1_lowerright.y; 
            ++src1_upperleft.y, ++src2_upperleft.y, ++dest_upperleft.y)
    {
        combineTwoLines(src1_upperleft.rowIterator(), 
                        src1_upperleft.rowIterator() + w, sa1, 
                        src2_upperleft.rowIterator(), sa2, 
                        dest_upperleft.rowIterator(), da, f);
    }
}
    
template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline
void
combineTwoImages(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                 pair<SrcImageIterator2, SrcAccessor2> src2,
                 pair<DestImageIterator, DestAccessor> dest,
                 Functor const & f)
{
    combineTwoImages(src1.first, src1.second, src1.third, 
                     src2.first, src2.second, 
                     dest.first, dest.second, f);
}

template <class T11, class S11,
          class T12, class S12,
          class T2, class S2,
          class Functor>
inline void
combineTwoImages(MultiArrayView<2, T11, S11> const & src1,
                 MultiArrayView<2, T12, S12> const & src2,
                 MultiArrayView<2, T2, S2> dest,
                 Functor const & f)
{
    vigra_precondition(src1.shape() == src2.shape() && src1.shape() == dest.shape(),
        "combineTwoImages(): shape mismatch between inputs and/or output.");
    combineTwoImages(srcImageRange(src1), 
                     srcImage(src2), 
                     destImage(dest), f);
}

/********************************************************/
/*                                                      */
/*                  combineTwoImagesIf                  */
/*                                                      */
/********************************************************/

/** \brief Combine ROI of two source images into destination image.

    The transformation given by the functor is applied to all source 
    pixels in the ROI (i.e. whenever the corresponding value of the mask array
    is non-zero)
    and the result written into the corresponding destination pixel.
    This is typically used for operations like add and subtract.
    Note that the binary functors of the STL can be used in addition to
    the functors specifically defined in \ref CombineFunctor.
    Creation of new functors is easiest by using \ref FunctorExpressions.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T11, class S11,
                  class T12, class S12,
                  class TM, class SM,
                  class T2, class S2,
                  class Functor>
        void
        combineTwoImagesIf(MultiArrayView<2, T11, S11> const & src1,
                           MultiArrayView<2, T12, S12> const & src2,
                           MultiArrayView<2, TM, SM> const & mask,
                           MultiArrayView<2, T2, S2> dest,
                           Functor const & f);
    }
    \endcode
    
    \deprecatedAPI{combineTwoImagesIf}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class MaskImageIterator, class MaskAccessor,
              class DestImageIterator, clas DestAccessor,
              class Functor>
        void
        combineTwoImagesIf(SrcImageIterator1 src1_upperleft, 
                   SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                   MaskImageIterator mask_upperleft, MaskAccessor ma,
                   DestImageIterator dest_upperleft, DestAccessor da,
                   Functor const & f)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class MaskImageIterator, class MaskAccessor,
              class DestImageIterator, clas DestAccessor,
              class Functor>
        void
        combineTwoImagesIf(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                   pair<SrcImageIterator2, SrcAccessor2> src2,
                   pair<MaskImageIterator, MaskAccessor> mask,
                   pair<DestImageIterator, DestAccessor> dest,
                   Functor const & f)
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/combineimages.hxx\><br>
    Namespace: vigra

    \code
    #include <functional>     // for plus
    MultiArray<2, float> src1(width, height), src2(width, height), mask(width, height),
                         dest(width, height);
    ... // fill source and mask images 

    combineTwoImagesIf(src1, src2, mask, dest,
                       std::plus<float>());
    \endcode

    \deprecatedUsage{combineTwoImagesIf}
    \code
    #include <functional>     // for plus
    FImage src1(width, height), src2(width, height), mask(width, height),
           dest(width, height);
    ... // fill source and mask images 

    vigra::combineTwoImagesIf(
                srcIterRange(src1.upperLeft(), src1.lowerRight()), 
                srcIter(src2.upperLeft()), 
                maskIter(mask.upperLeft()), 
                destIter(dest.upperLeft()),
                std::plus<SrcValueType>());
    
    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator1 src1_upperleft, src1_lowerright;
    SrcImageIterator2 src2_upperleft;
    MaskImageIterator mask_upperleft;
    DestImageIterator dest_upperleft;
    SrcImageIterator1::row_iterator sx1 = src1_upperleft.rowIterator();
    SrcImageIterator2::row_iterator sx2 = src2_upperleft.rowIterator();
    MaskImageIterator::row_iterator mx = mask_upperleft.rowIterator();
    DestImageIterator::row_iterator dx = dest_upperleft.rowIterator();
    
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    MaskAccessor mask_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;
    
    if(mask_accessor(mx))
       dest_accessor.set(
          functor(src1_accessor(sx1), src2_accessor(sx2)), 
          dx);

    \endcode
    \deprecatedEnd
    
    \see TransformFunctor, MultiMathModule, \ref FunctorExpressions
*/
doxygen_overloaded_function(template <...> void combineTwoImagesIf)

template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
void
combineTwoImagesIf(SrcImageIterator1 src1_upperleft, 
                   SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                   MaskImageIterator mask_upperleft, MaskAccessor ma,
                   DestImageIterator dest_upperleft, DestAccessor da,
                   Functor const & f)
{
    int w = src1_lowerright.x - src1_upperleft.x;
    
    for(; src1_upperleft.y < src1_lowerright.y;
          ++src1_upperleft.y, ++src2_upperleft.y, 
          ++dest_upperleft.y, ++mask_upperleft.y)
    {
        combineTwoLinesIf(src1_upperleft.rowIterator(), 
                          src1_upperleft.rowIterator() + w, sa1, 
                          src2_upperleft.rowIterator(), sa2, 
                          mask_upperleft.rowIterator(), ma, 
                          dest_upperleft.rowIterator(), da, f);
    }
}
    
template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class MaskImageIterator, class MaskAccessor,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline
void
combineTwoImagesIf(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                   pair<SrcImageIterator2, SrcAccessor2> src2,
                   pair<MaskImageIterator, MaskAccessor> mask,
                   pair<DestImageIterator, DestAccessor> dest,
                   Functor const & f)
{
    combineTwoImagesIf(src1.first, src1.second, src1.third, 
                       src2.first, src2.second, 
                       mask.first, mask.second, 
                       dest.first, dest.second, f);
}
    
template <class T11, class S11,
          class T12, class S12,
          class TM, class SM,
          class T2, class S2,
          class Functor>
inline void
combineTwoImagesIf(MultiArrayView<2, T11, S11> const & src1,
                   MultiArrayView<2, T12, S12> const & src2,
                   MultiArrayView<2, TM, SM> const & mask,
                   MultiArrayView<2, T2, S2> dest,
                   Functor const & f)
{
    vigra_precondition(src1.shape() == src2.shape() && src1.shape() == mask.shape() && src1.shape() == dest.shape(),
        "combineTwoImagesIf(): shape mismatch between inputs and/or output.");
    combineTwoImagesIf(srcImageRange(src1), 
                       srcImage(src2), 
                       maskImage(mask), 
                       destImage(dest), f);
}

/********************************************************/
/*                                                      */
/*                  combineThreeImages                  */
/*                                                      */
/********************************************************/

/** \brief Combine three source images into destination image.

    After the introduction of arithmetic and algebraic \ref MultiMathModule "array expressions",
    this function is rarely needed. Moreover, \ref combineThreeMultiArrays() provides the 
    same functionality for arbitrary dimensional arrays.

    The transformation given by the functor is applied to the source 
    pixels and the result written into the corresponding destination pixel.
    Creation of new functors is easiest by using \ref FunctorExpressions.
    
    <b> Declarations:</b>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T11, class S11,
                  class T12, class S12,
                  class T13, class S13,
                  class T2, class S2,
                  class Functor>
        void
        combineThreeImages(MultiArrayView<2, T11, S11> const & src1,
                           MultiArrayView<2, T12, S12> const & src2,
                           MultiArrayView<2, T13, S13> const & src3,
                           MultiArrayView<2, T2, S2> dest,
                           Functor const & f);
    }
    \endcode
    
    \deprecatedAPI{combineThreeImages}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class SrcImageIterator3, class SrcAccessor3,
              class DestImageIterator, class DestAccessor,
              class Functor>
        void
        combineThreeImages(SrcImageIterator1 src1_upperleft, 
                   SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                   SrcImageIterator3 src2_upperleft, SrcAccessor3 sa3,
                   DestImageIterator dest_upperleft, DestAccessor da,
                   Functor const & f)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator1, class SrcAccessor1,
              class SrcImageIterator2, class SrcAccessor2,
              class SrcImageIterator3, class SrcAccessor3,
              class DestImageIterator, class DestAccessor,
              class Functor>
        void
        combineThreeImages(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                 pair<SrcImageIterator2, SrcAccessor2> src2,
                 pair<SrcImageIterator3, SrcAccessor3> src3,
                 pair<DestImageIterator, DestAccessor> dest,
                 Functor const & f)
    }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/combineimages.hxx\><br>
    Namespace: vigra

    \code
    #include <vigra/functorexpression.hxx>
    
    MultiArray<2, float> src1(width, height), src2(width, height), src3(width, height),
                         dest(width, height);
    ... // fill source images 
    
    using namespace vigra::functor; // activate VIGRA's lambda library
    
    combineThreeImages(src1, src2, src3, dest,
                       Arg1()*exp(-abs(Arg2()-Arg3())));
    \endcode

    \deprecatedUsage{combineThreeImages}
    \code
    FImage src1(width, height), src2(width, height), src3(width, height),
           dest(width, height);
    ... // fill source images 
    
    vigra::combineThreeImages(
                srcIterRange(src1.upperLeft(), src1.lowerRight()), 
                srcIter(src2.upperLeft()), 
                srcIter(src3.upperLeft()), 
                destIter(dest.upperLeft()),  
                SomeThreeArgumentFunctor());
    
    \endcode
    <b> Required Interface:</b>
    \code
    SrcImageIterator1 src1_upperleft, src1_lowerright;
    SrcImageIterator2 src2_upperleft;
    SrcImageIterator3 src3_upperleft;
    DestImageIterator dest_upperleft;
    SrcImageIterator1::row_iterator sx1 = src1_upperleft.rowIterator();
    SrcImageIterator2::row_iterator sx2 = src2_upperleft.rowIterator();
    SrcImageIterator3::row_iterator sx3 = src3_upperleft.rowIterator();
    DestImageIterator::row_iterator dx = dest_upperleft.rowIterator();
    
    SrcAccessor1 src1_accessor;
    SrcAccessor2 src2_accessor;
    SrcAccessor3 src3_accessor;
    DestAccessor dest_accessor;
    
    Functor functor;

    dest_accessor.set(
          functor(src1_accessor(sx1), 
                  src2_accessor(sx2), 
                  src3_accessor(sx3)), 
          dx);

    \endcode
    \deprecatedEnd
    
    \see TransformFunctor, MultiMathModule, \ref FunctorExpressions
*/
doxygen_overloaded_function(template <...> void combineThreeImages)

template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class SrcImageIterator3, class SrcAccessor3,
          class DestImageIterator, class DestAccessor,
          class Functor>
void
combineThreeImages(SrcImageIterator1 src1_upperleft, 
                   SrcImageIterator1 src1_lowerright, SrcAccessor1 sa1,
                   SrcImageIterator2 src2_upperleft, SrcAccessor2 sa2,
                   SrcImageIterator3 src3_upperleft, SrcAccessor3 sa3,
                   DestImageIterator dest_upperleft, DestAccessor da,
                   Functor const & f)
{
    int w = src1_lowerright.x - src1_upperleft.x;
    
    for(; src1_upperleft.y < src1_lowerright.y;
        ++src1_upperleft.y, ++src2_upperleft.y, ++src3_upperleft.y, 
        ++dest_upperleft.y)
    {
        combineThreeLines(src1_upperleft.rowIterator(), 
                          src1_upperleft.rowIterator() + w, sa1, 
                          src2_upperleft.rowIterator(), sa2, 
                          src3_upperleft.rowIterator(), sa3, 
                          dest_upperleft.rowIterator(), da, f);
    }
}
    
template <class SrcImageIterator1, class SrcAccessor1,
          class SrcImageIterator2, class SrcAccessor2,
          class SrcImageIterator3, class SrcAccessor3,
          class DestImageIterator, class DestAccessor,
          class Functor>
inline
void
combineThreeImages(triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                   pair<SrcImageIterator2, SrcAccessor2> src2,
                   pair<SrcImageIterator3, SrcAccessor3> src3,
                   pair<DestImageIterator, DestAccessor> dest,
                   Functor const & f)
{
    combineThreeImages(src1.first, src1.second, src1.third, 
                       src2.first, src2.second, 
                       src3.first, src3.second, 
                       dest.first, dest.second, f);
}

template <class T11, class S11,
          class T12, class S12,
          class T13, class S13,
          class T2, class S2,
          class Functor>
inline void
combineThreeImages(MultiArrayView<2, T11, S11> const & src1,
                   MultiArrayView<2, T12, S12> const & src2,
                   MultiArrayView<2, T13, S13> const & src3,
                   MultiArrayView<2, T2, S2> dest,
                   Functor const & f)
{
    vigra_precondition(src1.shape() == src2.shape() && src1.shape() == src3.shape() && src1.shape() == dest.shape(),
        "combineThreeImages(): shape mismatch between inputs and/or output.");
    combineThreeImages(srcImageRange(src1), 
                       srcImage(src2), 
                       srcImage(src3), 
                       destImage(dest), f);
}

    
//@}

/** \addtogroup CombineFunctor Functors to Combine Images

    Common functors with several arguments
*/
//@{

/********************************************************/
/*                                                      */
/*                    MagnitudeFunctor                  */
/*                                                      */
/********************************************************/

/** Calculate the magnitude from two arguments.
    Can be used in conjunction with \ref gradientBasedTransform().
    
    If the gradient is represented by a vector-valued image instead of 
    a pair of scalar images, use \ref vigra::VectorNormFunctor.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isBinaryFunctor</tt> are true (<tt>VigraTrueType</tt>)    
*/
template <class ValueType>
class MagnitudeFunctor
{
  public:
        /** the functor's first argument type
        */
    typedef ValueType first_argument_type;
    
        /** the functor's second argument type
        */
    typedef ValueType second_argument_type;
    
        /** the functor's result type
        */
    typedef typename SquareRootTraits<typename NormTraits<ValueType>::SquaredNormType>::SquareRootResult result_type;
    
        /** \deprecated use first_argument_type, second_argument_type, result_type
        */
    typedef ValueType value_type;
    
        /** calculate transform '<TT>sqrt(squaredNorm(v1) + squaredNorm(v2))</TT>'. 
            
        */
    result_type operator()(first_argument_type const & v1, second_argument_type const & v2) const
    {
        return VIGRA_CSTD::sqrt(squaredNorm(v1) + squaredNorm(v2));
    }
};

template <class T>
class FunctorTraits<MagnitudeFunctor<T> >
: public FunctorTraitsBase<MagnitudeFunctor<T> >
{
public:
    typedef VigraTrueType isBinaryFunctor;
};

/********************************************************/
/*                                                      */
/*             RGBGradientMagnitudeFunctor              */
/*                                                      */
/********************************************************/


/** Calculate the gradient magnitude from RGB arguments.
    Can be used in conjunction with \ref gradientBasedTransform().

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isBinaryFunctor</tt> are true (<tt>VigraTrueType</tt>)    
*/
template <class ValueType>
class RGBGradientMagnitudeFunctor
{
  public:
        /** the functor's first argument type
        */
    typedef RGBValue<ValueType> first_argument_type;
    
        /** the functor's second argument type
        */
    typedef RGBValue<ValueType> second_argument_type;
    
        /** the functor's result type
        */
    typedef typename NumericTraits<ValueType>::RealPromote result_type;
    
        /** \deprecated use first_argument_type, second_argument_type, result_type
        */
    typedef ValueType value_type;
    
        /** Calculate the gradient magnitude form given RGB components.
            The function returns 
            
            \f[ \sqrt{|\nabla red|^2 + |\nabla green|^2 + |\nabla blue|^2}
            \f]
            
            where \f$|\nabla red|^2\f$ is defined by <TT>gx.red()*gx.red() + gy.red()*gy.red()</TT>.
            
            <TT>ValueType</TT> (the RGB's component type) must support addition, multiplication, 
            abd <TT>sqrt()</TT>.
        */
    result_type 
    operator()(first_argument_type const & gx, second_argument_type const & gy) const
    {
        return VIGRA_CSTD::sqrt(gx.red()*gx.red() + gx.green()*gx.green() +
                    gx.blue()*gx.blue() + gy.red()*gy.red() + 
                    gy.green()*gy.green() + gy.blue()*gy.blue());
    }
};

template <class T>
class FunctorTraits<RGBGradientMagnitudeFunctor<T> >
: public FunctorTraitsBase<RGBGradientMagnitudeFunctor<T> >
{
public:
    typedef VigraTrueType isBinaryFunctor;
};

//@}

} // namespace vigra

#endif // VIGRA_COMBINEIMAGES_HXX
