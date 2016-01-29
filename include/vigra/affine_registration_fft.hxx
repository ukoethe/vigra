/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2014 by Benjamin Seppke                 */
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

#ifndef VIGRA_AFFINE_REGISTRATION_FFT_HXX
#define VIGRA_AFFINE_REGISTRATION_FFT_HXX

#include "mathutil.hxx"
#include "matrix.hxx"
#include "linear_solve.hxx"
#include "tinyvector.hxx"
#include "splineimageview.hxx"
#include "imagecontainer.hxx"
#include "multi_shape.hxx"
#include "affinegeometry.hxx"
#include "correlation.hxx"

#include <cmath>

namespace vigra {

/** \addtogroup Registration Image Registration
*/
//@{

/********************************************************/
/*                                                      */
/*               transformToPolarCoordinates            */
/*                                                      */
/********************************************************/

/** \brief Transforms a given image to its (image-centered) polar coordinates representation.

 This algorithm transforms a given image (by means of an spline image view) to its
 image-centered polar coordinates reprensentation. The sampling of the polar coordinate system
 is determined by the shape of the dest. image.

 <b> Declarations:</b>

 <b>\#include</b> \<vigra/affine_registration_fft.hxx\><br>
 Namespace: vigra

 pass 2D array views:
 \code
 namespace vigra {
     template <class SplineImage,
               class T1, class S1>
     void
     transformToPolarCoordinates(SplineImage const & src,
                                 MultiArrayView<2, T1, S1> dest);
 }
 \endcode

 \deprecatedAPI{estimateTranslation}
     pass \ref ImageIterators and \ref DataAccessors :
     \code
     namespace vigra {
         template <class SplineImage,
                   class DestIterator, class DestAccessor>
         void
         transformToPolarCoordinates(SplineImage const & src,
                                     DestIterator dul, DestIterator dlr, DestAccessor dest)
     }
     \endcode
     use argument objects in conjunction with \ref ArgumentObjectFactories :
     \code
         namespace vigra {
             template <class SplineImage,
                       class DestIterator, class DestAccessor>
             void
             transformToPolarCoordinates(SplineImage const & src,
                                 triple<DestIterator, DestIterator, DestAccessor> dest)
        }
     \endcode
 \deprecatedEnd
*/

template <class SplineImage,
          class DestIterator, class DestAccessor>
void
transformToPolarCoordinates(SplineImage const & src,
                            DestIterator d_ul, DestIterator d_lr, DestAccessor d_acc)
{
    typename DestIterator::difference_type d_shape = (d_lr - d_ul);

    int s_w = src.width(),
        s_h = src.height();

    int s_size = min(s_w, s_h);

    int d_w = d_shape.x,
        d_h = d_shape.y;

    double r_max = s_size / 2.0;

    DestIterator yd = d_ul;
    DestIterator xd = yd;

    for (int t_step = 0; t_step < d_h; t_step++, yd.y++)
    {
        xd = yd;
        for (int r_step = 0; r_step < d_w; r_step++, xd.x++)
        {
            double theta = 2.0 * M_PI * double(t_step) / double(d_h);
            double r = r_max * double(r_step) / double(d_w);
            double u = r * cos(theta) + r_max;
            double v = r * -sin(theta) + r_max;

            if (   u >= 0 && u < s_size
                && v >= 0 && v < s_size)
            {
                d_acc.set(src(u, v), xd);
            }
        }
    }
}

template <class SplineImage,
          class DestIterator, class DestAccessor>
inline void
transformToPolarCoordinates(SplineImage const & src,
                            vigra::triple<DestIterator, DestIterator, DestAccessor> dest)
{
    transformToPolarCoordinates(src, dest.first, dest.second, dest.third);
}

template <class SplineImage,
          class T1, class S1>
void
transformToPolarCoordinates(SplineImage const & src,
                            MultiArrayView<2, T1, S1> dest)
{
    transformToPolarCoordinates(src, srcImageRange(dest));
}




namespace detail
{
    template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
    void
    maximumFastNCC(SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,
                   DestIterator d_ul, DestIterator d_lr, DestAccessor d_acc,
                   TinyVector<int,2> & position,
                   double & correlation_coefficent)
    {
        typename DestIterator::difference_type s_shape = s_lr - s_ul;
        typename DestIterator::difference_type d_shape = d_lr - d_ul;

        MultiArray<2, float> src(s_shape.x, s_shape.y), dest(d_shape.x, d_shape.y), ncc(d_shape.x, d_shape.y);
        BasicImageView<float> src_view = makeBasicImageView(src);
        BasicImageView<float> dest_view = makeBasicImageView(dest);

        copyImage(srcIterRange(s_ul, s_lr, s_acc), destImage(src_view));
        copyImage(srcIterRange(d_ul, d_lr, d_acc), destImage(dest_view));

        fastNormalizedCrossCorrelation(dest, src, ncc);

        int max_x = 0;
        int max_y = 0;
        float val = 0.0;
        float max_val = -1.0;

        for (int y = 0; y < ncc.height()-s_shape.y; y++)
        {
            for (int x = 0; x < ncc.width()-s_shape.x; x++)
            {
                val = ncc(x+s_shape.x/2, y+s_shape.y/2);

                if (val > max_val)
                {
                    max_x = x;
                    max_y = y;
                    max_val = val;
                }
            }
        }

        position[0] = max_x;
        position[1] = max_y;

        correlation_coefficent = max_val;
    }

    template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
    inline void
    maximumFastNCC(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   triple<DestIterator, DestIterator, DestAccessor> dest,
                   TinyVector<int,2> & position,
                   double & correlation_coefficent)
    {
        maximumFastNCC(src.first, src.second, src.third,
                       dest.first, dest.second, dest.third,
                       position,
                       correlation_coefficent);
    }

    template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
    void fourierLogAbsSpectrumInPolarCoordinates(SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,
                                                 DestIterator d_ul, DestIterator d_lr, DestAccessor d_acc)
    {
        using namespace vigra;

        typename SrcIterator::difference_type shape = s_lr - s_ul;

        FFTWComplexImage fourier(shape);
        FImage           fft_mag(shape);

        fourierTransform(srcIterRange(s_ul, s_lr, s_acc), destImage(fourier));
        moveDCToCenter(srcImageRange(fourier, FFTWMagnitudeAccessor<>()), destImage(fft_mag));

        vigra::SplineImageView<2, double> spl(srcImageRange(fft_mag));

        transformToPolarCoordinates(spl,
                                    destIterRange(d_ul, d_lr, d_acc));

        transformImage(srcIterRange(d_ul,d_lr,d_acc), destIter(d_ul,d_acc), log(abs(functor::Arg1())));
    }

    template <class SrcIterator, class SrcAccessor,
              class DestIterator, class DestAccessor>
    void fourierLogAbsSpectrumInPolarCoordinates(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                                 triple<DestIterator, DestIterator, DestAccessor> dest)
    {
        fourierLogAbsSpectrumInPolarCoordinates(src.first, src.second, src.third, dest.first, dest.second, dest.third);
    }
} //namespace detail




/********************************************************/
/*                                                      */
/*                 estimateGlobalRotation               */
/*                                                      */
/********************************************************/

/** \brief Estimate the rotation between two images by means of a normalized cross correlation matching of the FFT spectra.

 This algorithm uses the fast normalized cross correlation to determine a global rotation
 between two images (from image2 to image1). To derive the rotation, the algorithm performs the following steps:
 <ol>
    <li>Transforming both images to the frequency domain using FFTW</li>
    <li>Create LogAbs PolarCoordinate representations for each spectrum.</li>
    <li>Determining the final Rotation by performing a fast normalized cross correlation
        based on the polar representations.</li>
 </ol>
 The images are cropped to the corresponding images center-squared before the estimation
 takes place.

 <b> Declarations:</b>

 <b>\#include</b> \<vigra/affine_registration_fft.hxx\><br>
 Namespace: vigra

 pass 2D array views:
 \code
 namespace vigra {
     template <class T1, class S1,
               class T2, class S2>
     void
     estimateGlobalRotation(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            Matrix<double> & affineMatrix,
                            double & correlation_coefficent);
 }
 \endcode

 \deprecatedAPI{estimateGlobalRotation}
     pass \ref ImageIterators and \ref DataAccessors :
     \code
     namespace vigra {
         template <class SrcIterator, class SrcAccessor,
                   class DestIterator, class DestAccessor>
         void
         estimateGlobalRotation(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                DestIterator dul, DestIterator dlr, DestAccessor dest,
                                Matrix<double> & affineMatrix,
                                double & correlation_coefficent)
     }
     \endcode
     use argument objects in conjunction with \ref ArgumentObjectFactories :
     \code
         namespace vigra {
             template <class SrcIterator, class SrcAccessor,
                       class DestIterator, class DestAccessor>
             void
             estimateGlobalRotation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                    triple<DestIterator, DestIterator, DestAccessor> dest,
                                    Matrix<double> & affineMatrix,
                                    double & correlation_coefficent)
        }
     \endcode
 \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void estimateGlobalRotation)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
estimateGlobalRotation(SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,
                        DestIterator d_ul, DestIterator d_lr, DestAccessor d_acc,
                       Matrix<double> & affineMatrix,
                       double & correlation_coefficient)
{
    //determine squared centers of both images without consuming any additional memory!
    typename SrcIterator::difference_type s_shape = s_lr - s_ul;
    Diff2D s2_shape(min(s_shape.x, s_shape.y),min(s_shape.x, s_shape.y));
    Diff2D s2_offset = (s_shape-s2_shape)/2;

    typename DestIterator::difference_type d_shape = d_lr - d_ul;
    Diff2D d2_shape(min(d_shape.x, d_shape.y),min(d_shape.x, d_shape.y));
    Diff2D d2_offset = (d_shape-d2_shape)/2;

    //Determine Shape for united polar coordinate representation
    Diff2D mean_shape = (s_shape + d_shape)/2;

    int size = min(mean_shape.x, mean_shape.y);
    if(size %2 == 0)
        size++;

    const int pc_w = size,
              pc_h = size*6+1;


    DImage s_polCoords(pc_w, pc_h/2),
           d_polCoords(pc_w, pc_h),
           ncc(pc_w, pc_h);

    detail::fourierLogAbsSpectrumInPolarCoordinates(srcIterRange(s_ul+s2_offset, s_ul+s2_offset+s2_shape, s_acc),
                                                    destImageRange(d_polCoords));
    copyImage(srcIterRange(d_polCoords.upperLeft(), d_polCoords.upperLeft() + vigra::Diff2D(pc_w, pc_h/2), d_polCoords.accessor()),
              destImage(s_polCoords));

    detail::fourierLogAbsSpectrumInPolarCoordinates(srcIterRange(d_ul+d2_offset, d_ul+d2_offset+d2_shape, d_acc),
                                                    destImageRange(d_polCoords));

    //Basic Cross correlation is assumed to be faster here, as there are only pc_h comparisons possible...
    normalizedCrossCorrelation(srcImageRange(d_polCoords), srcImageRange(s_polCoords), destImage(ncc));

    int max_idx = 0;
    double max_val = -1;

    const int x=pc_w/2;
    double val;

    //Only look at a stripe for the maximum angle of rotation
    //at the image center, at find the best fitting angle...
    for (int y=0; y<pc_h/2; y++)
    {
        val = ncc(x,y+pc_h/4);

        if (val > max_val)
        {
            max_idx = y;
            max_val = val;
        }
    }

    double theta = double(max_idx) / pc_h * 360.0;

    affineMatrix = rotationMatrix2DDegrees(theta, vigra::TinyVector<double,2>(s_shape.x/2.0, s_shape.y/2.0));
    correlation_coefficient = max_val;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
estimateGlobalRotation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                        triple<DestIterator, DestIterator, DestAccessor> dest,
                        Matrix<double> & affineMatrix,
                        double & correlation_coefficient)
{
    estimateGlobalRotation(src.first, src.second, src.third,
                            dest.first, dest.second, dest.third,
                           affineMatrix,
                           correlation_coefficient);
}

template <class T1, class S1,
          class T2, class S2>
inline void
estimateGlobalRotation(MultiArrayView<2, T1, S1> const & src,
                       MultiArrayView<2, T2, S2> dest,
                       Matrix<double> & affineMatrix,
                       double & correlation_coefficent)
{
    estimateGlobalRotation(srcImageRange(src),
                           destImageRange(dest),
                           affineMatrix,
                           correlation_coefficent);
}

/********************************************************/
/*                                                      */
/*               estimateGlobalTranslation              */
/*                                                      */
/********************************************************/

/** \brief Estimate the translation between two images by means of a normalized cross correlation matching.

 This algorithm uses the fast normalized cross correlation to determine a global translation
 between two images (from image2 to image1). To derive the translation, the algorithm consists of differents steps:
 <ol>
    <li>Separation of the second image<br/>
        The second image (the one, for which the translation shall be determined) is cut into five
        subregions: UpperLeft, UpperRight, Center, LowerLeft and LowerRight, each of 1/4 the size of
        the original image. Using a border > 0 results in (all) overlapping regions.</li>
    <li>Cross-Correlation of the subimages to the first image<br/>
        The subimages are normalized cross-correlated to the (complete) first image.
        The resulting maximum-likelihood translations and the correlation coefficients are stored for the next step.</li>
    <li>Determining the final Translation by voting<br/>
        Each correlation vector gets one vote at the beginning. For each equality of derived motion vectors,
        the voting to these vectors is incremented. If the maximum number of votes is larger than 1, the vector with the
        most votes is chosen. If the maximum number of votes is 1, the vector with the maximum likelihood is choosen.</li>
 </ol>
 <b> Declarations:</b>

 <b>\#include</b> \<vigra/affine_registration_fft.hxx\><br>
 Namespace: vigra

 pass 2D array views:
 \code
 namespace vigra {
     template <class T1, class S1,
               class T2, class S2>
     void
     estimateGlobalTranslation(MultiArrayView<2, T1, S1> const & src,
                               MultiArrayView<2, T2, S2> dest,
                               Matrix<double> & affineMatrix,
                               double & correlation_coefficent,
                               Diff2D border = Diff2D(0,0));
 }
 \endcode

 \deprecatedAPI{estimateGlobalTranslation}
     pass \ref ImageIterators and \ref DataAccessors :
     \code
     namespace vigra {
         template <class SrcIterator, class SrcAccessor,
                   class DestIterator, class DestAccessor>
         void
         estimateGlobalTranslation(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                   DestIterator dul, DestIterator dlr, DestAccessor dest,
                                   Matrix<double> & affineMatrix,
                                   double & correlation_coefficent,
                                   Diff2D border = Diff2D(0,0))
     }
     \endcode
     use argument objects in conjunction with \ref ArgumentObjectFactories :
     \code
         namespace vigra {
             template <class SrcIterator, class SrcAccessor,
                       class DestIterator, class DestAccessor>
             void
             estimateGlobalTranslation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                       triple<DestIterator, DestIterator, DestAccessor> dest,
                                       Matrix<double> & affineMatrix,
                                       double & correlation_coefficent,
                                       Diff2D border = Diff2D(0,0))
        }
     \endcode
 \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void estimateGlobalTranslation)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void estimateGlobalTranslation(SrcIterator    s_ul, SrcIterator  s_lr, SrcAccessor  s_acc,
                               DestIterator d_ul, DestIterator d_lr, DestAccessor d_acc,
                               Matrix<double> & affine_matrix,
                               double & correlation_coefficent,
                               Diff2D border = Diff2D(0,0))
{
    typename SrcIterator::difference_type s_shape = s_lr - s_ul;

    //determine matrix by using 5 quater-matches and a maximum likelihood decision:
    Diff2D q_shape = (s_shape - border - border)/2;
    if (q_shape.x % 2 == 0)        q_shape.x--;
    if (q_shape.y % 2 == 0)        q_shape.y--;

    Diff2D q_offsets[5];
    q_offsets[0] = border;
    q_offsets[1] = Diff2D(s_shape.x, 0)/2 + border;
    q_offsets[2] = s_shape/4;
    q_offsets[3] = Diff2D(0, s_shape.y)/2 + border;
    q_offsets[4] = s_shape/2 + border;

    TinyVector<int,2> translation_vectors[5];
    double translation_correlations[5] = {0.0,0.0,0.0,0.0,0.0};
    int translation_votes[5] = {1,1,1,1,1};

    int max_corr_idx=0;

    for (int q=0; q!=5; q++)
    {
        Diff2D offset = q_offsets[q];

        //we are searching a transformation from img2 ->  transformed image1, thus we switch dest and tmp
        detail::maximumFastNCC(srcIterRange(d_ul+offset, d_ul+offset+q_shape, d_acc),
                               srcIterRange(s_ul, s_lr, s_acc),
                               translation_vectors[q],
                               translation_correlations[q]);

        translation_vectors[q] = translation_vectors[q] - TinyVector<int,2>(offset);

        if(translation_correlations[q] > translation_correlations[max_corr_idx])
        {
            max_corr_idx=q;
        }

        for (int q_old=0; q_old!=q; q_old++)
        {
            if (translation_vectors[q] == translation_vectors[q_old])
            {
                translation_votes[q_old]++;
            }
        }
    }

    int max_votes_idx=0;

    for (int q=0; q!=5; q++)
    {
        if(translation_votes[q] > translation_votes[max_votes_idx])
        {
            max_votes_idx=q;
        }
    }

    int best_idx=0;
    if(translation_votes[max_votes_idx] > 1)
    {
        best_idx = max_votes_idx;
    }
    else
    {
        best_idx = max_corr_idx;
    }

    affine_matrix = translationMatrix2D(translation_vectors[best_idx]);
    correlation_coefficent = translation_correlations[best_idx];
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
estimateGlobalTranslation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          triple<DestIterator, DestIterator, DestAccessor> dest,
                          Matrix<double> & affineMatrix,
                          double & correlation_coefficent,
                          Diff2D border = Diff2D(0,0))
{
    estimateGlobalTranslation(src.first, src.second, src.third,
                              dest.first, dest.second, dest.third,
                              affineMatrix,
                              correlation_coefficent,
                              border);
}

template <class T1, class S1,
          class T2, class S2>
inline void
estimateGlobalTranslation(MultiArrayView<2, T1, S1> const & src,
                          MultiArrayView<2, T2, S2> dest,
                          Matrix<double> & affineMatrix,
                          double & correlation_coefficent,
                          Diff2D border = Diff2D(0,0))
{
    estimateGlobalTranslation(srcImageRange(src),
                              destImageRange(dest),
                              affineMatrix,
                              correlation_coefficent,
                              border);
}

/********************************************************/
/*                                                      */
/*              estimateGlobalRotationTranslation       */
/*                                                      */
/********************************************************/

/** \brief Estimate the (global) rotation and translation between two images by means a normalized cross correlation matching.

 This algorithm use the functions \ref estimateGlobalRotation() and
 \ref estimateGlobalTranslation() to estimate a matrix which describes
 the global rotation and translation from the second to the first image.

 <b> Declarations:</b>

 <b>\#include</b> \<vigra/affine_registration_fft.hxx\><br>
 Namespace: vigra

 pass 2D array views:
 \code
 namespace vigra {
     template <class T1, class S1,
               class T2, class S2>
     void
     estimateGlobalRotationTranslation(MultiArrayView<2, T1, S1> const & src,
                                       MultiArrayView<2, T2, S2> dest,
                                       Matrix<double> & affineMatrix,
                                       double & rotation_correlation,
                                       double & translation_correlation,
                                       Diff2D border = Diff2D(0,0));
 }
 \endcode

 \deprecatedAPI{estimateGlobalRotationTranslation}
     pass \ref ImageIterators and \ref DataAccessors :
     \code
     namespace vigra {
         template <class SrcIterator, class SrcAccessor,
                   class DestIterator, class DestAccessor>
         void
         estimateGlobalRotationTranslation(SrcIterator sul, SrcIterator slr, SrcAccessor src,
                                           DestIterator dul, DestIterator dlr, DestAccessor dest,
                                           Matrix<double> & affineMatrix,
                                           double & rotation_correlation,
                                           double & translation_correlation,
                                           Diff2D border = Diff2D(0,0))
     }
     \endcode
     use argument objects in conjunction with \ref ArgumentObjectFactories :
     \code
         namespace vigra {
             template <class SrcIterator, class SrcAccessor,
                       class DestIterator, class DestAccessor>
             void
             estimateGlobalRotationTranslation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                               triple<DestIterator, DestIterator, DestAccessor> dest,
                                               Matrix<double> & affineMatrix,
                                               double & rotation_correlation,
                                               double & translation_correlation,
                                               Diff2D border = Diff2D(0,0))
        }
     \endcode
 \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void estimateGlobalRotationTranslation)
template <class SrcIterator, class SrcAccessor,
       class DestIterator, class DestAccessor>
void
estimateGlobalRotationTranslation(SrcIterator s_ul, SrcIterator s_lr, SrcAccessor s_acc,
                                  DestIterator d_ul, DestIterator d_lr, DestAccessor d_acc,
                                  Matrix<double> & affineMatrix,
                                  double & rotation_correlation,
                                  double & translation_correlation,
                                  Diff2D border = Diff2D(0,0))
{
    typename DestIterator::difference_type d_shape = d_lr - d_ul;

    //First step: Estimate rotation from img2 -> img1.
    Matrix<double> rotation_matrix;
    estimateGlobalRotation(srcIterRange(s_ul+border, s_lr-border, s_acc),
                           srcIterRange(d_ul+border, d_lr-border, d_acc),
                           rotation_matrix,
                           rotation_correlation);

    //Second step: correct image according to the estimated rotation:
    FImage tmp(d_shape);
    SplineImageView<3, double> spl(srcIterRange(s_ul, s_lr, s_acc));
    affineWarpImage(spl, destImageRange(tmp), rotation_matrix);

    //Third step: find rotation between temp image (of step 2) and dest:
    Matrix<double> translation_matrix;
    estimateGlobalTranslation(srcImageRange(tmp),
                              srcIterRange(d_ul, d_lr, d_acc),
                              translation_matrix,
                              translation_correlation,
                              border);

    affineMatrix = rotation_matrix * translation_matrix;
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
estimateGlobalRotationTranslation(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                  triple<DestIterator, DestIterator, DestAccessor> dest,
                                  Matrix<double> & affineMatrix,
                                  double & rotation_correlation,
                                  double & translation_correlation,
                                  Diff2D border = Diff2D(0,0))
{
    estimateGlobalRotationTranslation(src.first, src.second, src.third,
                                      dest.first, dest.second, dest.third,
                                      affineMatrix,
                                      rotation_correlation,
                                      translation_correlation,
                                      border);
}

template <class T1, class S1,
          class T2, class S2>
inline void
estimateGlobalRotationTranslation(MultiArrayView<2, T1, S1> const & src,
                                  MultiArrayView<2, T2, S2> dest,
                                  Matrix<double> & affineMatrix,
                                  double & rotation_correlation,
                                  double & translation_correlation,
                                  Diff2D border = Diff2D(0,0))
{
    estimateGlobalRotationTranslation(srcImageRange(src),
                                      destImageRange(dest),
                                      affineMatrix,
                                      rotation_correlation,
                                      translation_correlation,
                                      border);
}

//@}

} // namespace vigra


#endif /* VIGRA_AFFINE_REGISTRATION_FFT_HXX */
