// -*- indent-tabs-mode: nil -*-
#ifndef VIGRA_GABORFILTER_HXX
#define VIGRA_GABORFILTER_HXX

#include "vigra/imagecontainer.hxx"
#include "vigra/config.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/combineimages.hxx"
#include "vigra/utilities.hxx"
#include "vigra/fftw.hxx"

#include <functional>
#include <vector>
#include <cmath>

namespace vigra {

/** \addtogroup GaborFilter Gabor Filter
    Functions to create or apply gabor filter (latter based on FFTW).
*/
//@{

/********************************************************/
/*                                                      */
/*                   createGaborFilter                  */
/*                                                      */
/********************************************************/

/** \brief Create a gabor filter in frequency space.

    The orientation is given in radians, the other parameters are the
    center frequency (for example 0.375 or smaller) and the two
    angular and radial sigmas of the gabor filter. (See \ref
    angularGaborSigma() for an explanation of possible values.)

    The energy of the filter is explicitely normalized to 1.0.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class DestImageIterator, class DestAccessor>
        void createGaborFilter(DestImageIterator destUpperLeft,
                               DestImageIterator destLowerRight,
                               DestAccessor da,
                               double orientation, double centerFrequency,
                               double angularSigma, double radialSigma)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class DestImageIterator, class DestAccessor>
        inline
        void createGaborFilter(triple<DestImageIterator,
                                      DestImageIterator,
                                      DestAccessor> dest,
                               double orientation, double centerFrequency,
                               double angularSigma, double radialSigma)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="gaborfilter_8hxx-source.html">vigra/gaborfilter.hxx</a>"<br>

    Namespace: vigra

    \code
    vigra::FImage gabor(w,h);

    vigra::createGaborFilter(destImageRange(gabor), orient, freq,
                             angularGaborSigma(directionCount, freq)
                             radialGaborSigma(freq));
    \endcode
*/

template <class DestImageIterator, class DestAccessor>
void createGaborFilter(DestImageIterator destUpperLeft,
                       DestImageIterator destLowerRight, DestAccessor da,
                       double orientation, double centerFrequency,
                       double angularSigma, double radialSigma)
{
    int w= destLowerRight.x - destUpperLeft.x;
    int h= destLowerRight.y - destUpperLeft.y;

    double squaredSum = 0.0;
    double cosTheta= VIGRA_CSTD::cos(orientation);
    double sinTheta= VIGRA_CSTD::sin(orientation);

    double radialSigma2 = radialSigma*radialSigma;
    double angularSigma2 = angularSigma*angularSigma;

    double wscale = w % 1 ?
                    1.0f / (w-1) :
                    1.0f / w;
    double hscale = h % 1 ?
                    1.0f / (h-1) :
                    1.0f / h;

    int dcX= (w+1)/2, dcY= (h+1)/2;

    double u, v;
    for ( int y=0; y<h; y++, destUpperLeft.y++ )
    {
        typename DestImageIterator::row_iterator dix = destUpperLeft.rowIterator();

        v = hscale * ((h - (y - dcY))%h - dcY);
        for ( int x=0; x<w; x++, dix++ )
        {
            u= wscale*((x - dcX + w)%w - dcX);

            double uu = cosTheta*u + sinTheta*v - centerFrequency;
            double vv = -sinTheta*u + cosTheta*v;
            double gabor;

            gabor = VIGRA_CSTD::exp(-0.5*(uu*uu / radialSigma2 + vv*vv / angularSigma2));
            squaredSum += gabor * gabor;
            da.set( gabor, dix );
        }
    }
    destUpperLeft.y -= h;

    // clear out DC value and remove it from the squared sum
    double dcValue = da(destUpperLeft);
    squaredSum -= dcValue * dcValue;
    da.set( 0.0, destUpperLeft );

    // normalize energy to one
    double factor = VIGRA_CSTD::sqrt(squaredSum);
    for ( int y=0; y<h; y++, destUpperLeft.y++ )
    {
        typename DestImageIterator::row_iterator dix = destUpperLeft.rowIterator();

        for ( int x=0; x<w; x++, dix++ )
        {
            da.set( da(dix) / factor, dix );
        }
    }
}

template <class DestImageIterator, class DestAccessor>
inline
void createGaborFilter(triple<DestImageIterator, DestImageIterator,
                              DestAccessor> dest,
                              double orientation, double centerFrequency,
                              double angularSigma, double radialSigma)
{
    createGaborFilter(dest.first, dest.second, dest.third,
                      orientation, centerFrequency,
                      angularSigma, radialSigma);
}

/********************************************************/
/*                                                      */
/*                   radialGaborSigma                   */
/*                                                      */
/********************************************************/

/** \brief Calculate sensible radial sigma for given parameters.

    For a brief introduction what is meant with "sensible" sigmas, see
    \ref angularGaborSigma().

    <b> Declaration:</b>

    \code
    namespace vigra {
        inline
        double radialGaborSigma(double centerFrequency)
    }
    \endcode
 */

inline double radialGaborSigma(double centerFrequency)
{
    static double sfactor = 3.0 * VIGRA_CSTD::sqrt(VIGRA_CSTD::log(4.0));
    return centerFrequency / sfactor;
}

/********************************************************/
/*                                                      */
/*                   angularGaborSigma                  */
/*                                                      */
/********************************************************/

/** \brief Calculate sensible angular sigma for given parameters.

    "Sensible" means: If you use a range of gabor filters for feature
    detection, you are interested in minimal redundance. This is hard
    to define but one possible try is to arrange the filters in
    frequency space, so that the half-peak-magnitude ellipses touch
    each other.

    To do so, you must know the number of directions (first parameter
    for the angular sigma function) and the center frequency of the
    filter you want to calculate the sigmas for.

    The exact formulas are:
    \code
    sigma_radial= 1/sqrt(ln(4)) * centerFrequency/3
    \endcode

    \code
    sigma_angular= 1/sqrt(ln(4)) * tan(pi/(directions*2))
                   * sqrt(8/9) * centerFrequency
    \endcode

    <b> Declaration:</b>

    \code
    namespace vigra {
        inline
        double angularGaborSigma(int directionCount, double centerFrequency)
    }
    \endcode
 */

inline double angularGaborSigma(int directionCount, double centerFrequency)
{
    return VIGRA_CSTD::tan(M_PI/directionCount/2.0) * centerFrequency
        * VIGRA_CSTD::sqrt(8.0 / (9 * VIGRA_CSTD::log(4.0)));
}

/********************************************************/
/*                                                      */
/*                   GaborFilterFamily                  */
/*                                                      */
/********************************************************/

/** \brief Family of gabor filters of different scale and direction.

    A GaborFilterFamily can be used to quickly create a whole family
    of gabor filters in frequency space. Especially useful in
    conjunction with \ref applyFourierFilterFamily, since it's derived
    from \ref ImageArray.

    The filter parameters are chosen to make the center frequencies
    decrease in octaves with increasing scale indices, and to make the
    half-peak-magnitude ellipses touch each other to somewhat reduce
    redundancy in the filter answers. This is done by using \ref
    angularGaborSigma() and \ref radialGaborSigma(), you'll find more
    information there.

    The template parameter ImageType should be a scalar image type suitable for filling in

    <b>\#include</b> "<a href="gaborfilter_8hxx-source.html">vigra/gaborfilter.hxx</a>"

    Namespace: vigra
*/
template <class ImageType>
class GaborFilterFamily : public ImageArray<ImageType>
{
    typedef ImageArray<ImageType> ParentClass;
    int scaleCount_, directionCount_;
    double maxCenterFrequency_;

protected:
    void initFilters()
    {
        for(int direction= 0; direction<directionCount_; direction++)
            for(int scale= 0; scale<scaleCount_; scale++)
            {
                double angle = direction * M_PI / directionCount();
                double centerFrequency = maxCenterFrequency_ / VIGRA_CSTD::pow(2,scale);
                createGaborFilter(destImageRange(images_[filterIndex(direction, scale)]),
                                  angle, centerFrequency,
                                  angularGaborSigma(directionCount(), centerFrequency),
                                  radialGaborSigma(centerFrequency));
            }
    }

public:
    enum { stdFilterSize= 128, stdDirectionCount= 6, stdScaleCount= 4 };

        /** Constructs a family of gabor filters in frequency
            space. The filters will be calculated on construction, so
            it makes sense to provide good parameters right now
            although they can be changed later, too. If you leave them
            out, the defaults are a \ref directionCount of 6, a \ref
            scaleCount of 4 and a \ref maxCenterFrequency of
            3/8(=0.375).
        */
    GaborFilterFamily(const Diff2D & size,
                      int directionCount = stdDirectionCount, int scaleCount = stdScaleCount,
                      double maxCenterFrequency = 3.0/8.0)
        : ParentClass(directionCount*scaleCount, size),
          scaleCount_(scaleCount),
          directionCount_(directionCount),
          maxCenterFrequency_(maxCenterFrequency)
    {
        initFilters();
    }

        /** Convenience variant of the above constructor taking width
            and height separately. Also, this one serves as default
            constructor constructing 128x128 pixel filters.
         */
    GaborFilterFamily(int width= stdFilterSize, int height= -1,
                      int directionCount = stdDirectionCount, int scaleCount = stdScaleCount,
                      double maxCenterFrequency = 3.0/8.0)
        : ParentClass(directionCount*scaleCount, Diff2D(width, height > 0 ? height : width)),
          scaleCount_(scaleCount),
          directionCount_(directionCount),
          maxCenterFrequency_(maxCenterFrequency)
    {
        initFilters();
    }

        /** Return the index of the filter with the given direction and
            scale in this ImageArray. direction must in the range
            0..directionCount()-1 and scale in the range
            0..rangeCount()-1. This is useful for example if you used
            \ref applyFourierFilterFamily() and got a resulting
            ImageArray which still has the same order of images, but no
            \ref getFilter() method anymore.
         */
    int filterIndex(int direction, int scale) const
    {
        return scale*directionCount()+direction;
    }

        /** Return the filter with the given direction and
            scale. direction must in the range 0..directionCount()-1
            and scale in the range 0..rangeCount()-1.
            <tt>filters.getFilter(direction, scale)</tt> is the same as
            <tt>filters[filterIndex(direction, scale)]</tt>.
         */
    ImageType const & getFilter(int direction, int scale) const
    {
        return images_[filterIndex(direction, scale)];
    }

        /** Resize all filters (causing their recalculation).
         */
    virtual void resizeImages(const Diff2D &newSize)
    {
        ParentClass::resizeImages(newSize);
        initFilters();
    }

        /** Query the number of filter scales available.
         */
    int scaleCount() const
        { return scaleCount_; }

        /** Query the number of filter directions available.
         */
    int directionCount() const
        { return directionCount_; }

        /** Change the number of directions / scales. This causes the
            recalculation of all filters.
         */
    void setDirectionScaleCounts(int directionCount, int scaleCount)
    {
        resize(directionCount * scaleCount);
        scaleCount_ = scaleCount;
        directionCount_ = directionCount;
        initFilters();
    }

        /** Return the center frequency of the filter(s) with
            scale==0. Filters with scale>0 will have a center frequency
            reduced in octaves:
            <tt>centerFrequency= maxCenterFrequency / 2.0^scale</tt>
        */
    double maxCenterFrequency()
        { return maxCenterFrequency_; }

        /** Change the center frequency of the filter(s) with
            scale==0. See \ref maxCenterFrequency().
         */
    void setMaxCenterFrequency(double maxCenterFrequency)
    {
        maxCenterFrequency_ = maxCenterFrequency;
        initFilters();
    }
};

//@}

} // namespace vigra

#endif // VIGRA_GABORFILTER_HXX
