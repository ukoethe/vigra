/************************************************************************/
/*                                                                      */
/*               Copyright 2012 Christoph Spiel                         */
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

#ifndef VIGRA_IMPEXALPHA_HXX
#define VIGRA_IMPEXALPHA_HXX

#include <vector>

#include "imageinfo.hxx"
#include "impex.hxx"
#include "impexbase.hxx"
#include "multi_shape.hxx"

namespace vigra
{
/** \addtogroup VigraImpex
 * @{
*/
    namespace detail
    {
        template <class ValueType,
                  class ImageIterator, class ImageAccessor,
                  class AlphaIterator, class AlphaAccessor>
        void
        read_image_band_and_alpha(Decoder* decoder,
                                  ImageIterator image_iterator, ImageAccessor image_accessor,
                                  AlphaIterator alpha_iterator, AlphaAccessor alpha_accessor)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;
            typedef typename AlphaIterator::row_iterator AlphaRowIterator;

            vigra_precondition(decoder->getNumExtraBands() == 1,
                               "vigra::detail::read_image_band_and_alpha: expecting exactly one alpha band");
            vigra_precondition(decoder->getNumBands() - decoder->getNumExtraBands() == 1,
                               "vigra::detail::read_image_band_and_alpha: expecting exactly one image band");

            const unsigned width(decoder->getWidth());
            const unsigned height(decoder->getHeight());
            const unsigned offset(decoder->getOffset());

            for (unsigned y = 0U; y != height; ++y)
            {
                decoder->nextScanline();

                const ValueType* scanline0 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));
                const ValueType* scanline1 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(1));

                ImageRowIterator is(image_iterator.rowIterator());
                const ImageRowIterator is_end(is + width);
                AlphaRowIterator as(alpha_iterator.rowIterator());

                while (is != is_end)
                {
                    image_accessor.set(*scanline0, is);
                    scanline0 += offset;
                    ++is;

                    alpha_accessor.set(*scanline1, as);
                    scanline1 += offset;
                    ++as;
                }

                ++image_iterator.y;
                ++alpha_iterator.y;
            }
        }


        template <class ValueType,
                  class ImageIterator, class ImageAccessor,
                  class AlphaIterator, class AlphaAccessor>
        void
        read_image_bands_and_alpha(Decoder* decoder,
                                   ImageIterator image_iterator, ImageAccessor image_accessor,
                                   AlphaIterator alpha_iterator, AlphaAccessor alpha_accessor)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;
            typedef typename AlphaIterator::row_iterator AlphaRowIterator;

            vigra_precondition(decoder->getNumExtraBands() == 1,
                               "vigra::detail::read_image_bands_and_alpha: expecting exactly one alpha band");
            vigra_precondition(decoder->getNumBands() - decoder->getNumExtraBands() == image_accessor.size(image_iterator),
                               "vigra::detail::read_image_bands_and_alpha: number of channels and image accessor do not match");

            const unsigned width(decoder->getWidth());
            const unsigned height(decoder->getHeight());
            const unsigned offset(decoder->getOffset());
            const unsigned accessor_size(image_accessor.size(image_iterator));

            // OPTIMIZATION: Specialization for the most common case
            // of an RGBA-image, i.e. three color channels plus one
            // alpha channel.
            if (accessor_size == 3U)
            {
                const ValueType* scanline_0;
                const ValueType* scanline_1;
                const ValueType* scanline_2;
                const ValueType* scanline_3; // alpha

                for (unsigned y = 0U; y != height; ++y)
                {
                    decoder->nextScanline();

                    scanline_0 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));
                    scanline_1 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(1));
                    scanline_2 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(2));
                    scanline_3 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(3));

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);
                    AlphaRowIterator as(alpha_iterator.rowIterator());

                    while (is != is_end)
                    {
                        image_accessor.setComponent(*scanline_0, is, 0);
                        image_accessor.setComponent(*scanline_1, is, 1);
                        image_accessor.setComponent(*scanline_2, is, 2);
                        alpha_accessor.set(*scanline_3, as);
                        scanline_0 += offset;
                        scanline_1 += offset;
                        scanline_2 += offset;
                        scanline_3 += offset;

                        ++is;
                        ++as;
                    }

                    ++image_iterator.y;
                    ++alpha_iterator.y;
                }
            }
            else
            {
                std::vector<const ValueType*> scanlines(accessor_size + 1U);

                for (unsigned y = 0U; y != height; ++y)
                {
                    decoder->nextScanline();

                    for (unsigned i = 0U; i != accessor_size + 1U; ++i)
                    {
                        scanlines[i] = static_cast<const ValueType*>(decoder->currentScanlineOfBand(i));
                    }

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);
                    AlphaRowIterator as(alpha_iterator.rowIterator());

                    while (is != is_end)
                    {
                        for (unsigned i = 0U; i != accessor_size; ++i)
                        {
                            image_accessor.setComponent(*scanlines[i], is, static_cast<int>(i));
                            scanlines[i] += offset;
                        }
                        ++is;

                        alpha_accessor.set(*scanlines[accessor_size], as);
                        scanlines[accessor_size] += offset;
                        ++as;
                    }

                    ++image_iterator.y;
                    ++alpha_iterator.y;
                }
            }
        }


        template <class ImageIterator, class ImageAccessor,
                  class AlphaIterator, class AlphaAccessor>
        void
        importImageAlpha(const ImageImportInfo& import_info,
                         ImageIterator image_iterator, ImageAccessor image_accessor,
                         AlphaIterator alpha_iterator, AlphaAccessor alpha_accessor,
                         /* isScalar? */ VigraTrueType)
        {
            VIGRA_UNIQUE_PTR<Decoder> decoder(vigra::decoder(import_info));

            switch (pixel_t_of_string(decoder->getPixelType()))
            {
            case UNSIGNED_INT_8:
                read_image_band_and_alpha<UInt8>(decoder.get(),
                                                 image_iterator, image_accessor,
                                                 alpha_iterator, alpha_accessor);
                break;
            case UNSIGNED_INT_16:
                read_image_band_and_alpha<UInt16>(decoder.get(),
                                                  image_iterator, image_accessor,
                                                  alpha_iterator, alpha_accessor);
                break;
            case UNSIGNED_INT_32:
                read_image_band_and_alpha<UInt32>(decoder.get(),
                                                  image_iterator, image_accessor,
                                                  alpha_iterator, alpha_accessor);
                break;
            case SIGNED_INT_16:
                read_image_band_and_alpha<Int16>(decoder.get(),
                                                 image_iterator, image_accessor,
                                                 alpha_iterator, alpha_accessor);
                break;
            case SIGNED_INT_32:
                read_image_band_and_alpha<Int32>(decoder.get(),
                                                 image_iterator, image_accessor,
                                                 alpha_iterator, alpha_accessor);
                break;
            case IEEE_FLOAT_32:
                read_image_band_and_alpha<float>(decoder.get(),
                                                 image_iterator, image_accessor,
                                                 alpha_iterator, alpha_accessor);
                break;
            case IEEE_FLOAT_64:
                read_image_band_and_alpha<double>(decoder.get(),
                                                  image_iterator, image_accessor,
                                                  alpha_iterator, alpha_accessor);
                break;
            default:
                vigra_fail("vigra::detail::importImageAlpha<scalar>: not reached");
            }

            decoder->close();
        }


        template <class ImageIterator, class ImageAccessor,
                  class AlphaIterator, class AlphaAccessor>
        void
        importImageAlpha(const ImageImportInfo& import_info,
                         ImageIterator image_iterator, ImageAccessor image_accessor,
                         AlphaIterator alpha_iterator, AlphaAccessor alpha_accessor,
                         /* isScalar? */ VigraFalseType)
        {
            VIGRA_UNIQUE_PTR<Decoder> decoder(vigra::decoder(import_info));

            switch (pixel_t_of_string(decoder->getPixelType()))
            {
            case UNSIGNED_INT_8:
                read_image_bands_and_alpha<UInt8>(decoder.get(),
                                                  image_iterator, image_accessor,
                                                  alpha_iterator, alpha_accessor);
                break;
            case UNSIGNED_INT_16:
                read_image_bands_and_alpha<UInt16>(decoder.get(),
                                                   image_iterator, image_accessor,
                                                   alpha_iterator, alpha_accessor);
                break;
            case UNSIGNED_INT_32:
                read_image_bands_and_alpha<UInt32>(decoder.get(),
                                                   image_iterator, image_accessor,
                                                   alpha_iterator, alpha_accessor);
                break;
            case SIGNED_INT_16:
                read_image_bands_and_alpha<Int16>(decoder.get(),
                                                  image_iterator, image_accessor,
                                                  alpha_iterator, alpha_accessor);
                break;
            case SIGNED_INT_32:
                read_image_bands_and_alpha<Int32>(decoder.get(),
                                                  image_iterator, image_accessor,
                                                  alpha_iterator, alpha_accessor);
                break;
            case IEEE_FLOAT_32:
                read_image_bands_and_alpha<float>(decoder.get(),
                                                  image_iterator, image_accessor,
                                                  alpha_iterator, alpha_accessor);
                break;
            case IEEE_FLOAT_64:
                read_image_bands_and_alpha<double>(decoder.get(),
                                                   image_iterator, image_accessor,
                                                   alpha_iterator, alpha_accessor);
                break;
            default:
                vigra_fail("vigra::detail::importImageAlpha<non-scalar>: not reached");
            }

            decoder->close();
        }
    } // end namespace detail


    /**
    
    \brief Read the image specified by the given \ref
    vigra::ImageImportInfo object including its alpha channel. 
    
    See \ref importImage() for more information.
    
    <B>Declarations</B>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        importImageAlpha(ImageImportInfo const & import_info,
                         MultiArrayView<2, T1, S1> image,
                         MultiArrayView<2, T2, S2> alpha);
    }
    \endcode
    
    \deprecatedAPI{importImageAlpha}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
        namespace vigra {
            template <class ImageIterator, class ImageAccessor,
                      class AlphaIterator, class AlphaAccessor>
            void
            importImageAlpha(const ImageImportInfo& importInfo,
                             ImageIterator imageIterator, ImageAccessor imageAccessor,
                             AlphaIterator alphaIterator, AlphaAccessor alphaAccessor)
        }
    \endcode
    Use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
        namespace vigra {
            template <class ImageIterator, class ImageAccessor,
                      class AlphaIterator, class AlphaAccessor>
            void
            importImageAlpha(const ImageImportInfo& importInfo,
                             const pair<ImageIterator, ImageAccessor>& image,
                             const pair<AlphaIterator, AlphaAccessor>& alpha)
        }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <B>\#include \<vigra/impexalpha.hxx\></B><br/>
    Namespace: vigra

    \code
    typedef UInt8 value_t;
    ImageImportInfo info("zorro.tif");

    if (info.isGrayscale())
    {
        MultiArray<2, value_t> alpha(info.shape());
        MultiArray<2, value_t> image(info.shape());

        importImageAlpha(info, image, alpha);
        ...
    }
    else
    {
        MultiArray<2, value_t>            alpha(info.shape());
        MultiArray<2, RGBValue<value_t> > image(info.shape());

        importImageAlpha(info, image, alpha);
        ...
    }
    \endcode

    \deprecatedUsage{importImageAlpha}
    \code
    typedef UInt8 value_t;
    ImageImportInfo info("zorro.tif");

    if (info.isGrayscale())
    {
        BasicImage<value_t> alpha(info.size());
        BasicImage<value_t> image(info.size());

        importImageAlpha(info,
                         image.upperLeft(), image.accessor(),
                         alpha.upperLeft(), alpha.accessor());
        ...
    }
    else
    {
        BasicImage<value_t> alpha(info.size());
        BasicImage<vigra::RGBValue<value_t> > image(info.size());

        importImageAlpha(info,
                         image.upperLeft(), image.accessor(),
                         alpha.upperLeft(), alpha.accessor());
        ...
    }
    \endcode
    \deprecatedEnd
    
    <B>Preconditions</B>
    
    - The same preconditions hold as for importImage(), however the
      only image formats that support alpha channels are
      + TIFF and
      + PNG.
      In particular, JPEG does <B>not</B> support alpha channels.
    - The alpha channel always is scalar-valued, i.e. comprises a
      single band.
*/
    doxygen_overloaded_function(template <...> void importImageAlpha)


    template <class ImageIterator, class ImageAccessor,
              class AlphaIterator, class AlphaAccessor>
    inline void
    importImageAlpha(const ImageImportInfo& import_info,
                     ImageIterator image_iterator, ImageAccessor image_accessor,
                     AlphaIterator alpha_iterator, AlphaAccessor alpha_accessor)
    {
        typedef typename ImageAccessor::value_type ImageValueType;
        typedef typename vigra::NumericTraits<ImageValueType>::isScalar is_scalar;

        detail::importImageAlpha(import_info,
                                 image_iterator, image_accessor,
                                 alpha_iterator, alpha_accessor,
                                 is_scalar());
    }


    template <class ImageIterator, class ImageAccessor,
              class AlphaIterator, class AlphaAccessor>
    inline void
    importImageAlpha(ImageImportInfo const & import_info,
                     pair<ImageIterator, ImageAccessor> image,
                     pair<AlphaIterator, AlphaAccessor> alpha)
    {
        importImageAlpha(import_info,
                         image.first, image.second,
                         alpha.first, alpha.second);
    }

    template <class T1, class S1,
              class T2, class S2>
    inline void
    importImageAlpha(ImageImportInfo const & import_info,
                     MultiArrayView<2, T1, S1> image,
                     MultiArrayView<2, T2, S2> alpha)
    {
        vigra_precondition(import_info.shape() == image.shape() && import_info.shape() == alpha.shape(),
            "importImageAlpha(): shape mismatch between input and output.");
        importImageAlpha(import_info, destImage(image), destImage(alpha));
    }


    namespace detail
    {
        template<class ValueType,
                 class ImageIterator, class ImageAccessor, class ImageScaler,
                 class AlphaIterator, class AlphaAccessor, class AlphaScaler>
        void
        write_image_band_and_alpha(Encoder* encoder,
                                   ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                                   const ImageScaler& image_scaler,
                                   AlphaIterator alpha_upper_left, AlphaAccessor alpha_accessor,
                                   const AlphaScaler& alpha_scaler)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;
            typedef typename ImageAccessor::value_type ImageValueType;

            typedef typename AlphaIterator::row_iterator AlphaRowIterator;
            typedef typename AlphaAccessor::value_type AlphaValueType;

            typedef detail::RequiresExplicitCast<ValueType> explicit_cast;

            vigra_precondition(image_lower_right.x >= image_upper_left.x,
                               "vigra::detail::write_image_band_and_alpha: negative width");
            vigra_precondition(image_lower_right.y >= image_upper_left.y,
                               "vigra::detail::write_image_band_and_alpha: negative height");

            const unsigned width(static_cast<unsigned>(image_lower_right.x - image_upper_left.x));
            const unsigned height(static_cast<unsigned>(image_lower_right.y - image_upper_left.y));

            encoder->setWidth(width);
            encoder->setHeight(height);
            encoder->setNumBands(1 + 1);
            encoder->finalizeSettings();

            const unsigned offset(encoder->getOffset()); // correct offset only _after_ finalizeSettings()

            // IMPLEMENTATION NOTE: We avoid calling the default constructor
            // to allow classes ImageIterator and AlphaIterator that do not
            // define one.
            ImageIterator image_iterator(image_upper_left);
            AlphaIterator alpha_iterator(alpha_upper_left);

            for (unsigned y = 0U; y != height; ++y)
            {
                ValueType* scanline0 = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));
                ValueType* scanline1 = static_cast<ValueType*>(encoder->currentScanlineOfBand(1));

                ImageRowIterator is(image_iterator.rowIterator());
                const ImageRowIterator is_end(is + width);
                AlphaRowIterator as(alpha_iterator.rowIterator());

                while (is != is_end)
                {
                    *scanline0 = explicit_cast::cast(image_scaler(image_accessor(is)));
                    scanline0 += offset;
                    ++is;

                    *scanline1 = explicit_cast::cast(alpha_scaler(alpha_accessor(as)));
                    scanline1 += offset;
                    ++as;
                }

                encoder->nextScanline();

                ++image_iterator.y;
                ++alpha_iterator.y;
            }
        }


        template<class ValueType,
                 class ImageIterator, class ImageAccessor, class ImageScaler,
                 class AlphaIterator, class AlphaAccessor, class AlphaScaler>
        void
        write_image_bands_and_alpha(Encoder* encoder,
                                    ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                                    const ImageScaler& image_scaler,
                                    AlphaIterator alpha_upper_left, AlphaAccessor alpha_accessor,
                                    const AlphaScaler& alpha_scaler)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;
            typedef typename AlphaIterator::row_iterator AlphaRowIterator;
            typedef detail::RequiresExplicitCast<ValueType> explicit_cast;

            vigra_precondition(image_lower_right.x >= image_upper_left.x,
                               "vigra::detail::write_image_bands_and_alpha: negative width");
            vigra_precondition(image_lower_right.y >= image_upper_left.y,
                               "vigra::detail::write_image_bands_and_alpha: negative height");

            const unsigned width(static_cast<unsigned>(image_lower_right.x - image_upper_left.x));
            const unsigned height(static_cast<unsigned>(image_lower_right.y - image_upper_left.y));
            const unsigned accessor_size(image_accessor.size(image_upper_left));

            encoder->setWidth(width);
            encoder->setHeight(height);
            encoder->setNumBands(accessor_size + 1U);
            encoder->finalizeSettings();

            const unsigned offset(encoder->getOffset()); // correct offset only _after_ finalizeSettings()

            // IMPLEMENTATION NOTE: We avoid calling the default constructor
            // to allow classes ImageIterator and AlphaIterator that do not
            // define one.
            ImageIterator image_iterator(image_upper_left);
            AlphaIterator alpha_iterator(alpha_upper_left);

            // OPTIMIZATION: Specialization for the most common case
            // of an RGBA-image, i.e. three color channels plus one
            // alpha channel.
            if (accessor_size == 3U)
            {
                ValueType* scanline_0;
                ValueType* scanline_1;
                ValueType* scanline_2;
                ValueType* scanline_3; // alpha

                for (unsigned y = 0U; y != height; ++y)
                {
                    scanline_0 = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));
                    scanline_1 = static_cast<ValueType*>(encoder->currentScanlineOfBand(1));
                    scanline_2 = static_cast<ValueType*>(encoder->currentScanlineOfBand(2));
                    scanline_3 = static_cast<ValueType*>(encoder->currentScanlineOfBand(3));

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);
                    AlphaRowIterator as(alpha_iterator.rowIterator());

                    while (is != is_end)
                    {
                        *scanline_0 = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, 0)));
                        *scanline_1 = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, 1)));
                        *scanline_2 = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, 2)));
                        *scanline_3 = explicit_cast::cast(alpha_scaler(alpha_accessor(as)));
                        scanline_0 += offset;
                        scanline_1 += offset;
                        scanline_2 += offset;
                        scanline_3 += offset;

                        ++is;
                        ++as;
                    }

                    encoder->nextScanline();

                    ++image_iterator.y;
                    ++alpha_iterator.y;
                }
            }
            else
            {
                std::vector<ValueType*> scanlines(accessor_size + 1U);

                for (unsigned y = 0U; y != height; ++y)
                {
                    for (unsigned i = 0U; i != accessor_size + 1U; ++i)
                    {
                        scanlines[i] = static_cast<ValueType*>(encoder->currentScanlineOfBand(i));
                    }

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);
                    AlphaRowIterator as(alpha_iterator.rowIterator());

                    while (is != is_end)
                    {
                        for (unsigned i = 0U; i != accessor_size; ++i)
                        {
                            *scanlines[i] = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, static_cast<int>(i))));
                            scanlines[i] += offset;
                        }
                        ++is;

                        *scanlines[accessor_size] = explicit_cast::cast(alpha_scaler(alpha_accessor(as)));
                        scanlines[accessor_size] += offset;
                        ++as;
                    }

                    encoder->nextScanline();

                    ++image_iterator.y;
                    ++alpha_iterator.y;
                }
            }
        }


        template <class ImageIterator, class ImageAccessor,
                  class AlphaIterator, class AlphaAccessor>
        void
        exportImageAlpha(ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                         AlphaIterator alpha_upper_left, AlphaAccessor alpha_accessor,
                         const ImageExportInfo& export_info,
                         /* isScalar? */ VigraTrueType)
        {
            typedef typename ImageAccessor::value_type ImageValueType;

            VIGRA_UNIQUE_PTR<Encoder> encoder(vigra::encoder(export_info));

            std::string pixel_type(export_info.getPixelType());
            const bool downcast(negotiatePixelType(encoder->getFileType(), TypeAsString<ImageValueType>::result(), pixel_type));
            const pixel_t type(pixel_t_of_string(pixel_type));

            encoder->setPixelType(pixel_type);

            const range_t image_source_range(find_source_value_range(export_info,
                                                                     image_upper_left, image_lower_right, image_accessor));
            const range_t alpha_source_range(find_source_value_range(export_info,
                                                                     alpha_upper_left,
                                                                     alpha_upper_left + (image_lower_right - image_upper_left),
                                                                     alpha_accessor));
            const range_t destination_range(find_destination_value_range(export_info, type));

            if ((downcast || export_info.hasForcedRangeMapping()) &&
                (image_source_range.first != destination_range.first || image_source_range.second != destination_range.second ||
                 alpha_source_range.first != destination_range.first || alpha_source_range.second != destination_range.second))
            {
                const linear_transform image_rescaler(image_source_range, destination_range);
                const linear_transform alpha_rescaler(alpha_source_range, destination_range);

                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_band_and_alpha<UInt8>(encoder.get(),
                                                      image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                      alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case UNSIGNED_INT_16:
                    write_image_band_and_alpha<UInt16>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                       alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case UNSIGNED_INT_32:
                    write_image_band_and_alpha<UInt32>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                       alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case SIGNED_INT_16:
                    write_image_band_and_alpha<Int16>(encoder.get(),
                                                      image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                      alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case SIGNED_INT_32:
                    write_image_band_and_alpha<Int32>(encoder.get(),
                                                      image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                      alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case IEEE_FLOAT_32:
                    write_image_band_and_alpha<float>(encoder.get(),
                                                      image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                      alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case IEEE_FLOAT_64:
                    write_image_band_and_alpha<double>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                       alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                default:
                    vigra_fail("vigra::detail::exportImageAlpha<scalar>: not reached");
                }
            }
            else
            {
                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_band_and_alpha<UInt8>(encoder.get(),
                                                      image_upper_left, image_lower_right, image_accessor, identity(),
                                                      alpha_upper_left, alpha_accessor, identity());
                    break;
                case UNSIGNED_INT_16:
                    write_image_band_and_alpha<UInt16>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, identity(),
                                                       alpha_upper_left, alpha_accessor, identity());
                    break;
                case UNSIGNED_INT_32:
                    write_image_band_and_alpha<UInt32>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, identity(),
                                                       alpha_upper_left, alpha_accessor, identity());
                    break;
                case SIGNED_INT_16:
                    write_image_band_and_alpha<Int16>(encoder.get(),
                                                      image_upper_left, image_lower_right, image_accessor, identity(),
                                                      alpha_upper_left, alpha_accessor, identity());
                    break;
                case SIGNED_INT_32:
                    write_image_band_and_alpha<Int32>(encoder.get(),
                                                      image_upper_left, image_lower_right, image_accessor, identity(),
                                                      alpha_upper_left, alpha_accessor, identity());
                    break;
                case IEEE_FLOAT_32:
                    write_image_band_and_alpha<float>(encoder.get(),
                                                      image_upper_left, image_lower_right, image_accessor, identity(),
                                                      alpha_upper_left, alpha_accessor, identity());
                    break;
                case IEEE_FLOAT_64:
                    write_image_band_and_alpha<double>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, identity(),
                                                       alpha_upper_left, alpha_accessor, identity());
                    break;
                default:
                    vigra_fail("vigra::detail::exportImageAlpha<scalar>: not reached");
                }
            }

            encoder->close();
        }


        template <class ImageIterator, class ImageAccessor,
                  class AlphaIterator, class AlphaAccessor>
        void
        exportImageAlpha(ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                         AlphaIterator alpha_upper_left, AlphaAccessor alpha_accessor,
                         const ImageExportInfo& export_info,
                         /* isScalar? */ VigraFalseType)
        {
            typedef typename ImageAccessor::value_type ImageBaseType;
            typedef typename ImageBaseType::value_type ImageValueType;

            VIGRA_UNIQUE_PTR<Encoder> encoder(vigra::encoder(export_info));

            std::string pixel_type(export_info.getPixelType());
            const bool downcast(negotiatePixelType(encoder->getFileType(), TypeAsString<ImageValueType>::result(), pixel_type));
            const pixel_t type(pixel_t_of_string(pixel_type));

            encoder->setPixelType(pixel_type);

            vigra_precondition(isBandNumberSupported(encoder->getFileType(), image_accessor.size(image_upper_left)),
                               "exportImageAlpha(): file format does not support requested number of bands (color channels)");

            const range_t image_source_range(find_source_value_range(export_info,
                                                                     image_upper_left, image_lower_right, image_accessor));
            const range_t alpha_source_range(find_source_value_range(export_info,
                                                                     alpha_upper_left,
                                                                     alpha_upper_left + (image_lower_right - image_upper_left),
                                                                     alpha_accessor));
            const range_t destination_range(find_destination_value_range(export_info, pixel_t_of_string(pixel_type)));

            if ((downcast || export_info.hasForcedRangeMapping()) &&
                (image_source_range.first != destination_range.first || image_source_range.second != destination_range.second ||
                 alpha_source_range.first != destination_range.first || alpha_source_range.second != destination_range.second))
            {
                const linear_transform image_rescaler(image_source_range, destination_range);
                const linear_transform alpha_rescaler(alpha_source_range, destination_range);

                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_bands_and_alpha<UInt8>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                       alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case UNSIGNED_INT_16:
                    write_image_bands_and_alpha<UInt16>(encoder.get(),
                                                        image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                        alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case UNSIGNED_INT_32:
                    write_image_bands_and_alpha<UInt32>(encoder.get(),
                                                        image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                        alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case SIGNED_INT_16:
                    write_image_bands_and_alpha<Int16>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                       alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case SIGNED_INT_32:
                    write_image_bands_and_alpha<Int32>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                       alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case IEEE_FLOAT_32:
                    write_image_bands_and_alpha<float>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                       alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                case IEEE_FLOAT_64:
                    write_image_bands_and_alpha<double>(encoder.get(),
                                                        image_upper_left, image_lower_right, image_accessor, image_rescaler,
                                                        alpha_upper_left, alpha_accessor, alpha_rescaler);
                    break;
                default:
                    vigra_fail("vigra::detail::exportImageAlpha<non-scalar>: not reached");
                }
            }
            else
            {
                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_bands_and_alpha<UInt8>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, identity(),
                                                       alpha_upper_left, alpha_accessor, identity());
                    break;
                case UNSIGNED_INT_16:
                    write_image_bands_and_alpha<UInt16>(encoder.get(),
                                                        image_upper_left, image_lower_right, image_accessor, identity(),
                                                        alpha_upper_left, alpha_accessor, identity());
                    break;
                case UNSIGNED_INT_32:
                    write_image_bands_and_alpha<UInt32>(encoder.get(),
                                                        image_upper_left, image_lower_right, image_accessor, identity(),
                                                        alpha_upper_left, alpha_accessor, identity());
                    break;
                case SIGNED_INT_16:
                    write_image_bands_and_alpha<Int16>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, identity(),
                                                       alpha_upper_left, alpha_accessor, identity());
                    break;
                case SIGNED_INT_32:
                    write_image_bands_and_alpha<Int32>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, identity(),
                                                       alpha_upper_left, alpha_accessor, identity());
                    break;
                case IEEE_FLOAT_32:
                    write_image_bands_and_alpha<float>(encoder.get(),
                                                       image_upper_left, image_lower_right, image_accessor, identity(),
                                                       alpha_upper_left, alpha_accessor, identity());
                    break;
                case IEEE_FLOAT_64:
                    write_image_bands_and_alpha<double>(encoder.get(),
                                                        image_upper_left, image_lower_right, image_accessor, identity(),
                                                        alpha_upper_left, alpha_accessor, identity());
                    break;
                default:
                    vigra_fail("vigra::detail::exportImageAlpha<non-scalar>: not reached");
                }
            }

            encoder->close();
        }
    } // end namespace detail


    /**
    \brief Write the image and its alpha channel to a file.
    
    See \ref exportImage() for more information.
    
    <B>Declarations</B>
    
    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        exportImageAlpha(MultiArrayView<2, T1, S1> const & image,
                         MultiArrayView<2, T2, S2> const & alpha,
                         ImageExportInfo const & export_info);

        template <class T1, class S1,
                  class T2, class S2>
        void
        exportImageAlpha(MultiArrayView<2, T1, S1> const & image,
                         MultiArrayView<2, T2, S2> const & alpha,
                         char const * filename)

        template <class T1, class S1,
                  class T2, class S2>
        void
        exportImageAlpha(MultiArrayView<2, T1, S1> const & image,
                         MultiArrayView<2, T2, S2> const & alpha,
                         std::string const & filename)
    }
    \endcode
    
    \deprecatedAPI{exportImageAlpha}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
        namespace vigra {
            template <class ImageIterator, class ImageAccessor,
                      class AlphaIterator, class AlphaAccessor>
            void
            exportImageAlpha(ImageIterator imageUpperLeft, ImageIterator imageLowerRight, ImageAccessor imageAccessor,
                             AlphaIterator alphaUpperLeft, AlphaAccessor alphaAccessor,
                             const ImageExportInfo& exportInfo)
        }
    \endcode
    Use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
        namespace vigra {
        template <class ImageIterator, class ImageAccessor,
                  class AlphaIterator, class AlphaAccessor>
        void
        exportImageAlpha(const triple<ImageIterator, ImageIterator, ImageAccessor>& image,
                         const pair<AlphaIterator, AlphaAccessor>& alpha,
                         const ImageExportInfo& exportInfo)
        }
    \endcode
    \deprecatedEnd
    
    <b> Usage:</b>
    
    <B>\#include \<vigra/impexalpha.hxx\></B><br/>
    Namespace: vigra
    
    \code
    typedef UInt8 value_t;

    MultiArray<2, value_t>            alpha(width, height);
    MultiArray<2, RGBValue<value_t> > image(width, height);

   ... // do some image processing

    // specify the output filename 
    exportImageAlpha(image, alpha, "zorro.tif");
    
    // use a ImageExportInfo if you need more control over the export
    exportImageAlpha(image, alpha, ImageExportInfo("zorro.tif").setPixelType("FLOAT"));
   \endcode

    \deprecatedUsage{exportImageAlpha}
    \code
    typedef UInt8 value_t;
    ImageExportInfo info("zorro.tif");

    if (info.isGrayscale())
    {
        BasicImage<value_t> alpha;
        BasicImage<value_t> image;

        ...

        exportImageAlpha(image.upperLeft(), image.lowerRight(), image.accessor(),
                         alpha.upperLeft(), alpha.accessor(),
                         info);
    }
    else
    {
        BasicImage<value_t> alpha;
        BasicImage<vigra::RGBValue<value_t> > image;

        ...

        exportImageAlpha(image.upperLeft(), image.lowerRight(), image.accessor(),
                         alpha.upperLeft(), alpha.accessor(),
                         info);
    }
    \endcode
    \deprecatedEnd
    
    <B>Preconditions</B>
    
    - The same preconditions hold as for exportImage(), however the
      only image formats that support alpha channels are
      + TIFF and
      + PNG.
      In particular, JPEG does <B>not</B> support alpha channels.
    - The alpha channel always is scalar-valued, i.e. comprises a
      single band.
*/
    doxygen_overloaded_function(template <...> void exportImageAlpha)


    template <class ImageIterator, class ImageAccessor,
              class AlphaIterator, class AlphaAccessor>
    inline void
    exportImageAlpha(ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                     AlphaIterator alpha_upper_left, AlphaAccessor alpha_accessor,
                     const ImageExportInfo& export_info)
    {
        typedef typename ImageAccessor::value_type ImageValueType;
        typedef typename vigra::NumericTraits<ImageValueType>::isScalar is_scalar;

        try
        {
            detail::exportImageAlpha(image_upper_left, image_lower_right, image_accessor,
                                     alpha_upper_left, alpha_accessor,
                                     export_info,
                                     is_scalar());
        }
        catch (Encoder::TIFFCompressionException&)
        {
            ImageExportInfo info(export_info);

            info.setCompression("");
            detail::exportImageAlpha(image_upper_left, image_lower_right, image_accessor,
                                     alpha_upper_left, alpha_accessor,
                                     info,
                                     is_scalar());
        }
    }


    template <class ImageIterator, class ImageAccessor,
              class AlphaIterator, class AlphaAccessor>
    inline void
    exportImageAlpha(triple<ImageIterator, ImageIterator, ImageAccessor> image,
                     pair<AlphaIterator, AlphaAccessor> alpha,
                     ImageExportInfo const & export_info)
    {
        exportImageAlpha(image.first, image.second, image.third,
                         alpha.first, alpha.second,
                         export_info);
    }

    template <class T1, class S1,
              class T2, class S2>
    inline void
    exportImageAlpha(MultiArrayView<2, T1, S1> const & image,
                     MultiArrayView<2, T2, S2> const & alpha,
                     ImageExportInfo const & export_info)
    {
        exportImageAlpha(srcImageRange(image),
                         srcImage(alpha),
                         export_info);
    }

    template <class T1, class S1,
              class T2, class S2>
    inline void
    exportImageAlpha(MultiArrayView<2, T1, S1> const & image,
                     MultiArrayView<2, T2, S2> const & alpha,
                     char const * name)
    {
        ImageExportInfo export_info(name);
        exportImageAlpha(srcImageRange(image),
                         srcImage(alpha),
                         export_info);
    }

    template <class T1, class S1,
              class T2, class S2>
    inline void
    exportImageAlpha(MultiArrayView<2, T1, S1> const & image,
                     MultiArrayView<2, T2, S2> const & alpha,
                     std::string const & name)
    {
        ImageExportInfo export_info(name.c_str());
        exportImageAlpha(srcImageRange(image),
                         srcImage(alpha),
                         export_info);
    }

/** @} */
    
} // end namespace vigra

#endif // VIGRA_IMPEXALPHA_HXX
