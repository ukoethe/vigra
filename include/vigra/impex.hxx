/************************************************************************/
/*                                                                      */
/*               Copyright 2001-2002 by Gunnar Kedenburg                */
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

/*!
  \file  impex.hxx
  \brief image import and export functions

  this file provides the declarations and implementations of importImage()
  and exportImage(). the matching implementation for the given datatype is
  selected by template metacode.
*/

#ifndef VIGRA_IMPEX_HXX
#define VIGRA_IMPEX_HXX

#if defined(_MSC_VER)
#pragma warning (disable: 4267)
#endif

#include "vigra/stdimage.hxx"
#include "vigra/tinyvector.hxx"
#include "vigra/imageinfo.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/codec.hxx"
#include "vigra/accessor.hxx"
#include "vigra/inspectimage.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/copyimage.hxx"

// TODO
// next refactoring: pluggable conversion algorithms

namespace vigra
{
/** \addtogroup VigraImpex
**/
//@{

    /*!
      \brief used for reading bands after the source data type has been figured out.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template< class ImageIterator, class Accessor, class SrcValueType >
            void read_bands( Decoder * dec, ImageIterator ys, Accessor a, SrcValueType )
        }
        \endcode

      \param dec decoder object through which the source data will be accessed
      \param ys  image iterator referencing the upper left pixel of the destination image
      \param a   image accessor for the destination image
    */
    template< class ImageIterator, class Accessor, class SrcValueType >
    void read_bands( Decoder * dec, ImageIterator ys, Accessor a, SrcValueType )
    {
        typedef unsigned int size_type;
        typedef typename ImageIterator::row_iterator DstRowIterator;
        typedef typename Accessor::value_type  AccessorValueType;
        typedef typename AccessorValueType::value_type DstValueType;

        const size_type width = dec->getWidth();
        const size_type height = dec->getHeight();
        const size_type num_bands = dec->getNumBands();

        SrcValueType const * scanline;
        DstRowIterator xs;

        // iterate
        for( size_type y = 0; y < height; ++y, ++ys.y ) {
            dec->nextScanline();
            for( size_type b = 0; b < num_bands; ++b ) {
                xs = ys.rowIterator();
                scanline = static_cast< SrcValueType const * >
                    (dec->currentScanlineOfBand(b));
                for( size_type x = 0; x < width; ++x, ++xs ) {
                    a.setComponent( *scanline, xs, b );
                    scanline += dec->getOffset();
                }
            }
        }
    } // read_bands()

    /*!
      \brief used for reading bands after the source data type has been figured out.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template< class ImageIterator, class Accessor, class SrcValueType >
            void read_band( Decoder * dec, ImageIterator ys, Accessor a, SrcValueType )
        }
        \endcode

      \param dec decoder object through which the source data will be accessed
      \param ys  image iterator referencing the upper left pixel of the destination image
      \param a   image accessor for the destination image
    */
    template< class ImageIterator, class Accessor, class SrcValueType >
    void read_band( Decoder * dec, ImageIterator ys, Accessor a, SrcValueType )
    {
        typedef unsigned int size_type;
        typedef typename ImageIterator::row_iterator DstRowIterator;
        typedef typename Accessor::value_type DstValueType;
        const size_type width = dec->getWidth();
        const size_type height = dec->getHeight();

        SrcValueType const * scanline;
        DstRowIterator xs;

        for( size_type y = 0; y < height; ++y, ++ys.y ) {
            dec->nextScanline();
            xs = ys.rowIterator();
            scanline = static_cast< SrcValueType const * >(dec->currentScanlineOfBand(0));
            for( size_type x = 0; x < width; ++x, ++xs )
                a.set( scanline[x], xs );
        }
    } // read_band()

    /*!
      \brief used for reading images of vector type, such as integer of float rgb.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template< class ImageIterator, class Accessor >
            void importVectorImage( const ImageImportInfo & info, ImageIterator iter, Accessor a )
        }
        \endcode

      \param ImageIterator the image iterator type for the destination image
      \param Accessor      the image accessor type for the destination image
      \param info          user supplied image import information
      \param iter          image iterator referencing the upper left pixel of the destination image
      \param a             image accessor for the destination image
    */
    template< class ImageIterator, class Accessor >
    void importVectorImage( const ImageImportInfo & info, ImageIterator iter, Accessor a )
    {
        std::auto_ptr<Decoder> dec = decoder(info);
        std::string pixeltype = dec->getPixelType();

        if ( pixeltype == "UINT8" )
            read_bands( dec.get(), iter, a, (unsigned char)0 );
        else if ( pixeltype == "INT16" )
            read_bands( dec.get(), iter, a, short() );
        else if ( pixeltype == "INT32" )
            read_bands( dec.get(), iter, a, int() );
        else if ( pixeltype == "FLOAT" )
            read_bands( dec.get(), iter, a, float() );
        else if ( pixeltype == "DOUBLE" )
            read_bands( dec.get(), iter, a, double() );
        else
            vigra_precondition( false, "invalid pixeltype" );

        // close the decoder
        dec->close();
    }

    /*!
      \brief used for reading images of  scalar type, such as integer and float grayscale.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template < class ImageIterator, class Accessor >
            void importScalarImage( const ImageImportInfo & info, ImageIterator iter, Accessor a )
        }
        \endcode

      \param ImageIterator the image iterator type for the destination image
      \param Accessor      the image accessor type for the destination image
      \param info          user supplied image import information
      \param iter          image iterator referencing the upper left pixel of the destination image
      \param a             image accessor for the destination image
    */
    template < class ImageIterator, class Accessor >
    void importScalarImage( const ImageImportInfo & info, ImageIterator iter, Accessor a )
    {
        std::auto_ptr<Decoder> dec = decoder(info);
        std::string pixeltype = dec->getPixelType();

        if ( pixeltype == "UINT8" )
            read_band( dec.get(), iter, a, (unsigned char)0 );
        else if ( pixeltype == "INT16" )
            read_band( dec.get(), iter, a, short() );
        else if ( pixeltype == "INT32" )
            read_band( dec.get(), iter, a, int() );
        else if ( pixeltype == "FLOAT" )
            read_band( dec.get(), iter, a, float() );
        else if ( pixeltype == "DOUBLE" )
            read_band( dec.get(), iter, a, double() );
        else
            vigra_precondition( false, "invalid pixeltype" );

        // close the decoder
        dec->close();
    }

    template < class ImageIterator, class Accessor >
    void importImage( const ImageImportInfo & info, ImageIterator iter, Accessor a, VigraFalseType )
    {
        importVectorImage( info, iter, a );
    }

    template < class ImageIterator, class Accessor >
    void importImage( const ImageImportInfo & info, ImageIterator iter, Accessor a, VigraTrueType )
    {
        importScalarImage( info, iter, a );
    }

/********************************************************/
/*                                                      */
/*                     importImage                      */
/*                                                      */
/********************************************************/

    /** \brief Read an image, given an \ref vigra::ImageImportInfo object.

        <b> Declarations:</b>

        pass arguments explicitly:
        \code
        namespace vigra {
            template <class ImageIterator, class Accessor>
            void
            importImage(ImageImportInfo const & image, ImageIterator iter, Accessor a)
        }
        \endcode

        use argument objects in conjuction with \ref ArgumentObjectFactories:
        \code
        namespace vigra {
            template <class ImageIterator, class Accessor>
            inline void
            importImage(ImageImportInfo const & image, pair<ImageIterator, Accessor> dest)
        }
        \endcode

        <b> Usage:</b>

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        \code

        vigra::ImageImportInfo info("myimage.gif");

        if(info.isGrayscale())
        {
            // create byte image of appropriate size
            vigra::BImage in(info.width(), info.height());

            vigra::importImage(info, destImage(in)); // read the image
            ...
        }
        else
        {
            // create byte RGB image of appropriate size
            vigra::BRGBImage in(info.width(), info.height());

            vigra::importImage(info, destImage(in)); // read the image
            ...
        }

        \endcode

        <b> Preconditions:</b>

        <UL>

        <LI> the image file must be readable
        <LI> the file type must be one of

                <DL>
                <DT>"BMP"<DD> Microsoft Windows bitmap image file.
                <DT>"GIF"<DD> CompuServe graphics interchange format; 8-bit color.
                <DT>"JPEG"<DD> Joint Photographic Experts Group JFIF format; compressed 24-bit color. (only available if libjpeg is installed)
                <DT>"PNG"<DD> Portable Network Graphic. (only available if libpng is installed)
                <DT>"PBM"<DD> Portable bitmap format (black and white).
                <DT>"PGM"<DD> Portable graymap format (gray scale).
                <DT>"PNM"<DD> Portable anymap.
                <DT>"PPM"<DD> Portable pixmap format (color).
                <DT>"SUN"<DD> SUN Rasterfile.
                <DT>"TIFF"<DD> Tagged Image File Format. (only available if libtiff is installed.)
                <DT>"VIFF"<DD> Khoros Visualization image file.
                </DL>
        </UL>
    **/
    template < class ImageIterator, class Accessor >
    void importImage( const ImageImportInfo & info, ImageIterator iter, Accessor a )
    {
        typedef typename NumericTraits<typename Accessor::value_type>::isScalar is_scalar;
        importImage( info, iter, a, is_scalar() );
    }

    template < class ImageIterator, class Accessor >
    void importImage( const ImageImportInfo & info, pair< ImageIterator, Accessor > dest )
    {
        importImage( info, dest.first, dest.second );
    }

    /*!
      \brief used for writing bands after the source data type has been figured out.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template< class ImageIterator, class Accessor, class DstValueType >
            void write_bands( Encoder * enc, ImageIterator ul, ImageIterator lr, Accessor a, DstValueType )
        }
        \endcode

      \param enc encoder object through which the destination data will be accessed
      \param ul  image iterator referencing the upper left pixel of the source image
      \param lr  image iterator referencing the lower right pixel of the source image
      \param a   image accessor for the source image
    */
    template< class ImageIterator, class Accessor, class DstValueType >
    void write_bands( Encoder * enc, ImageIterator ul, ImageIterator lr, Accessor a, DstValueType)
    {
        typedef unsigned int size_type;
        typedef typename ImageIterator::row_iterator SrcRowIterator;
        typedef typename Accessor::value_type  AccessorValueType;
        typedef typename AccessorValueType::value_type SrcValueType;

        // complete decoder settings
        const size_type width = lr.x - ul.x;
        const size_type height = lr.y - ul.y;
        enc->setWidth(width);
        enc->setHeight(height);
        const size_type num_bands = a(ul).size();
        enc->setNumBands(num_bands);
        enc->finalizeSettings();

        SrcRowIterator xs;
        DstValueType * scanline;

        // iterate
        ImageIterator ys(ul);
        for( size_type y = 0; y < height; ++y, ++ys.y ) {
            for( size_type b = 0; b < num_bands; ++b ) {
                xs = ys.rowIterator();
                scanline = static_cast< DstValueType * >
                    (enc->currentScanlineOfBand(b));
                for( size_type x = 0; x < width; ++x, ++xs ) {
                    *scanline = a.getComponent( xs, b );
                    scanline += enc->getOffset();
                }
            }
            enc->nextScanline();
        }
    } // write_bands()

    /*!
      \brief used for writing bands after the source data type has been figured out.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template< class ImageIterator, class Accessor, class DstValueType >
            void write_band( Encoder * enc, ImageIterator ul, ImageIterator lr, Accessor a, DstValueType )
        }
        \endcode

      \param enc encoder object through which the destination data will be accessed
      \param ul  image iterator referencing the upper left pixel of the source image
      \param lr  image iterator referencing the lower right pixel of the source image
      \param a   image accessor for the source image
    */
    template< class ImageIterator, class Accessor, class DstValueType >
    void write_band( Encoder * enc, ImageIterator ul, ImageIterator lr, Accessor a, DstValueType)
    {
        typedef unsigned int size_type;
        typedef typename ImageIterator::row_iterator SrcRowIterator;
        typedef typename Accessor::value_type SrcValueType;

        // complete decoder settings
        const size_type width = lr.x - ul.x;
        const size_type height = lr.y - ul.y;
        enc->setWidth(width);
        enc->setHeight(height);
        enc->setNumBands(1);
        enc->finalizeSettings();

        SrcRowIterator xs;
        DstValueType * scanline;

        // iterate
        ImageIterator ys(ul);
        size_type y;
        for(  y = 0; y < height; ++y, ++ys.y ) {
            xs = ys.rowIterator();
            scanline = static_cast< DstValueType * >(enc->currentScanlineOfBand(0));
            for( size_type x = 0; x < width; ++x, ++xs, ++scanline )
                *scanline = a(xs);
            enc->nextScanline();
        }
    } // write_band()

namespace detail {

    template <class VECTOR_ACCESSOR>
    struct VectorComponentSelector
    {
        VECTOR_ACCESSOR inner_;
        int index_;
        
        typedef typename VECTOR_ACCESSOR::component_type component_type;

        VectorComponentSelector(VECTOR_ACCESSOR const & a, int i)
        : inner_(a),
          index_(i)
        {}
        
        template <class ITERATOR>
        component_type const & operator()(ITERATOR const & i) const 
        { 
            return inner_.getComponent(i, index_); 
        }
    };
    
    template < class SrcIterator, class SrcAccessor,
               class DestIterator, class DestAccessor >
    void mapVectorImageToByteImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                    DestIterator dul, DestAccessor dget )
    {
        typedef typename SrcAccessor::value_type SrcValue;
        typedef typename SrcValue::value_type SrcComponent;
        typedef typename NumericTraits<SrcValue>::RealPromote PromoteValue;
        typedef typename PromoteValue::value_type PromoteComponent;
        typedef detail::VectorComponentSelector<SrcAccessor> Selector;
        
        vigra::FindMinMax<SrcComponent> minmax;
        for(int i=0; i<sget.size(sul); ++i)
        {
            vigra::inspectImage( sul, slr, Selector(sget, i), minmax );
        }
        const PromoteComponent scale = 255.0 / (minmax.max - minmax.min);
        PromoteValue offset;
        for(int i=0; i<sget.size(sul); ++i)
        {
            offset[i] = -minmax.min;
        }
        vigra::transformImage( sul, slr, sget, dul, dget,
                               linearIntensityTransform( scale, offset ) );
    }
    
    template < class SrcIterator, class SrcAccessor,
               class DestIterator, class DestAccessor >
    void mapScalarImageToByteImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                    DestIterator dul, DestAccessor dget )
    {
        typedef typename SrcAccessor::value_type SrcValue;
        typedef typename NumericTraits<SrcValue>::RealPromote PromoteValue;
        
        vigra::FindMinMax<SrcValue> minmax;
        vigra::inspectImage( sul, slr, sget, minmax );
        const PromoteValue scale = 255.0 / (minmax.max - minmax.min);
        const PromoteValue offset = -minmax.min;
        vigra::transformImage( sul, slr, sget, dul, dget,
                               linearIntensityTransform( scale, offset ) );
    }
    
    template < class SrcIterator, class SrcAccessor, class T >
    void exportImageWithConversion(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                   Encoder * enc, const char * type, T zero, VigraTrueType /* is_scalar */ )
    {
        if ( isPixelTypeSupported( enc->getFileType(), type ) ) {
            enc->setPixelType( type );
            write_band( enc, sul, slr, sget, zero );
        } else {
            // convert to unsigned char in the usual way
            enc->setPixelType( "UINT8" );
            vigra::BImage image(slr-sul);
            mapScalarImageToByteImage(sul, slr, sget, image.upperLeft(), image.accessor());
            write_band( enc, image.upperLeft(),
                        image.lowerRight(), image.accessor(), (unsigned char)0 );
        }
    }
        
    template < class SrcIterator, class SrcAccessor, class T >
    void exportImageWithConversion(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                   Encoder * enc, const char * type, T zero, VigraFalseType /* is_scalar */ )
    {
        if ( isPixelTypeSupported( enc->getFileType(), type ) ) {
            enc->setPixelType( type );
            write_bands( enc, sul, slr, sget, zero );
        } else {
            // convert to unsigned char in the usual way
            enc->setPixelType( "UINT8" );
            vigra::BRGBImage image(slr-sul);
            mapVectorImageToByteImage(sul, slr, sget, image.upperLeft(), image.accessor());
            write_bands( enc, image.upperLeft(),
                         image.lowerRight(), image.accessor(), (unsigned char)0 );
        }
    }
        
    template <class IsIntegral, int SIZE>
    struct ExportImageImpl;
    
    template <>
    struct ExportImageImpl<VigraTrueType, 1>
    {
        template < class SrcIterator, class SrcAccessor >
        static void exec(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         Encoder * enc, VigraTrueType /* is_scalar */)
        {
            enc->setPixelType( "UINT8" );
            write_band( enc, sul, slr, sget, (unsigned char)0 );
        }
        
        template < class SrcIterator, class SrcAccessor >
        static void exec(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         Encoder * enc, VigraFalseType /* is_scalar */)
        {
            enc->setPixelType( "UINT8" );
            write_bands( enc, sul, slr, sget, (unsigned char)0 );
        }
    };
    
    template <>
    struct ExportImageImpl<VigraTrueType, 2>
    {
        template < class SrcIterator, class SrcAccessor, class IsScalar >
        static void exec(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         Encoder * enc, IsScalar isScalar)
        {
            exportImageWithConversion( sul, slr, sget, enc, "INT16", short(), isScalar);
        }
    };
    
    template <>
    struct ExportImageImpl<VigraTrueType, 4>
    {
        template < class SrcIterator, class SrcAccessor, class IsScalar >
        static void exec(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         Encoder * enc, IsScalar isScalar)
        {
            exportImageWithConversion( sul, slr, sget, enc, "INT32", int(), isScalar );
        }
    };
    
    template <>
    struct ExportImageImpl<VigraFalseType, 4>
    {
        template < class SrcIterator, class SrcAccessor, class IsScalar >
        static void exec(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         Encoder * enc, IsScalar isScalar)
        {
            exportImageWithConversion( sul, slr, sget, enc, "FLOAT", float(), isScalar );
        }
    };
    
    template <>
    struct ExportImageImpl<VigraFalseType, 8>
    {
        template < class SrcIterator, class SrcAccessor, class IsScalar >
        static void exec(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         Encoder * enc, IsScalar isScalar)
        {
            exportImageWithConversion( sul, slr, sget, enc, "DOUBLE", double(), isScalar );
        }
    };
    
} // namespace detail
    

    /*!
      \brief used for writing images of floating point vector type, such as floating point rgb.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template < class SrcIterator, class SrcAccessor >
            void exportFloatingVectorImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                            const ImageExportInfo & info )
        }
        \endcode

      \param SrcIterator   the image iterator type for the source image
      \param SrcAccessor   the image accessor type for the source image
      \param sul           image iterator referencing the upper left pixel of the source image
      \param slr           image iterator referencing the lower right pixel of the source image
      \param sget          image accessor for the source image
      \param info          user supplied image export information
    */
    template < class SrcIterator, class SrcAccessor >
    void exportFloatingVectorImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                    const ImageExportInfo & info )
    {
        typedef typename SrcAccessor::value_type AccessorValueType;
        typedef typename AccessorValueType::value_type SrcValueType;

        std::auto_ptr<Encoder> enc = encoder(info);

        detail::ExportImageImpl<VigraFalseType, sizeof(SrcValueType)>::exec( sul, slr, sget, enc.get(), VigraFalseType() );

        // close the encoder
        enc->close();
    }

    /*!
      \brief used for writing images of integral vector type, such as integer rgb.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template < class SrcIterator, class SrcAccessor >
            void exportIntegralVectorImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                            const ImageExportInfo & info )
        }
        \endcode

      \param SrcIterator   the image iterator type for the source image
      \param SrcAccessor   the image accessor type for the source image
      \param sul           image iterator referencing the upper left pixel of the source image
      \param slr           image iterator referencing the lower right pixel of the source image
      \param sget          image accessor for the source image
      \param info          user supplied image export information
    */
    template < class SrcIterator, class SrcAccessor >
    void exportIntegralVectorImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                    const ImageExportInfo & info )
    {
        typedef typename SrcAccessor::value_type AccessorValueType;
        typedef typename AccessorValueType::value_type SrcValueType;

        std::auto_ptr<Encoder> enc = encoder(info);

        detail::ExportImageImpl<VigraTrueType, sizeof(SrcValueType)>::exec( sul, slr, sget, enc.get(), VigraFalseType() );

        // close the encoder
        enc->close();
    }

    /*!
      \brief used for writing images of floating point scalar type, such as floating point grayscale.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template < class SrcIterator, class SrcAccessor >
            void exportFloatingScalarImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                            const ImageExportInfo & info )
        }
        \endcode

      \param SrcIterator   the image iterator type for the source image
      \param SrcAccessor   the image accessor type for the source image
      \param sul           image iterator referencing the upper left pixel of the source image
      \param slr           image iterator referencing the lower right pixel of the source image
      \param sget          image accessor for the source image
      \param info          user supplied image export information
    */
    template < class SrcIterator, class SrcAccessor >
    void exportFloatingScalarImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                    const ImageExportInfo & info )
    {
        typedef typename SrcAccessor::value_type SrcValueType;

        std::auto_ptr<Encoder> enc = encoder(info);

        detail::ExportImageImpl<VigraFalseType, sizeof(SrcValueType)>::exec( sul, slr, sget, enc.get(), VigraTrueType() );

        // close the encoder
        enc->close();
    }

    /*!
      \brief used for writing images of integral scalar type, such as integer grayscale.

        <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
        Namespace: vigra

        <b> Declaration:</b>

        \code
        namespace vigra {
            template < class SrcIterator, class SrcAccessor >
            void exportIntegralScalarImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                            const ImageExportInfo & info )
        }
        \endcode

      \param SrcIterator   the image iterator type for the source image
      \param SrcAccessor   the image accessor type for the source image
      \param sul           image iterator referencing the upper left pixel of the source image
      \param slr           image iterator referencing the lower right pixel of the source image
      \param sget          image accessor for the source image
      \param info          user supplied image export information
    */
    template < class SrcIterator, class SrcAccessor >
    void exportIntegralScalarImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                                    const ImageExportInfo & info )
    {
        typedef typename SrcAccessor::value_type SrcValueType;

        std::auto_ptr<Encoder> enc = encoder(info);

        detail::ExportImageImpl<VigraTrueType, sizeof(SrcValueType)>::exec( sul, slr, sget, enc.get(), VigraTrueType() );

        enc->close();
    }

    template < class SrcIterator, class SrcAccessor >
    inline
    void exportImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                      const ImageExportInfo & info, VigraFalseType is_scalar)
    {
        typedef typename SrcAccessor::value_type AccessorValueType;
        typedef typename AccessorValueType::value_type SrcValueType;
        typedef typename NumericTraits<SrcValueType>::isIntegral is_integral;
        std::auto_ptr<Encoder> enc = encoder(info);
        detail::ExportImageImpl<is_integral, sizeof(SrcValueType)>::exec( sul, slr, sget, enc.get(), is_scalar);
        enc->close();
    }

    template < class SrcIterator, class SrcAccessor >
    inline
    void exportImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                      const ImageExportInfo & info, VigraTrueType is_scalar )
    {
        typedef typename SrcAccessor::value_type SrcValueType;
        typedef typename NumericTraits<SrcValueType>::isIntegral is_integral;
        std::auto_ptr<Encoder> enc = encoder(info);
        detail::ExportImageImpl<is_integral, sizeof(SrcValueType)>::exec( sul, slr, sget, enc.get(), is_scalar );
        enc->close();
    }

/********************************************************/
/*                                                      */
/*                     exportImage                      */
/*                                                      */
/********************************************************/

/** \brief Write an image, given an \ref vigra::ImageExportInfo object.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor>
        void exportImage(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         ImageExportInfo const & info)
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor>
        void exportImage(SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                         ImageExportInfo const & info)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="impex_8hxx-source.html">vigra/impex.hxx</a>"<br>
    Namespace: vigra

    \code


    vigra::BRGBImage out(w, h);
    ...

    // write as JPEG image, using compression quality 80
    vigra::exportImage(srcImageRange(out),
                      vigra::ImageExportInfo("myimage.jpg").setCompression("80"));

    \endcode

    <b> Preconditions:</b>

    <UL>

    <LI> the image file must be writable
    <LI> the file type must be one of

            <DL>
            <DT>"BMP"<DD> Microsoft Windows bitmap image file.
            <DT>"GIF"<DD> CompuServe graphics interchange format; 8-bit color.
            <DT>"JPEG"<DD> Joint Photographic Experts Group JFIF format; compressed 24-bit color. (only available if libjpeg is installed)
            <DT>"PNG"<DD> Portable Network Graphic. (only available if libpng is installed)
            <DT>"PBM"<DD> Portable bitmap format (black and white).
            <DT>"PGM"<DD> Portable graymap format (gray scale).
            <DT>"PNM"<DD> Portable anymap.
            <DT>"PPM"<DD> Portable pixmap format (color).
            <DT>"SUN"<DD> SUN Rasterfile.
            <DT>"TIFF"<DD> Tagged Image File Format. (only available if libtiff is installed.)
            <DT>"VIFF"<DD> Khoros Visualization image file.
            </DL>

    </UL>
**/
    template < class SrcIterator, class SrcAccessor >
    inline
    void exportImage( SrcIterator sul, SrcIterator slr, SrcAccessor sget,
                      const ImageExportInfo & info )
    {
        typedef typename NumericTraits<typename SrcAccessor::value_type>::isScalar is_scalar;
        exportImage( sul, slr, sget, info, is_scalar() );
    }

    template < class SrcIterator, class SrcAccessor >
    inline
    void exportImage( triple<SrcIterator, SrcIterator, SrcAccessor> src,
                      const ImageExportInfo & info )
    {
        exportImage( src.first, src.second, src.third, info );
    }

//@}

} // namespace vigra

#endif /* VIGRA_IMPEX_HXX */
