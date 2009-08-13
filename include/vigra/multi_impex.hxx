/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Gunnar Kedenburg                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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


#ifndef VIGRA_MULTI_IMPEX_HXX
#define VIGRA_MULTI_IMPEX_HXX

#include <memory>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>

#include "config.hxx"
#include "basicimageview.hxx"
#include "impex.hxx"
#include "multi_array.hxx"
#include "multi_pointoperators.hxx"

#ifdef _MSC_VER
# include <direct.h>
#else
# include <unistd.h>
#endif

namespace vigra {

/** \addtogroup VolumeImpex Import/export of volume data.
*/

//@{

/** \brief Argument object for the function importVolume().

    See \ref importVolume() for usage example. This object can be used
    to define the properties of a volume data set to be read from disk.
    Sorry, no \ref detailedDocumentation() available yet.

    <b>\#include</b> \<<a href="imageinfo_8hxx-source.html">vigra/multi_impex.hxx</a>\><br>
    Namespace: vigra
**/
class VolumeImportInfo
{
  public:
    typedef ImageImportInfo::PixelType PixelType;

        /// type of volume size returned by shape()
    typedef MultiArrayShape<3>::type   ShapeType;

        /// provided for backwards-compatibility (deprecated)
    typedef ShapeType                  size_type;

        /// 3D resolution type returned by resolution()
    typedef TinyVector<float, 3>       Resolution;

    VIGRA_EXPORT VolumeImportInfo(const std::string &filename);
    VIGRA_EXPORT VolumeImportInfo(const std::string &baseName, const std::string &extension);

    VIGRA_EXPORT ShapeType shape() const;

        /** Get width of the volume.
         **/
    VIGRA_EXPORT MultiArrayIndex width() const;

        /** Get height of the volume.
         **/
    VIGRA_EXPORT MultiArrayIndex height() const;

        /** Get depth of the volume.
         **/
    VIGRA_EXPORT MultiArrayIndex depth() const;

        /**
         * resolution() contains the alignment and resolution of the
         * volume.  resolution()[0] is the x increment in a left-handed
         * world coordinate system of one unstrided step in the volume
         * memory.  The [1] and [2] elements contain the y resp. z
         * increments of the strided row resp. slice steps in the
         * volume.
         *
         * EXAMPLES: (1.f, 1.f, 4.f) means that the slices are four
         * times thicker than the x/y resolution.
         * (1.f, -1.f, 1.f) means that the volume coordinate system is
         * right-handed.
         */
    VIGRA_EXPORT Resolution resolution() const;

        /** Query the pixel type of the image.

            Possible values are:
            <DL>
            <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"INT16"<DD> 16-bit signed integer (short)
            <DT>"UINT16"<DD> 16-bit unsigned integer (unsigned short)
            <DT>"INT32"<DD> 32-bit signed integer (long)
            <DT>"UINT32"<DD> 32-bit unsigned integer (unsigned long)
            <DT>"FLOAT"<DD> 32-bit floating point (float)
            <DT>"DOUBLE"<DD> 64-bit floating point (double)
            </DL>
         **/
    VIGRA_EXPORT const char * getPixelType() const;

        /** Query the pixel type of the image.

            Same as getPixelType(), but the result is returned as a 
            ImageImportInfo::PixelType enum. This is useful to implement
            a switch() on the pixel type.

            Possible values are:
            <DL>
            <DT>UINT8<DD> 8-bit unsigned integer (unsigned char)
            <DT>INT16<DD> 16-bit signed integer (short)
            <DT>UINT16<DD> 16-bit unsigned integer (unsigned short)
            <DT>INT32<DD> 32-bit signed integer (long)
            <DT>UINT32<DD> 32-bit unsigned integer (unsigned long)
            <DT>FLOAT<DD> 32-bit floating point (float)
            <DT>DOUBLE<DD> 64-bit floating point (double)
            </DL>
         **/
    VIGRA_EXPORT PixelType pixelType() const;

    VIGRA_EXPORT MultiArrayIndex numBands() const;
    VIGRA_EXPORT bool isGrayscale() const;
    VIGRA_EXPORT bool isColor() const;

    // get base file name without path, image index, and extension
    VIGRA_EXPORT const std::string &name() const;

    VIGRA_EXPORT const std::string &description() const;

    template <class T, class Stride>
    void importImpl(MultiArrayView <3, T, Stride> &volume) const;

  protected:
    void getVolumeInfoFromFirstSlice(const std::string &filename);

    size_type shape_;
    Resolution resolution_;
    //PixelType pixelType_;
    int numBands_;

    std::string path_, name_, description_, pixelType_;

    std::string rawFilename_;
    std::string baseName_, extension_;
    std::vector<std::string> numbers_;
};

/********************************************************/
/*                                                      */
/*                   VolumeExportInfo                    */
/*                                                      */
/********************************************************/

/** \brief Argument object for the function exportVolume().

    See \ref exportVolume() for usage example. This object must be used
    to define the properties of a volume to be written to disk.

    <b>\#include</b> \<<a href="imageinfo_8hxx-source.html">vigra/imageinfo.hxx</a>\><br>
    Namespace: vigra
**/
class VolumeExportInfo
{
  public:
        /** Construct VolumeExportInfo object.

            The volume will be stored in a by-slice manner, where the number of slices 
			equals the depth of the volume. The file names will be enumerated like
			<tt>name_base+"000"+name_ext</tt>, <tt>name_base+"001"+name_ext</tt> etc.
			(the actual number of zeros depends on the depth). If the target image type
			does not support the source voxel type, all slices will be mapped 
			simultaneously to the appropriate target range.
            The file type will be guessed from the extension unless overridden
            by \ref setFileType(). Recognized extensions: '.bmp', '.gif',
            '.jpeg', '.jpg', '.p7', '.png', '.pbm', '.pgm', '.pnm', '.ppm', '.ras',
            '.tif', '.tiff', '.xv', '.hdr'.
            JPEG support requires libjpeg, PNG support requires libpng, and
            TIFF support requires libtiff.
         **/
    VIGRA_EXPORT VolumeExportInfo( const char * name_base, const char * name_ext );
    VIGRA_EXPORT ~VolumeExportInfo();

        /** Set volume file name base.

		**/
    VIGRA_EXPORT VolumeExportInfo & setFileNameBase(const char * name_base);
        /** Set volume file name extension.

            The file type will be guessed from the extension unless overridden
            by \ref setFileType(). Recognized extensions: '.bmp', '.gif',
            '.jpeg', '.jpg', '.p7', '.png', '.pbm', '.pgm', '.pnm', '.ppm', '.ras',
            '.tif', '.tiff', '.xv', '.hdr'.
            JPEG support requires libjpeg, PNG support requires libpng, and
            TIFF support requires libtiff.
		**/
    VIGRA_EXPORT VolumeExportInfo & setFileNameExt(const char * name_ext);
    VIGRA_EXPORT const char * getFileNameBase() const;
    VIGRA_EXPORT const char * getFileNameExt() const;

        /** Store volume as given file type.

            This will override any type guessed
            from the file name's extension. Recognized file types:

            <DL>
            <DT>"BMP"<DD> Microsoft Windows bitmap image file.
            <DT>"GIF"<DD> CompuServe graphics interchange format; 8-bit color.
            <DT>"JPEG"<DD> Joint Photographic Experts Group JFIF format;
            compressed 24-bit color (only available if libjpeg is installed).
            <DT>"PNG"<DD> Portable Network Graphic
            (only available if libpng is installed).
            <DT>"PBM"<DD> Portable bitmap format (black and white).
            <DT>"PGM"<DD> Portable graymap format (gray scale).
            <DT>"PNM"<DD> Portable anymap.
            <DT>"PPM"<DD> Portable pixmap format (color).
            <DT>"SUN"<DD> SUN Rasterfile.
            <DT>"TIFF"<DD> Tagged Image File Format.
            (only available if libtiff is installed.)
            <DT>"VIFF"<DD> Khoros Visualization image file.
            </DL>

            With the exception of TIFF, VIFF, PNG, and PNM all file types store
            1 byte (gray scale and mapped RGB) or 3 bytes (RGB) per
            pixel.

            PNG can store UInt8 and UInt16 values, and supports 1 and 3 channel
            images. One additional alpha channel is also supported.

            PNM can store 1 and 3 channel images with UInt8, UInt16 and UInt32
            values in each channel.

            TIFF and VIFF are aditionally able to store short and long
            integers (2 or 4 bytes) and real values (32 bit float and
            64 bit double) without conversion. So you will need to use
            TIFF or VIFF if you need to store images with high
            accuracy (the appropriate type to write is automatically
            derived from the image type to be exported). However, many
            other programs using TIFF (e.g. ImageMagick) have not
            implemented support for those pixel types.  So don't be
            surprised if the generated TIFF is not readable in some
            cases.  If this happens, export the image as 'unsigned
            char' or 'RGBValue\<unsigned char\>' by calling
            \ref ImageExportInfo::setPixelType().

            Support to reading and writing ICC color profiles is
            provided for TIFF, JPEG, and PNG images.
         **/
    VIGRA_EXPORT VolumeExportInfo & setFileType( const char * );
    VIGRA_EXPORT const char * getFileType() const;

        /** Set compression type.

            Recognized strings: "" (no compression), "LZW",
            "RunLength", "1" ... "100". A number is interpreted as the
            compression quality for JPEG compression. JPEG compression is
            supported by the JPEG and TIFF formats. "LZW" is only available
            if libtiff was installed with LZW enabled. By default, libtiff came
            with LZW disabled due to Unisys patent enforcement. In this case,
            VIGRA stores the image uncompressed.

                Valid Compression for TIFF files:
                  JPEG    jpeg compression, call setQuality as well!
                  RLE     runlength compression
                  LZW     lzw compression
                  DEFLATE deflate compression
         **/
    VIGRA_EXPORT VolumeExportInfo & setCompression( const char * );
    VIGRA_EXPORT const char * getCompression() const;

        /** Set the pixel type of the volume file(s). Possible values are:
            <DL>
            <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"INT16"<DD> 16-bit signed integer (short)
            <DT>"UINT16"<DD> 16-bit unsigned integer (unsigned short)
            <DT>"INT32"<DD> 32-bit signed integer (long)
            <DT>"UINT32"<DD> 32-bit unsigned integer (unsigned long)
            <DT>"FLOAT"<DD> 32-bit floating point (float)
            <DT>"DOUBLE"<DD> 64-bit floating point (double)
            </DL>
         **/
    VIGRA_EXPORT VolumeExportInfo & setPixelType( const char * );

        /** Get the pixel type of the images in the volume. Possible values are:
            <DL>
            <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"INT16"<DD> 16-bit signed integer (short)
            <DT>"INT32"<DD> 32-bit signed integer (long)
            <DT>"FLOAT"<DD> 32-bit floating point (float)
            <DT>"DOUBLE"<DD> 64-bit floating point (double)
            </DL>
         **/
    VIGRA_EXPORT const char * getPixelType() const;
    
    VIGRA_EXPORT VolumeExportInfo & setForcedRangeMapping(double fromMin, double fromMax,
                                                     double toMin, double toMax);    
    VIGRA_EXPORT bool hasForcedRangeMapping() const;
    VIGRA_EXPORT double getFromMin() const;
    VIGRA_EXPORT double getFromMax() const;
    VIGRA_EXPORT double getToMin() const;
    VIGRA_EXPORT double getToMax() const;
    
        /** Set the volume resolution in horizontal direction
         **/
    VIGRA_EXPORT VolumeExportInfo & setXResolution( float );
    VIGRA_EXPORT float getXResolution() const;

        /** Set the image resolution in vertical direction
         **/
    VIGRA_EXPORT VolumeExportInfo & setYResolution( float );
    VIGRA_EXPORT float getYResolution() const;

        /** Set the image resolution in depth direction
         **/
    VIGRA_EXPORT VolumeExportInfo & setZResolution( float );
    VIGRA_EXPORT float getZResolution() const;

		/** Set the position of the upper Left corner on a global
            canvas.

            Currently only supported by TIFF and PNG files.

            The offset is encoded in the XPosition and YPosition TIFF tags.

            @param pos     position of the upper left corner in pixels
                           (must be >= 0)
         **/
	// FIXME: mhanselm: we might want to support 3D positions
    VIGRA_EXPORT VolumeExportInfo & setPosition(const Diff2D & pos);

        /** Get the position of the upper left corner on
            a global canvas.
         **/
	// FIXME: mhanselm: we might want to support 3D positions
    VIGRA_EXPORT Diff2D getPosition() const;

        /**
          ICC profiles (handled as raw data so far).
          see getICCProfile()/setICCProfile()
         **/
    typedef ArrayVector<unsigned char> ICCProfile;

        /** Returns a reference to the ICC profile.
         */
    VIGRA_EXPORT const ICCProfile & getICCProfile() const;

        /** Sets the ICC profile.
            ICC profiles are currently supported by TIFF, PNG and JPEG images.
            (Otherwise, the profile data is silently ignored.)
         **/
    VIGRA_EXPORT VolumeExportInfo & setICCProfile(const ICCProfile & profile);

  private:
	float m_x_res, m_y_res, m_z_res;

	std::string m_filetype, m_filename_base, m_filename_ext, m_pixeltype, m_comp;
    Diff2D m_pos;
    ICCProfile m_icc_profile;
    double fromMin_, fromMax_, toMin_, toMax_;
};

namespace detail {

template <class DestIterator, class Shape, class T>
inline void
readVolumeImpl(DestIterator d, Shape const & shape, std::ifstream & s, ArrayVector<T> & buffer, MetaInt<0>)
{
    s.read((char*)buffer.begin(), shape[0]*sizeof(T));

    DestIterator dend = d + shape[0];
	int k = 0;
    for(; d < dend; ++d, k++)
    {
        *d = buffer[k];
    }
}

template <class DestIterator, class Shape, class T, int N>
void
readVolumeImpl(DestIterator d, Shape const & shape, std::ifstream & s, ArrayVector<T> & buffer, MetaInt<N>)
{
    DestIterator dend = d + shape[N];
    for(; d < dend; ++d)
    {
        readVolumeImpl(d.begin(), shape, s, buffer, MetaInt<N-1>());
    }
}

} // namespace detail

template <class T, class Stride>
void VolumeImportInfo::importImpl(MultiArrayView <3, T, Stride> &volume) const
{
	vigra_precondition(this->shape() == volume.shape(), "importVolume(): Volume must be shaped according to VolumeImportInfo.");

    if(rawFilename_.size())
    {
        std::string dirName, baseName;
        char oldCWD[2048];

#ifdef _MSC_VER
        if(_getcwd(oldCWD, 2048) == 0)
        {
            perror("getcwd");
            vigra_fail("VolumeImportInfo: Unable to query current directory (getcwd).");
        }
        if(_chdir(path_.c_str()))
        {
            perror("chdir");
            vigra_fail("VolumeImportInfo: Unable to change to new directory (chdir).");
        }
#else
        if(getcwd(oldCWD, 2048) == 0)
        {
            perror("getcwd");
            vigra_fail("VolumeImportInfo: Unable to query current directory (getcwd).");
        }
        if(chdir(path_.c_str()))
        {
            perror("chdir");
            vigra_fail("VolumeImportInfo: Unable to change to new directory (chdir).");
        }
#endif

        std::ifstream s(rawFilename_.c_str(), std::ios::binary);
	    vigra_precondition(s.good(), "RAW file could not be opened");

		ArrayVector<T> buffer(shape_[0]);
		detail::readVolumeImpl(volume.traverser_begin(), shape_, s, buffer, vigra::MetaInt<2>());

        //vigra_precondition(s.good(), "RAW file could not be opened");
        //s.read((char*)volume.data(), shape_[0]*shape_[1]*shape_[2]*sizeof(T));

#ifdef _MSC_VER
        if(_chdir(oldCWD))
            perror("chdir");
#else
        if(chdir(oldCWD))
            perror("chdir");
#endif

        vigra_postcondition(
            volume.shape() == shape(), "imported volume has wrong size");
    }
    else
    {
        for (unsigned int i = 0; i < numbers_.size(); ++i)
        {
            // build the filename
            std::string name = baseName_ + numbers_[i] + extension_;

            // import the image
            ImageImportInfo info (name.c_str ());

            // generate a basic image view to the current layer
            MultiArrayView <2, T, Stride> view (volume.bindOuter (i));
            vigra_precondition(view.shape() == info.shape(),
                "importVolume(): the images have inconsistent sizes.");

            importImage (info, destImage(view));
        }
    }
}


VIGRA_EXPORT void findImageSequence(const std::string &name_base,
                       const std::string &name_ext,
                       std::vector<std::string> & numbers);

/********************************************************/
/*                                                      */
/*                    importVolume                      */
/*                                                      */
/********************************************************/

/** \brief Function for importing a 3D volume.

    The data are expected to be stored in a by-slice manner,
    where the slices are enumerated from <tt>name_base+"[0-9]+"+name_ext</tt>.
    <tt>name_base</tt> may contain a path. All slice files with the same name base and
    extension are considered part of the same volume. Slice numbers must be non-negative,
    but can otherwise start anywhere and need not be successive. Slices will be read
    in ascending numerical (not lexicographic) order. All slices must have the
    same size. The <tt>volume</tt> will be reshaped to match the count and
    size of the slices found.

    <b>\#include</b>
    \<<a href="multi__impex_8hxx-source.html">vigra/multi_impex.hxx</a>\>

    Namespace: vigra
*/
template <class T, class Allocator>
void importVolume (MultiArray <3, T, Allocator> & volume,
                   const std::string &name_base,
                   const std::string &name_ext)
{
    VolumeImportInfo info(name_base, name_ext);
	volume.reshape(info.shape());

    info.importImpl(volume);
}


/** \brief Function for importing a 3D volume.

    The data can be given in two ways:

    <UL>
    <LI> If the volume is stored in a by-slice manner (e.g. one image per slice),
         the <tt>filename</tt> can refer to an arbitrary image from the set. <tt>importVolume()</tt>
         then assumes that the slices are enumerated like <tt>name_base+"[0-9]+"+name_ext</tt>,
         where <tt>name_base</tt>, the index, and <tt>name_ext</tt> are determined automatically.
         All slice files with the same name base and extension are considered part of the same
         volume. Slice numbers must be non-negative, but can otherwise start anywhere and need
         not be successive. Slices will be read in ascending numerical (not lexicographic) order.
         All slices must have the same size.
    <li> Otherwise, <tt>importVolume()</tt> will try to read <tt>filename</tt> as an
         info text file with the following key-value pairs:
         <UL>
         <LI> name = [short descriptive name of the volume] (optional)
         <LI> filename = [absolute or relative path to raw voxel data file] (required)
         <li> gradfile =  [absolute or relative path to gradient data file] (currently ignored)
         <li> description =  [arbitrary description of the data set] (optional)
         <li> width = [positive integer] (required)
         <li> height = [positive integer] (required)
         <li> depth = [positive integer] (required)
         <li> datatype = [UNSIGNED_CHAR | UNSIGNED_BYTE] (default: UNSIGNED_CHAR)
         </UL>
         The voxel type is currently assumed to be binary compatible to the <tt>value_type T</TT>
         of the <tt>MuliArray</tt>. Lines starting with "#" are ignored.
    </UL>

    In either case, the <tt>volume</tt> will be reshaped to match the count and
    size of the slices found.

    <b>\#include</b>
    \<<a href="multi__impex_8hxx-source.html">vigra/multi_impex.hxx</a>\>

    Namespace: vigra
*/
template <class T, class Allocator>
void importVolume(MultiArray <3, T, Allocator> &volume,
                  const std::string &filename)
{
    VolumeImportInfo info(filename);
	volume.reshape(info.shape());

    info.importImpl(volume);
}

/** \brief Function for importing a 3D volume.

    Read the volume data set <tt>info</tt> refers to. Explicit construction
    of the info object allows to allocate a <tt>volume</tt> object type whose
    <tt>value_type</tt> matches the voxel type of the stored data.
    The <tt>volume</tt> will be reshaped to match the count and
    size of the slices found.

    <b>\#include</b>
    \<<a href="multi__impex_8hxx-source.html">vigra/multi_impex.hxx</a>\>

    Namespace: vigra
*/
template <class T, class Stride>
void importVolume(VolumeImportInfo const & info, MultiArrayView <3, T, Stride> &volume)
{
    info.importImpl(volume);
}

namespace detail {

template <class T>
void setRangeMapping(std::string const & pixeltype,
                     FindMinMax<T> const & minmax, ImageExportInfo & info)
{
    if(pixeltype == "UINT8")
        info.setForcedRangeMapping(minmax.min, minmax.max,
                                   (double)NumericTraits<Int8>::min(),
                                   (double)NumericTraits<Int8>::max());
    else if(pixeltype == "INT16")
        info.setForcedRangeMapping(minmax.min, minmax.max,
                                   (double)NumericTraits<Int16>::min(),
                                   (double)NumericTraits<Int16>::max());
    else if(pixeltype == "UINT16")
        info.setForcedRangeMapping(minmax.min, minmax.max,
                                   (double)NumericTraits<UInt16>::min(),
                                   (double)NumericTraits<UInt16>::max());
    else if(pixeltype == "INT32")
        info.setForcedRangeMapping(minmax.min, minmax.max,
                                   (double)NumericTraits<Int32>::min(),
                                   (double)NumericTraits<Int32>::max());
    else if(pixeltype == "UINT32")
        info.setForcedRangeMapping(minmax.min, minmax.max,
                                   (double)NumericTraits<UInt32>::min(),
                                   (double)NumericTraits<UInt32>::max());
    else if(pixeltype == "FLOAT")
        info.setForcedRangeMapping(minmax.min, minmax.max, 0.0, 1.0);
    else if(pixeltype == "DOUBLE")
        info.setForcedRangeMapping(minmax.min, minmax.max, 0.0, 1.0);
}

template <class T, class Tag>
void setRangeMapping(MultiArrayView <3, T, Tag> const & volume,
                     ImageExportInfo & info, VigraTrueType /* isScalar */)
{
    std::string pixeltype = info.getPixelType();
    bool downcast = negotiatePixelType(getEncoderType(info.getFileName(), info.getFileType()),
                                       TypeAsString<T>::result(), pixeltype);

    if(downcast)
    {
        FindMinMax<T> minmax;
        inspectMultiArray(srcMultiArrayRange(volume), minmax);
        setRangeMapping(pixeltype, minmax, info);
    }
}

template <class T, class Tag>
void setRangeMapping(MultiArrayView <3, T, Tag> const & volume,
                     ImageExportInfo & info, VigraFalseType /* isScalar */)
{
    typedef typename T::value_type SrcComponent;
    std::string pixeltype = info.getPixelType();
    bool downcast = negotiatePixelType(getEncoderType(info.getFileName(), info.getFileType()),
                                       TypeAsString<SrcComponent>::result(), pixeltype);

    if(downcast)
    {
        unsigned int bands = volume(0,0,0).size();
        FindMinMax<SrcComponent> minmax;
        for(unsigned int i=0; i<bands; ++i)
        {
            VectorComponentValueAccessor<T> band(i);
            inspectMultiArray(srcMultiArrayRange(volume, band), minmax );
        }
        setRangeMapping(pixeltype, minmax, info);
    }
}

} // namespace detail

/********************************************************/
/*                                                      */
/*                    exportVolume                      */
/*                                                      */
/********************************************************/

/** \brief Function for exporting a 3D volume.

    The volume is exported in a by-slice manner, where the number of slices equals
    the depth of the volume. The file names will be enumerated like
    <tt>name_base+"000"+name_ext</tt>, <tt>name_base+"001"+name_ext</tt> etc.
    (the actual number of zeros depends on the depth). If the target image type
    does not support the source voxel type, all slices will be mapped simultaneously
    to the appropriate target range.

    <b>\#include</b>
    \<<a href="multi__impex_8hxx-source.html">vigra/multi_impex.hxx</a>\>

    Namespace: vigra
*/
template <class T, class Tag>
void exportVolume (MultiArrayView <3, T, Tag> const & volume,
                   const VolumeExportInfo & volinfo)
{
	std::string name = std::string(volinfo.getFileNameBase()) + std::string(volinfo.getFileNameExt());
    ImageExportInfo info(name.c_str());
	info.setCompression(volinfo.getCompression());
	info.setPixelType(volinfo.getPixelType());
    detail::setRangeMapping(volume, info, typename NumericTraits<T>::isScalar());

    const unsigned int depth = volume.shape (2);
    int numlen = static_cast <int> (std::ceil (std::log10 ((double)depth)));
    for (unsigned int i = 0; i < depth; ++i)
    {

        // build the filename
        std::stringstream stream;
        stream << std::setfill ('0') << std::setw (numlen) << i;
        std::string name_num;
        stream >> name_num;
        std::string name = std::string(volinfo.getFileNameBase()) + name_num + std::string(volinfo.getFileNameExt());

        MultiArrayView <2, T, Tag> view (volume.bindOuter (i));
        //MultiArrayView <3, T, Tag> view (volume.bindAt (2, i));

        // export the image
        info.setFileName(name.c_str ());
        exportImage(srcImageRange(view), info);
    }
}

//@}

} // namespace vigra

#endif // VIGRA_MULTI_IMPEX_HXX
