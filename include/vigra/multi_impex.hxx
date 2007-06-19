/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Gunnar Kedenburg                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

class VolumeImportInfo
{
  public:
    typedef ImageImportInfo::PixelType PixelType;

    typedef TinyVector<int, 3>         size_type;

    typedef TinyVector<float, 3>       Resolution;

    VIGRA_EXPORT inline VolumeImportInfo(const std::string &filename);
    VIGRA_EXPORT inline VolumeImportInfo(const std::string &baseName, const std::string &extension);

    VIGRA_EXPORT size_type size() const { return size_; }

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
    VIGRA_EXPORT Resolution resolution() const { return resolution_; }

    VIGRA_EXPORT PixelType pixelType() const { return pixelType_; }

    VIGRA_EXPORT int numBands() const { return numBands_; }
    VIGRA_EXPORT bool isGrayscale() const { return numBands_ == 1; }
    VIGRA_EXPORT bool isColor() const { return numBands_ > 1; }

    // get base file name without path, image index, and extension
    VIGRA_EXPORT const std::string &name() const { return name_; }
    
    VIGRA_EXPORT const std::string &description() const { return description_; }

    template <class T, class Allocator>
    void importImpl(MultiArray <3, T, Allocator> &volume) const;

  protected:
    void getVolumeInfoFromFirstSlice(const std::string &filename);
    
    size_type size_;
    Resolution resolution_;
    PixelType pixelType_;
    int numBands_;

    std::string path_, name_, description_;

    std::string rawFilename_;
    std::string baseName_, extension_;
    std::vector<std::string> numbers_;
};

template <class T, class Allocator>
void VolumeImportInfo::importImpl(MultiArray <3, T, Allocator> &volume) const
{
    volume.reshape(this->size());
    
    if(rawFilename_.size())
    {
        std::string dirName, baseName;
        char oldCWD[2048];

#ifdef _MSC_VER
        _getcwd(oldCWD, 2048);
        if(_chdir(path_.c_str()))
            perror("chdir");
#else
        getcwd(oldCWD, 2048);
        if(chdir(path_.c_str()))
            perror("chdir");
#endif

        std::ifstream s(rawFilename_.c_str(), std::ios::binary);
        vigra_precondition(s.good(), "RAW file could not be opened");
        s.read((char*)volume.begin(), size_[0]*size_[1]*size_[2]*sizeof(T));

#ifdef _MSC_VER
        _chdir(oldCWD);
#else
        chdir(oldCWD);
#endif

        vigra_postcondition(
            volume.size() == size(), "imported volume has wrong size");
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
            MultiArrayView <2, T> array_view (volume.bindOuter (i));
            BasicImageView <T> view = makeBasicImageView (array_view);
            vigra_precondition(view.size() == info.size(),
                "importVolume(): image size mismatch.");

            importImage (info, destImage(view));
        }
    }
}


VIGRA_EXPORT void findImageSequence(const std::string &name_base,
                       const std::string &name_ext,
                       std::vector<std::string> & numbers);

/** \addtogroup VolumeImpex Import/export of volume data.
*/

//@{

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
    "<a href="multi__impex_8hxx-source.html">vigra/multi_impex.hxx</a>"

    Namespace: vigra
*/
template <class T, class Allocator>
void importVolume (MultiArray <3, T, Allocator> & volume,
                   const std::string &name_base,
                   const std::string &name_ext)
{
    VolumeImportInfo info(name_base, name_ext);

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
    "<a href="multi__impex_8hxx-source.html">vigra/multi_impex.hxx</a>"

    Namespace: vigra
*/
template <class T, class Allocator>
void importVolume(MultiArray <3, T, Allocator> &volume,
                  const std::string &filename)
{
    VolumeImportInfo info(filename);

    info.importImpl(volume);
}

/** \brief Function for importing a 3D volume.

    Read the volume data set <tt>info</tt> refers to. Explicit construction
    of the info object allows to allocate a <tt>volume</tt> object type whose
    <tt>value_type</tt> matches the voxel type of the stored data.
    The <tt>volume</tt> will be reshaped to match the count and
    size of the slices found.

    <b>\#include</b>
    "<a href="multi__impex_8hxx-source.html">vigra/multi_impex.hxx</a>"

    Namespace: vigra
*/
template <class T, class Allocator>
void importVolume(VolumeImportInfo const & info, MultiArray <3, T, Allocator> &volume)
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
    std::auto_ptr<Encoder> enc = encoder(info);
    bool downcast = negotiatePixelType(enc->getFileType(),
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
    typedef typename T::calue_type SrcComponent;
    std::string pixeltype = info.getPixelType();
    std::auto_ptr<Encoder> enc = encoder(info);
    bool downcast = negotiatePixelType(enc->getFileType(),
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
    "<a href="multi__impex_8hxx-source.html">vigra/multi_impex.hxx</a>"

    Namespace: vigra
*/
template <class T, class Tag>
void exportVolume (MultiArrayView <3, T, Tag> const & volume,
                   const std::string &name_base,
                   const std::string &name_ext)
{
    std::string name = name_base + name_ext;
    ImageExportInfo info(name.c_str ());
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
        std::string name = name_base + name_num + name_ext;

        // generate a basic image view to the current layer
        MultiArrayView <2, T, Tag> array_view (volume.bindOuter (i));
        BasicImageView <T> view = makeBasicImageView (array_view);

        // export the image
        info.setFileName(name.c_str ());
        exportImage(srcImageRange(view), info);
    }
}

//@}

} // namespace vigra

#endif // VIGRA_MULTI_IMPEX_HXX
