/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyimpex_PyArray_API
//#define NO_IMPORT_ARRAY

#include <Python.h>
#include <iostream>
#include <cstring>
#include <cstdio>
#include "vigra/numpy_array.hxx"
#include "vigra/impex.hxx"
#include "vigra/multi_impex.hxx"
#include <vigra/numpy_array_converters.hxx>

#ifdef HasHDF5
# include "vigra/hdf5impex.hxx"
#endif

namespace python = boost::python;

namespace vigra {

namespace detail {

template <class T>
NumpyAnyArray readImageImpl(ImageImportInfo const & info, std::string order = "")
{
    typedef UnstridedArrayTag Stride;
    
    if(order == "")
        order = detail::defaultOrder();

    switch(info.numBands())
    {
      case 1:
      {
        NumpyArray<2, Singleband<T>, Stride> res(MultiArrayShape<2>::type(info.width(), info.height()), order);
        importImage(info, destImage(res));
        return res;
      }
      case 2:
      {
        NumpyArray<2, TinyVector<T, 2>, Stride> res(MultiArrayShape<2>::type(info.width(), info.height()), order);
        importImage(info, destImage(res));
        return res;
      }
      case 3:
      {
        NumpyArray<2, RGBValue<T>, Stride> res(MultiArrayShape<2>::type(info.width(), info.height()), order);
        importImage(info, destImage(res));
        return res;
      }
      case 4:
      {
        NumpyArray<2, TinyVector<T, 4>, Stride> res(MultiArrayShape<2>::type(info.width(), info.height()), order);
        importImage(info, destImage(res));
        return res;
      }
      default:
      {
        NumpyArray<3, Multiband<T> > res(MultiArrayShape<3>::type(info.width(), info.height(), info.numBands()), order);
        importImage(info, destImage(res));
        return res;
      }
    }
}

std::string numpyTypeIdToImpexString(NPY_TYPES typeID)
{
    switch(typeID)
    {
        case NPY_BOOL: return std::string("UINT8");
        case NPY_INT8: return std::string("INT8");
        case NPY_UINT8: return std::string("UINT8");
        case NPY_INT16: return std::string("INT16");
        case NPY_UINT16: return std::string("UINT16");
        case NPY_INT32: return std::string("INT32");
        case NPY_UINT32: return std::string("UINT32");
        case NPY_INT64: return std::string("DOUBLE");
        case NPY_UINT64: return std::string("DOUBLE");
        case NPY_FLOAT32: return std::string("FLOAT");
        case NPY_FLOAT64: return std::string("DOUBLE");
        default: return std::string("UNKNOWN");
    }
}

} // namespace detail

NumpyAnyArray 
readImage(const char * filename, python::object import_type, std::string order = "")
{
    ImageImportInfo info(filename);
    std::string importType(info.getPixelType());
    
    if(python::extract<std::string>(import_type).check())
    {
        std::string type = python::extract<std::string>(import_type)();        
        if(type != "" && type != "NATIVE")
            importType = type;
    }
    else if(python::extract<NPY_TYPES>(import_type).check())
    {
        importType = detail::numpyTypeIdToImpexString(python::extract<NPY_TYPES>(import_type)());
    }
    else if(import_type)
        vigra_precondition(false, "readImage(filename, import_type, order): import_type must be a string or a numpy dtype.");

    // FIXME: support all types, at least via a type cast at the end?
    if(importType == "FLOAT")
        return detail::readImageImpl<float>(info, order);
    if(importType == "UINT8")
        return detail::readImageImpl<UInt8>(info, order);
    if(importType == "INT16")
        return detail::readImageImpl<Int16>(info, order);
    if(importType == "UINT16")
        return detail::readImageImpl<UInt16>(info, order);
    if(importType == "INT32")
        return detail::readImageImpl<Int32>(info, order);
    if(importType == "UINT32")
        return detail::readImageImpl<UInt32>(info, order);
    if(importType == "DOUBLE")
        return detail::readImageImpl<double>(info, order);
    vigra_fail("readImage(filename, import_type, order): import_type specifies an unknown pixel type.");
    return NumpyAnyArray();
}

// when export_type == "", writeImage() will figure out the best compromise
// between the input pixel type and the capabilities of the given file format
// (see negotiatePixelType())
template <class T>
void writeImage(NumpyArray<3, Multiband<T> > const & image,
                    const char * filename, 
                    python::object export_type,  
                    const char * compression = "")
{
    ImageExportInfo info(filename);
    
    if(python::extract<std::string>(export_type).check())
    {
        std::string type = python::extract<std::string>(export_type)();        
        if(type == "NBYTE")
        {
            info.setForcedRangeMapping(0.0, 0.0, 0.0, 255.0);
            info.setPixelType("UINT8");
        }
        else if(type != "" && type != "NATIVE")
        {
            info.setPixelType(type.c_str());
        }
    }
    else if(python::extract<NPY_TYPES>(export_type).check())
    {
        info.setPixelType(detail::numpyTypeIdToImpexString(python::extract<NPY_TYPES>(export_type)()).c_str());
    }
    else if(export_type)
        vigra_precondition(false, "writeImage(filename, export_type): export_type must be a string or a numpy dtype.");
        
    if(std::string(compression) == "RunLength")
        info.setCompression("RLE");
    else if(std::string(compression) != "")
        info.setCompression(compression);
    exportImage(srcImageRange(image), info);
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywriteImage, writeImage)

namespace detail {

template <class T>
NumpyAnyArray readVolumeImpl(VolumeImportInfo const & info, std::string order = "")
{
    typedef UnstridedArrayTag Stride;
    
    if(order == "")
        order = detail::defaultOrder();

    switch(info.numBands())
    {
      case 1:
      {
        NumpyArray<3, Singleband<T> > volume(info.shape(), order);
        importVolume(info, volume);
        return volume;
      }
      case 2:
      {
        NumpyArray<3, TinyVector<T, 2> > volume(info.shape(), order);
        importVolume(info, volume);
        return volume;
      }
      case 3:
      {
        NumpyArray<3, RGBValue<T> > volume(info.shape(), order);
        importVolume(info, volume);
        return volume;
      }
      case 4:
      {
        NumpyArray<3, TinyVector<T, 4> > volume(info.shape(), order);
        importVolume(info, volume);
        return volume;
      }
      //FIXME not yet supported
      /*default:
      {
        NumpyArray<4, Multiband<T> > volume(MultiArrayShape<4>::type(info.width(), info.height(), info.depth(), info.numBands()));
        importVolume(info, volume);
        return volume;
      }*/
      default:
      {
        NumpyArray<3, RGBValue<T> > volume(info.shape(), order);
        importVolume(info, volume);
        return volume;
      }
    }
}

} // namespace detail

NumpyAnyArray 
readVolume(const char * filename, python::object import_type, std::string order = "")
{
    VolumeImportInfo info(filename);
    std::string importType(info.getPixelType());

    if(python::extract<std::string>(import_type).check())
    {
        std::string type = python::extract<std::string>(import_type)();        
        if(type != "" && type != "NATIVE")
            importType = type;
    }
    else if(python::extract<NPY_TYPES>(import_type).check())
    {
        importType = detail::numpyTypeIdToImpexString(python::extract<NPY_TYPES>(import_type)());
    }
    else if(import_type)
        vigra_precondition(false, "readVolume(filename, import_type, order): import_type must be a string or a numpy dtype.");

    if(importType == "FLOAT")
        return detail::readVolumeImpl<float>(info, order);
    if(importType == "UINT8")
        return detail::readVolumeImpl<UInt8>(info, order);
    if(importType == "INT16")
        return detail::readVolumeImpl<Int16>(info, order);
    if(importType == "UINT16")
        return detail::readVolumeImpl<UInt16>(info, order);
    if(importType == "INT32")
        return detail::readVolumeImpl<Int32>(info, order);
    if(importType == "UINT32")
        return detail::readVolumeImpl<UInt32>(info, order);
    if(importType == "DOUBLE")
        return detail::readVolumeImpl<double>(info, order);
    vigra_fail("readVolume(filename, import_type, order): import_type specifies an unknown pixel type.");
    return NumpyAnyArray();
}

template <class T>
void writeVolume(NumpyArray<3, Multiband<T> > const & volume,
                    const char * filename_base, 
                    const char * filename_ext, 
                    python::object export_type,  
                    const char * compression = "")
{
    VolumeExportInfo info(filename_base, filename_ext);
    
    if(python::extract<std::string>(export_type).check())
    {
        std::string type = python::extract<std::string>(export_type)();        
        if(type == "NBYTE")
        {
            info.setForcedRangeMapping(0.0, 0.0, 0.0, 255.0);
            info.setPixelType("UINT8");
        }
        else if(type != "" && type != "NATIVE")
        {
            info.setPixelType(type.c_str());
        }
    }
    else if(python::extract<NPY_TYPES>(export_type).check())
    {
        info.setPixelType(detail::numpyTypeIdToImpexString(python::extract<NPY_TYPES>(export_type)()).c_str());
    }
    else if(export_type)
        vigra_precondition(false, "writeVolume(filename, export_type): export_type must be a string or a numpy dtype.");
        
    if(std::string(compression) == "RunLength")
        info.setCompression("RLE");
    else if(std::string(compression) != "")
        info.setCompression(compression);
    exportVolume(volume, info);
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywriteVolume, writeVolume)

#ifdef HasHDF5

namespace detail {
template <class T>
NumpyAnyArray readImageHDF5Impl(HDF5ImportInfo const & info, std::string order = "")
{
    typedef UnstridedArrayTag Stride;
    
    if(order == "")
        order = detail::defaultOrder();

    switch(info.numDimensions())
    {
      case 2:
      {
        NumpyArray<2, Singleband<T> > res(MultiArrayShape<2>::type(info.shapeOfDimension(0), 
                                                                    info.shapeOfDimension(1)),
                                           order);
        readHDF5(info, res);
        return res;
      }
      case 3:
      {
      
        if(info.shapeOfDimension(0) == 3)
        {
            NumpyArray<2, RGBValue<T> > res(MultiArrayShape<2>::type(info.shapeOfDimension(1), 
                                                                      info.shapeOfDimension(2)),
                                            order);
            readHDF5(info, res);
            return res;
        }
        else
        {
            NumpyArray<3, Multiband<T> > res(MultiArrayShape<3>::type(info.shapeOfDimension(1), 
                                                                       info.shapeOfDimension(2), 
                                                                       info.shapeOfDimension(0)),
                                             order);
            MultiArrayShape<3>::type permute(2,0,1);
            readHDF5(info, res.permuteDimensions(permute));
            return res;
        }
      }
      default:
      {
        vigra_precondition(false, "readImageFromHDF5(filename, datasetname, import_type): dataset has wrong number of dimensions (must be 2 or 3).");
        return NumpyAnyArray();
      }
    }
}

} // namespace detail

NumpyAnyArray 
readImageFromHDF5(const char * filePath, const char * pathInFile, 
                  python::object import_type, std::string order = "")
{
    HDF5ImportInfo info(filePath, pathInFile);
    std::string importType(info.getPixelType());
    
    if(python::extract<std::string>(import_type).check())
    {
        std::string type = python::extract<std::string>(import_type)();        
        if(type != "" && type != "NATIVE")
            importType = type;
    }
    else if(python::extract<NPY_TYPES>(import_type).check())
    {
        importType = detail::numpyTypeIdToImpexString(python::extract<NPY_TYPES>(import_type)());
    }
    else if(import_type)
        vigra_precondition(false, "readImageFromHDF5(filename, datasetname, import_type, order): import_type must be a string or a numpy dtype.");

    // FIXME: support all types, at least via a type cast at the end?
    if(importType == "FLOAT")
        return detail::readImageHDF5Impl<npy_float32>(info, order);
    if(importType == "UINT8")
        return detail::readImageHDF5Impl<npy_uint8>(info, order);
    if(importType == "INT16")
        return detail::readImageHDF5Impl<npy_int16>(info, order);
    if(importType == "UINT16")
        return detail::readImageHDF5Impl<npy_uint16>(info, order);
    if(importType == "INT32")
        return detail::readImageHDF5Impl<npy_int32>(info, order);
    if(importType == "UINT32")
        return detail::readImageHDF5Impl<npy_uint32>(info, order);
    if(importType == "DOUBLE")
        return detail::readImageHDF5Impl<npy_float64>(info, order);
    vigra_fail("readImageFromHDF5(filePath, pathInFile, import_type, order): import_type specifies an unknown pixel type.");
    return NumpyAnyArray();
}

template <class T>
void writeImageToHDF5(NumpyArray<3, Multiband<T> > const & image,
                    const char * filePath, 
                    const char * pathInFile, 
                    python::object export_type)  
{
    // write the data
    // if scalar image
    if(image.shape(2) == 1)
    {
        writeHDF5(filePath, pathInFile, image.bindOuter(0));
    }
    else
    {
        MultiArrayShape<3>::type permute(2,0,1);
        writeHDF5(filePath, pathInFile, image.permuteDimensions(permute));
    }
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywriteImageToHDF5, writeImageToHDF5)

namespace detail {
template <class T>
NumpyAnyArray 
readVolumeHDF5Impl(HDF5ImportInfo const & info, std::string order = "")
{
    typedef UnstridedArrayTag Stride;
    
    if(order == "")
        order = detail::defaultOrder();

    switch(info.numDimensions())
    {
      case 3:
      {
        NumpyArray<3, Singleband<T> > res(MultiArrayShape<3>::type(info.shapeOfDimension(0), 
                                                                    info.shapeOfDimension(1), 
                                                                    info.shapeOfDimension(2)),
                                           order);
        readHDF5(info, res);
        return res;
      }
      case 4:
      {
        NumpyArray<4, Multiband<T> > res(MultiArrayShape<4>::type(info.shapeOfDimension(1), 
                                                                   info.shapeOfDimension(2), 
                                                                   info.shapeOfDimension(3), 
                                                                   info.shapeOfDimension(0)),
                                         order);
        MultiArrayShape<4>::type permute(3, 0, 1, 2);
        readHDF5(info, res.permuteDimensions(permute));
        return res;
      }
      default:
      {
        vigra_precondition(false, "readVolumeFromHDF5(filename, datasetname, import_type): dataset has wrong number of dimensions (must be 3 or 4).");
        return NumpyAnyArray();
      }
    }
}
} // namespace detail

NumpyAnyArray 
readVolumeFromHDF5(const char * filePath, const char * pathInFile, 
                   python::object import_type, std::string order = "")
{
    HDF5ImportInfo info(filePath, pathInFile);
    std::string importType(info.getPixelType());
    
    if(python::extract<std::string>(import_type).check())
    {
        std::string type = python::extract<std::string>(import_type)();        
        if(type != "" && type != "NATIVE")
            importType = type;
    }
    else if(python::extract<NPY_TYPES>(import_type).check())
    {
        importType = detail::numpyTypeIdToImpexString(python::extract<NPY_TYPES>(import_type)());
    }
    else if(import_type)
        vigra_precondition(false, "readVolumeFromHDF5(filename, datasetname, import_type, order): import_type must be a string or a numpy dtype.");

    // FIXME: support all types, at least via a type cast at the end?
    if(importType == "FLOAT")
        return detail::readVolumeHDF5Impl<float>(info, order);
    if(importType == "UINT8")
        return detail::readVolumeHDF5Impl<UInt8>(info, order);
    if(importType == "INT16")
        return detail::readVolumeHDF5Impl<Int16>(info, order);
    if(importType == "UINT16")
        return detail::readVolumeHDF5Impl<UInt16>(info, order);
    if(importType == "INT32")
        return detail::readVolumeHDF5Impl<Int32>(info, order);
    if(importType == "UINT32")
        return detail::readVolumeHDF5Impl<UInt32>(info, order);
    if(importType == "DOUBLE")
        return detail::readVolumeHDF5Impl<double>(info, order);
    vigra_fail("readVolumeHDF5(filePath, pathInFile, import_type, order): import_type specifies an unknown pixel type.");
    return NumpyAnyArray();
}

template <class T>
void writeVolumeToHDF5(NumpyArray<4, Multiband<T> > const & volume,
                    const char * filePath, 
                    const char * pathInFile, 
                    python::object export_type) 
{
    // write the data
    // if scalar volume
    if(volume.shape(3) == 1)
    {
        writeHDF5(filePath, pathInFile, volume.bindOuter(0));
    }
    else
    {
        MultiArrayShape<4>::type permute(3,0,1,2);
        writeHDF5(filePath, pathInFile, volume.permuteDimensions(permute));
    }
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pywriteVolumeToHDF5, writeVolumeToHDF5)

#endif // HasHDF5

void defineImpexFunctions()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);
    
    // FIXME: add an order parameter to the import functions
    def("readVolume", &readVolume, (arg("filename"), arg("dtype") = "FLOAT", arg("order") = ""),
        "Read a 3D volume from a directory::\n"
        "\n"
        "   readVolume(filename, dtype='FLOAT', order='') -> Volume\n"
        "\n"
        "If the volume is stored in a by-slice manner (e.g. one image per\n"
        "slice), the 'filename' can refer to an arbitrary image from the set.\n"
        "readVolume() then assumes that the slices are enumerated like::\n"
        "\n"
        "   name_base+[0-9]+name_ext\n"
        "\n"
        "where name_base, the index, and name_ext\n"
        "are determined automatically. All slice files with the same name base\n"
        "and extension are considered part of the same volume. Slice numbers\n"
        "must be non-negative, but can otherwise start anywhere and need not\n"
        "be successive. Slices will be read in ascending numerical (not\n"
        "lexicographic) order. All slices must have the same size.\n"
        "\n"
        "Otherwise, readVolume() will try to read 'filename' as an info text\n"
        "file with the following key-value pairs::\n"
        "\n"
        "    name = [short descriptive name of the volume] (optional)\n"
        "    filename = [absolute or relative path to raw voxel data file] (required)\n"
        "    gradfile =  [abs. or rel. path to gradient data file] (currently ignored)\n"
        "    description =  [arbitrary description of the data set] (optional)\n"
        "    width = [positive integer] (required)\n"
        "    height = [positive integer] (required)\n"
        "    depth = [positive integer] (required)\n"
        "    datatype = [UNSIGNED_CHAR | UNSIGNED_BYTE] (default: UNSIGNED_CHAR)\n"
        "\n"
        "Lines starting with # are ignored.\n"
        "When import_type is 'UINT8', 'INT16', 'UINT16', 'INT32', 'UINT32',\n"
        "'FLOAT', 'DOUBLE', or one of the corresponding numpy dtypes (numpy.uint8\n"
        "etc.), the returned volume will have the requested pixel type. \n"
        "\n"
        "The order parameter determines the axis ordering of the resulting array\n"
        "(allowed values: 'C', 'F', 'V'). When order == '' (the default), " 
        "vigra.VigraArray.defaultOrder is used.\n"
        "\n"
        "For details see the help for :func:`readImage`.\n");
        
    multidef("writeVolume", pywriteVolume<Int8, UInt64, Int64, UInt16, Int16, UInt32, Int32, double, float, UInt8>(), 
       (arg("volume"), arg("filename_base"), arg("filename_ext"), arg("dtype") = "", arg("compression") = ""),
       "Wrtie a volume as a sequence of images::\n\n"
       "   writeVolume(volume, filename_base, filename_ext, dtype = '', compression = '')\n\n"
       "The resulting image sequence will be enumerated in the form::\n\n"
       "    filename_base+[0-9]+filename_ext\n\n"
       "Parameters 'dtype' and 'compression' will be handled as in :func:`writeImage`.\n\n");
    
    def("readImage", &readImage, (arg("filename"), arg("dtype") = "FLOAT", arg("order") = ""),
        "Read an image from a file::\n"
        "\n"
        "   readImage(filename, dtype='FLOAT', order='') -> Image\n"
        "\n"
        "When import_type is 'UINT8', 'INT16', 'UINT16', 'INT32', 'UINT32',\n"
        "'FLOAT', 'DOUBLE', or one of the corresponding numpy dtypes (numpy.uint8\n"
        "etc.), the returned image will have the requested pixel type. If\n"
        "dtype is 'NATIVE' or '' (empty string), the image is imported with\n"
        "the original type of the data in the file. Caution: If the requested \n"
        "dtype is smaller than the original type in the file, values will be\n"
        "clipped at the bounds of the representable range, which may not be the\n"
        "desired behavior.\n"
        "\n"
        "The order parameter determines the axis ordering of the resulting array\n"
        "(allowed values: 'C', 'F', 'V'). When order == '' (the default), " 
        "vigra.VigraArray.defaultOrder is used.\n"
        "\n"
        "Supported file formats are listed by the function vigra.impexListFormats().\n"
        "When filename does not refer to a recognized image file format, an\n"
        "exception is raised. The file can be checked beforehand with the function\n"
        ":func:`isImage`(filename).\n");
        
    multidef("writeImage", pywriteImage<Int8, UInt64, Int64, UInt16, Int16, UInt32, Int32, double, float, UInt8>(), 
       (arg("image"), arg("filename"), arg("dtype") = "", arg("compression") = ""),
        "Save an image to a file::\n"
        "\n"
        "   writeImage(image, filename, dtype = '', compression = '')\n"
        "\n"
        "Parameters:\n\n"
        " image:\n"
        "    the image to be saved\n"
        " filename:\n"
        "    the file name to save to. The file type will be deduced\n"
        "    from the file name extension (see vigra.impexListExtensions()\n"
        "    for a list of supported extensions).\n"
        " dtype:\n"
        "    the pixel type written to the file. Possible values:\n\n"
        "     '' or 'NATIVE':\n"
        "        save with original pixel type, or convert automatically\n"
        "        when this type is unsupported by the target file format\n"
        "     'UINT8', 'INT16', 'UINT16', 'INT32', 'UINT32', 'FLOAT', 'DOUBLE':\n"
        "        save as specified, or raise exception when this type is not \n"
        "        supported by the target file format (see list below)\n"
        "     'NBYTE':\n"
        "        normalize to range 0...255 and then save as 'UINT8'\n"
        "     numpy.uint8, numpy.int16 etc.:\n"
        "        behaves like the corresponding string argument\n\n" 
        " compression:\n"
        "     how to compress the data (ignored when compression type is unsupported \n"
        "     by the file format). Possible values:\n\n"
        "     '' or not given:\n"
        "        save with the native compression of the target file format\n"
        "     'RLE', 'RunLength':\n"
        "        use run length encoding (native in BMP, supported by TIFF)\n"
        "     'DEFLATE':\n"
        "        use deflate encoding (only supported by TIFF)\n"
        "     'LZW':\n"
        "        use LZW algorithm (only supported by TIFF with LZW enabled)\n"
        "     'ASCII':\n"
        "        write as ASCII rather than binary file (only supported by PNM)\n"
        "     '1' ... '100':\n"
        "        use this JPEG compression level (only supported by JPEG and TIFF)\n\n"
        "Supported file formats are listed by the function vigra.impexListFormats().\n"
        "The different file formats support the following pixel types:\n\n"
        "   BMP:\n"
        "       Microsoft Windows bitmap image file (pixel type: UINT8 as gray and RGB).\n" 
        "   GIF:\n"
        "       CompuServe graphics interchange format; 8-bit color\n"
        "       (pixel type: UINT8 as gray and RGB).\n" 
        "   JPEG:\n"
        "       Joint Photographic Experts Group JFIF format; compressed 24-bit color\n"
        "       (pixel types: UINT8 as gray and RGB). (only available if libjpeg is installed)\n" 
        "   PNG:\n"
        "       Portable Network Graphic (pixel types: UINT8 and UINT16 with\n"
        "       up to 4 channels). (only available if libpng is installed)\n"
        "   PBM:\n"
        "       Portable bitmap format (black and white).\n" 
        "   PGM:\n"
        "       Portable graymap format (pixel types: UINT8, INT16, INT32 as gray scale)).\n" 
        "   PNM:\n"
        "       Portable anymap (pixel types: UINT8, INT16, INT32, gray and RGB)\n" 
        "   PPM:\n"
        "       Portable pixmap format (pixel types: UINT8, INT16, INT32 as RGB)\n"
        "   SUN:\n"
        "       SUN Rasterfile (pixel types: UINT8 as gray and RGB).\n" 
        "   TIFF:\n"
        "       Tagged Image File Format (pixel types: UINT8, INT16, INT32, FLOAT, DOUBLE\n"
        "       with up to 4 channels). (only available if libtiff is installed.)\n"
        "   VIFF:\n"
        "       Khoros Visualization image file (pixel types: UINT8, INT16\n"
        "       INT32, FLOAT, DOUBLE with arbitrary many channels).\n\n");
        
    def("listFormats", &impexListFormats, 
        "Ask for the image file formats that vigra.impex understands::\n\n"
        "    listFormats() -> string\n\n"
        "This function returns a string containing the supported image file "
        "formats for reading and writing with the functions :func:`readImage` and "
        ":func:`writeImage`.\n");
        
    def("listExtensions", &impexListExtensions,
        "Ask for the image file extensions that vigra.impex understands::\n\n"
        "    listExtensions() -> string\n\n"
        "This function returns a string containing the supported image file "
        "extensions for reading and writing with the functions :func:`readImage` and "
        ":func:`writeImage`.\n");
        
    def("isImage", &isImage, args("filename"),
        "Check whether the given file name contains image data::\n\n"
        "   isImage(filename) -> bool\n\n"
        "This function tests whether a file has a supported image format. "
        "It checks the first few bytes of the file and compares them with "
        "the \"magic strings\" of each recognized image format. If the "
        "image format is supported it returns True otherwise False.\n");

#ifdef HasHDF5
    def("readImageFromHDF5", &readImageFromHDF5, 
         (arg("filepath"), arg("pathInFile"), arg("dtype") = "FLOAT", arg("order") = ""),
        "Read an image from a HDF5 file::\n"
        "\n"
        "   readImageFromHDF5(filepath, pathInFile, dtype='FLOAT', order='') -> image\n"
        "\n"
        "The order parameter determines the axis ordering of the resulting array\n"
        "(allowed values: 'C', 'F', 'V'). When order == '' (the default), " 
        "vigra.VigraArray.defaultOrder is used.\n"
        "\n"
        "If the file contains 3-dimensional data, the innermost (last)\n"
        "index is interpreted as a channel dimension.\n\n");
        
    def("readVolumeFromHDF5", &readVolumeFromHDF5, 
         (arg("filepath"), arg("pathInFile"), arg("dtype") = "FLOAT", arg("order") = ""),
        "Read a volume from a file::\n"
        "\n"
        "   readVolumeFromHDF5(filepath, pathInFile, dtype='FLOAT', order='') -> volume\n"
        "\n"
        "The order parameter determines the axis ordering of the resulting array\n"
        "(allowed values: 'C', 'F', 'V'). When order == '' (the default), " 
        "vigra.VigraArray.defaultOrder is used.\n"
        "\n"
        "If the file contains 4-dimensional data, the innermost (last)\n"
        "index is interpreted as a channel dimension.\n\n");
        
    multidef("writeImageToHDF5", pywriteImageToHDF5<Int8, UInt64, Int64, UInt16, Int16, UInt32, Int32, double, float, UInt8>(), 
       (arg("image"), arg("filepath"), arg("pathInFile"), arg("dtype") = ""),
        "Save an image to an HDF5 file::\n"
        "\n"
        "   writeImageToHDF5(image, filepath, pathInFile, dtype='')\n"
        "\n"
        "Argument 'dtype' is currently ignored.\n"
        "The resulting HDF5 dataset should be identical to the one created by\n"
        "the Python module `h5py <http://h5py.alfven.org/>`_ as follows::\n\n"
        "   import vigra, h5py\n"
        "   image = vigra.readImage(filename)\n"
        "   vigra.writeImageToHDF5(image, 'vigra_export.h5', 'MyImage')\n"
        "   h5py_file = h5py.File('h5py_export.h5', 'w')\n"
        "   h5py_file.create_dataset('MyImage', data=image.transposeToNumpyOrder())\n"
        "   h5py_file.close()\n"
        "\n"
        "(note the axes transposition which accounts for the VIGRA indexing convention).\n\n");
        
    multidef("writeVolumeToHDF5", pywriteVolumeToHDF5<Int8, UInt64, Int64, UInt16, Int16, UInt32, Int32, double, float, UInt8>(), 
       (arg("volume"), arg("filepath"), arg("pathInFile"), arg("dtype") = ""),
        "Save a volume to an HDF5 file::\n"
        "\n"
        "   writeVolumeToHDF5(volume, filepath, pathInFile, dtype='')\n"
        "\n"
        "Argument 'dtype' is currently ignored.\n"
        "The resulting HDF5 dataset should be identical to the one created by\n"
        "the Python module `h5py <http://h5py.alfven.org/>`_ as follows::\n\n"
        "   import vigra, h5py\n"
        "   volume = vigra.readVolume(filename)\n"
        "   vigra.writeVolumeToHDF5(volume, 'vigra_export.h5', 'MyVolume')\n"
        "   h5py_file = h5py.File('h5py_export.h5', 'w')\n"
        "   h5py_file.create_dataset('MyVolume', data=volume.transposeToNumpyOrder())\n"
        "   h5py_file.close()\n"
        "\n"
        "(note the axes transposition which accounts for the VIGRA indexing convention).\n\n");

#endif // HasHDF5
}

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(impex)
{
    import_vigranumpy();
    defineImpexFunctions();
}
