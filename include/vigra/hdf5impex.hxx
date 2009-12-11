/************************************************************************/
/*                                                                      */
/*       Copyright 2009 by Michael Hanselmann and Ullrich Koethe        */
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

#ifndef VIGRA_HDF5IMPEX_HXX
#define VIGRA_HDF5IMPEX_HXX

#include <string>
#include <hdf5.h>

#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
# ifndef H5Gopen
#   define H5Gopen(a, b, c) H5Gopen(a, b)
# endif
# ifndef H5Gcreate
#  define H5Gcreate(a, b, c, d, e) H5Gcreate(a, b, 1)
# endif
# ifndef H5Dopen
#  define H5Dopen(a, b, c) H5Dopen(a, b)
# endif
# ifndef H5Dcreate
#  define H5Dcreate(a, b, c, d, e, f, g) H5Dcreate(a, b, c, d, f)
# endif
# include <H5LT.h>
#else
# include <hdf5_hl.h>
#endif

#include "impex.hxx"
#include "multi_array.hxx"
#include "multi_impex.hxx"

namespace vigra {

class HDF5Handle
{
public:
    typedef herr_t (*Destructor)(hid_t);
    
private:
    hid_t handle_;
    Destructor destructor_;
    
public:

    HDF5Handle()
    : handle_( 0 ),
      destructor_(0)
    {}

    HDF5Handle(hid_t h, Destructor destructor, const char * error_message)
    : handle_( h ),
      destructor_(destructor)
    {
        if(handle_ < 0)
            vigra_fail(error_message);
    }

    HDF5Handle(HDF5Handle const & h)
    : handle_( h.handle_ ),
      destructor_(h.destructor_)
    {
        const_cast<HDF5Handle &>(h).handle_ = 0;
    }
    
    HDF5Handle & operator=(HDF5Handle const & h)
    {
        if(h.handle_ != handle_)
        {
            close();
            handle_ = h.handle_;
            destructor_ = h.destructor_;
            const_cast<HDF5Handle &>(h).handle_ = 0;
        }
        return *this;
    }

    ~HDF5Handle()
    {
        close();
    }
    
    herr_t close()
    {
        herr_t res = 1;
        if(handle_ && destructor_)
            res = (*destructor_)(handle_);
        handle_ = 0;
        return res;
    }

    hid_t get() const
    {
        return handle_;
    }

    operator hid_t() const
    {
        return handle_;
    }

    bool operator==(HDF5Handle const & h) const
    {
        return handle_ == h.handle_;
    }

    bool operator==(hid_t h) const
    {
        return handle_ == h;
    }

    bool operator!=(HDF5Handle const & h) const
    {
        return handle_ != h.handle_;
    }

    bool operator!=(hid_t h) const
    {
        return handle_ != h;
    }
};


/********************************************************/
/*                                                      */
/*                   HDF5ImportInfo                     */
/*                                                      */
/********************************************************/

/**
Namespace: vigra
**/
class HDF5ImportInfo
{
  public:
    enum PixelType { UINT8, UINT16, UINT32, UINT64, 
	   				 INT8, INT16, INT32, INT64,
					 FLOAT, DOUBLE };

        /** Construct HDF5ImageImportInfo object.

            The dataset in the given HDF5 file is accessed and the properties 
            are set accordingly.
         **/
    VIGRA_EXPORT HDF5ImportInfo( const char* filePath, const char* pathInFile );

    VIGRA_EXPORT ~HDF5ImportInfo();

    VIGRA_EXPORT const std::string& getFilePath() const;

    VIGRA_EXPORT const std::string& getPathInFile() const;

    VIGRA_EXPORT const hid_t getH5FileHandle() const;

    VIGRA_EXPORT const hid_t getDatasetHandle() const;

    VIGRA_EXPORT MultiArrayIndex numDimensions() const;

    VIGRA_EXPORT ArrayVector<hsize_t> const & shape() const
    {
        return m_dims;
    }

    VIGRA_EXPORT MultiArrayIndex shapeOfDimension(const int dim) const;

        /** Query the pixel type of the dataset.

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

        /** Query the pixel type of the dataset.

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

  private:
    HDF5Handle m_file_handle, m_dataset_handle;
    std::string m_filename, m_path, m_pixeltype;
    hssize_t m_dimensions;
    ArrayVector<hsize_t> m_dims;
};

namespace detail {

template<class type>
inline hid_t getH5DataType()
{
	std::runtime_error("getH5DataType(): invalid type");
	return 0;
}

#define VIGRA_H5_DATATYPE(type, h5type) \
template<> \
inline hid_t getH5DataType<type>() \
{ return h5type;}
VIGRA_H5_DATATYPE(char, H5T_NATIVE_CHAR)
VIGRA_H5_DATATYPE(Int8, H5T_NATIVE_INT8)
VIGRA_H5_DATATYPE(Int16, H5T_NATIVE_INT16)
VIGRA_H5_DATATYPE(Int32, H5T_NATIVE_INT32)
VIGRA_H5_DATATYPE(Int64, H5T_NATIVE_INT64)
VIGRA_H5_DATATYPE(UInt8, H5T_NATIVE_UINT8)
VIGRA_H5_DATATYPE(UInt16, H5T_NATIVE_UINT16)
VIGRA_H5_DATATYPE(UInt32, H5T_NATIVE_UINT32)
VIGRA_H5_DATATYPE(UInt64, H5T_NATIVE_UINT64)
VIGRA_H5_DATATYPE(float, H5T_NATIVE_FLOAT)
VIGRA_H5_DATATYPE(double, H5T_NATIVE_DOUBLE)
VIGRA_H5_DATATYPE(long double, H5T_NATIVE_LDOUBLE)

#undef VIGRA_H5_DATATYPE

} // namespace detail

namespace detail {

template <class Shape>
inline void
selectHyperslabs(HDF5Handle & mid1, HDF5Handle & mid2, Shape const & shape, int & counter, const int elements)
{
    // select hyperslab in HDF5 file
    hsize_t shapeHDF5[2];
    shapeHDF5[0] = 1;
    shapeHDF5[1] = elements;
    hsize_t startHDF5[2];
    startHDF5[0] = 0;
    startHDF5[1] = counter * shape[0];    
    hsize_t strideHDF5[2];
    strideHDF5[0] = 1;
    strideHDF5[1] = 1;                        
    hsize_t countHDF5[2];
    countHDF5[0] = 1;
    countHDF5[1] = shape[0];
    hsize_t blockHDF5[2];
    blockHDF5[0] = 1;
    blockHDF5[1] = 1;
    mid1 = HDF5Handle(H5Screate_simple(2, shapeHDF5, NULL),
                      &H5Sclose, "unable to create hyperslabs."); 
    H5Sselect_hyperslab(mid1, H5S_SELECT_SET, startHDF5, strideHDF5, countHDF5, blockHDF5);
    // select hyperslab in input data object
    hsize_t shapeData[2];
    shapeData[0] = 1;
    shapeData[1] = shape[0];
    hsize_t startData[2];
    startData[0] = 0;
    startData[1] = 0;
    hsize_t strideData[2];
    strideData[0] = 1;
    strideData[1] = 1;
    hsize_t countData[2];
    countData[0] = 1;
    countData[1] = shape[0];
    hsize_t blockData[2];
    blockData[0] = 1;
    blockData[1] = 1;
    mid2 = HDF5Handle(H5Screate_simple(2, shapeData, NULL),
                      &H5Sclose, "unable to create hyperslabs."); 
    H5Sselect_hyperslab(mid2, H5S_SELECT_SET, startData, strideData, countData, blockData);
}

template <class DestIterator, class Shape, class T>
inline void
readHDF5Impl(DestIterator d, Shape const & shape, hid_t dataset_id, ArrayVector<T> & buffer, int & counter, const int elements, MetaInt<0>)
{
    HDF5Handle mid1, mid2;

    // select hyperslabs
    selectHyperslabs(mid1, mid2, shape, counter, elements);

    // read from hdf5
    H5Dread(dataset_id, detail::getH5DataType<T>(), mid2, mid1, H5P_DEFAULT, buffer.data());

    // increase counter
    counter++;

    DestIterator dend = d + shape[0];
    int k = 0;
    for(; d < dend; ++d, k++)
    {
        *d = buffer[k];
        //std::cout << buffer[k] << " ";
    }

}

template <class DestIterator, class Shape, class T, int N>
void
readHDF5Impl(DestIterator d, Shape const & shape, hid_t dataset_id, ArrayVector<T> & buffer, int & counter, const int elements, MetaInt<N>)
{
    DestIterator dend = d + shape[N];
    for(; d < dend; ++d)
    {
        readHDF5Impl(d.begin(), shape, dataset_id, buffer, counter, elements, MetaInt<N-1>());
    }
}

} // namespace detail

template<unsigned int N, class T, class Tag>
void loadFromHDF5File(const HDF5ImportInfo &info, MultiArrayView<N, T, Tag> array, const bool rowMajorOrder = false) 
{
    //std::cout << N << " vs. " << info.numDimensions() << std::endl;
    vigra_precondition((N == info.numDimensions()),// || (N == 1 + info.numDimensions())),
        "loadFromHDF5File(): Array dimension disagrees with HDF5ImportInfo.numDimensions().");

    typename MultiArrayShape<N>::type shape;
    for(unsigned int k=0; k<N; ++k)
        shape[k] = (MultiArrayIndex)info.shapeOfDimension(k);

    vigra_precondition(shape == array.shape(), 
         "loadFromHDF5File(): Array shape disagrees with HDF5ImportInfo.");

    //Get the data
    int counter = 0;
    int elements = 1;
    for(int i=0;i<N;++i)
        elements *= shape[i];
    if(rowMajorOrder)
    {
        ArrayVector<T> buffer(shape[0]);
        detail::readHDF5Impl(array.traverser_begin(), shape, info.getDatasetHandle(), buffer, counter, elements, vigra::MetaInt<N-1>());
    } else {
        /*
        MultiArrayView<N, T, StridedArrayTag> arrayTransposed = array.permuteStridesDescending();
        ArrayVector<T> buffer(arrayTransposed.shape(0));
        detail::readHDF5Impl(arrayTransposed.traverser_begin(), arrayTransposed.shape(), info.getDatasetHandle(), buffer, counter, elements, vigra::MetaInt<N-1>());
        */
        vigra::TinyVector<int,N> strideNew;
        vigra::TinyVector<int,N> shapeNew;
        for(unsigned int k=0; k<N; ++k)
        {
            //std::cout << "StrideOld[" << k << "]=" << array.stride(k) << std::endl;
            strideNew[k] = array.stride(N-1-k);
            shapeNew[k] = array.shape(N-1-k);
            //std::cout << "StrideNew[" << k << "]=" << strideNew[k] << std::endl;
            //std::cout << "ShapeNew[" << k << "]=" << shapeNew[k] << std::endl;
        }
        MultiArrayView<N, T, StridedArrayTag> arrayNew (shapeNew, strideNew, array.data());
        ArrayVector<T> buffer(arrayNew.shape(0));
        detail::readHDF5Impl(arrayNew.traverser_begin(), arrayNew.shape(), info.getDatasetHandle(), buffer, counter, elements, vigra::MetaInt<N-1>());
    }

    /*vigra_postcondition(H5Dread(info.getDatasetHandle(), detail::getH5DataType<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data()) >= 0,
           "loadFromHDF5File(): Unable to transfer data.");
    */
}

inline hid_t createAllGroups(hid_t parent, std::string group_name)
{
    //std::cout << group_name << std::endl;
    std::string::size_type begin = 0, end = group_name.find('/');
    int a = 0;
    while (end != std::string::npos)
    {
        std::string group(group_name.begin()+begin, group_name.begin()+end);
        //std::cout << "createAllGroups(1): " << group.c_str() << std::endl;
        // FIXME: also handle relative paths correctly
        
        hid_t to_close = parent;

        if(H5LTfind_dataset(parent, group.c_str()) == 0)
        {
            //std::cout << "180 release, exists=false" << std::endl;
            //std::cout << "parent_a=" << parent << std::endl;
            parent = H5Gcreate(parent, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            //std::cout << "parent_b=" << parent << ", " << group.c_str() << std::endl;
        } else {
            //std::cout << "180 release, exists=true" << std::endl;
            //std::cout << "parent_a=" << parent << std::endl;
            parent = H5Gopen(parent, group.c_str(), H5P_DEFAULT);
            //std::cout << "parent_b=" << parent << ", " << group.c_str() << std::endl;
        }
        if(a > 0)
            H5Gclose(to_close); 
        ++a;
        if(parent < 0)
            return parent;
        begin = end + 1;
        end = group_name.find('/', begin);
    }
    std::string group(group_name.begin()+begin, group_name.end());
    //std::cout << "createAllGroups(2): " << group.c_str() << std::endl;
    if(H5LTfind_dataset(parent, group.c_str()) == 0)
    {
        //std::cout << "parent_a=" << parent << std::endl;
        parent = H5Gcreate(parent, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //std::cout << "parent_b=" << parent << ", " << group.c_str() << std::endl;
        return parent;
    } else {
        //std::cout << "parent_a=" << parent << std::endl;
        parent = H5Gopen(parent, group.c_str(), H5P_DEFAULT);
        //std::cout << "parent_b=" << parent << ", " << group.c_str() << std::endl;
        return parent;
    }
}

namespace detail {

template <class DestIterator, class Shape, class T>
inline void
writeHDF5Impl(DestIterator d, Shape const & shape, hid_t file_id, hid_t dataset_id, ArrayVector<T> & buffer, int & counter, const int elements, MetaInt<0>)
{
    DestIterator dend = d + (typename DestIterator::difference_type)shape[0];
    int k = 0;
    for(; d < dend; ++d, k++)
    {
        buffer[k] = *d;
        //std::cout << buffer[k] << " ";
    }
    HDF5Handle mid1, mid2;

    // select hyperslabs
    selectHyperslabs(mid1, mid2, shape, counter, elements);

    // write to hdf5
    H5Dwrite(dataset_id, detail::getH5DataType<T>(), mid2, mid1, H5P_DEFAULT, buffer.data());
    // increase counter
    counter++;
}

template <class DestIterator, class Shape, class T, int N>
void
writeHDF5Impl(DestIterator d, Shape const & shape, hid_t file_id, hid_t dataset_id, ArrayVector<T> & buffer, int & counter, const int elements, MetaInt<N>)
{
    DestIterator dend = d + (typename DestIterator::difference_type)shape[N];
    for(; d < dend; ++d)
    {
        writeHDF5Impl(d.begin(), shape, file_id, dataset_id, buffer, counter, elements, MetaInt<N-1>());
    }
}

} // namespace detail


template<unsigned int N, class T, class Tag>
void writeToHDF5File(const char* filePath, const char* pathInFile, const MultiArrayView<N, T, Tag> & array, const bool rowMajorOrder = false)
{
    /*
    std::cout << "Values (0,0), (0,1), (1,0): " << array(0,0) << " " << array(0,1) << " " << array(1,0) << " " << std::endl;
    std::cout << "Shape  (0), (1): " << array.shape(0) << " " << array.shape(1) << " " << std::endl;
    std::cout << "Stride (0), (1): " << array.stride(0) << " " << array.stride(1) << " " << std::endl;
    */

    // check if file already exists
    HDF5Handle file_id;
    FILE * pFile;
    pFile = fopen ( filePath, "r" );
    if ( pFile == NULL )
    {
        // if not, create the file
        file_id = HDF5Handle(H5Fcreate(filePath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT),
                             &H5Fclose, "writeToHDF5File(): output file could not be created.");
    } else {
        fclose( pFile );
        pFile = 0;
        // open the file
        file_id = HDF5Handle(H5Fopen(filePath, H5F_ACC_RDWR, H5P_DEFAULT),
                             &H5Fclose, "writeToHDF5File(): unable to open output file.");
    }

    std::string path_name(pathInFile), group_name, data_set_name, message;
    std::string::size_type delimiter = path_name.rfind('/');
    
    if(delimiter == std::string::npos)
    {
        data_set_name = path_name;
    }
    else
    {
        group_name = std::string(path_name.begin(), path_name.begin()+delimiter);
        data_set_name = std::string(path_name.begin()+delimiter+1, path_name.end());
    }

    HDF5Handle group;
    if(group_name != "")
    {
        group = HDF5Handle(createAllGroups(file_id, group_name), &H5Gclose, "writeToHDF5File(): Unable to create and open group.");
    }
    else
    {
        group = HDF5Handle(file_id, 0, "");
    }

    hsize_t shape[N];
    for(unsigned int k=0; k<N; ++k)
        shape[k] = array.shape(k);

    // create dataspace
    HDF5Handle dataspace_handle(H5Screate_simple(N, shape, NULL),
                                &H5Sclose, "writeToHDF5File(): unable to create dataspace.");

    // delete existing data and create new dataset
    if(H5LTfind_dataset(group, data_set_name.c_str()))
    {
        //std::cout << "dataset already exists" << std::endl;
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
		if(H5Gunlink(group, data_set_name.c_str()) < 0)
        {
            vigra_postcondition(false, "writeToHDF5File(): Unable to delete existing data.");
        }
#else
		if(H5Ldelete(group, data_set_name.c_str(), H5P_DEFAULT ) < 0)
        {
            vigra_postcondition(false, "writeToHDF5File(): Unable to delete existing data.");
        }
#endif
    } /*else {
        std::cout << "dataset does not exist so far" << std::endl;
    }*/
   	detail::getH5DataType<char>();
   	detail::getH5DataType<Int8>();
    HDF5Handle dataset_handle(H5Dcreate(group, data_set_name.c_str(), detail::getH5DataType<T>(), dataspace_handle, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                              &H5Dclose, "writeToHDF5File(): unable to create dataset.");
    
    // Write the data to the HDF5 dataset
    //dataset.write( array.data(), GetH5DataType<T>() ); // old version without support for strided arrays
    int elements = 1;
    for(int i=0;i<int(N);++i)
        elements *= (int)shape[i];
    int counter = 0;

    if(rowMajorOrder)
    {
        ArrayVector<T> buffer((int)shape[0]);
        detail::writeHDF5Impl(array.traverser_begin(), shape, file_id, dataset_handle, buffer, counter, elements, vigra::MetaInt<N-1>());
    } else {
        // for column major order we have to reverse the shape and strides before calling the write function
        vigra::TinyVector<int,N> strideNew;
        vigra::TinyVector<int,N> shapeNew;
        for(unsigned int k=0; k<N; ++k)
        {
            strideNew[k] = array.stride(N-1-k);
            shapeNew[k] = array.shape(N-1-k);
            //std::cout << "StrideNew[" << k << "]=" << strideNew[k] << std::endl;
            //std::cout << "ShapeNew[" << k << "]=" << shapeNew[k] << std::endl;
        }
        MultiArrayView<N, T, StridedArrayTag> arrayNew (shapeNew, strideNew, array.data());
        ArrayVector<T> buffer((int)arrayNew.shape(0));
        detail::writeHDF5Impl(arrayNew.traverser_begin(), arrayNew.shape(), file_id, dataset_handle, buffer, counter, elements, vigra::MetaInt<N-1>());
    }

    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
}


} // namespace vigra

#endif // VIGRA_HDF5IMPEX_HXX
