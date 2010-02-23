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
# ifndef H5Acreate
#  define H5Acreate(a, b, c, d, e, f) H5Acreate(a, b, c, d, e)
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
selectHyperslabs(HDF5Handle & mid1, HDF5Handle & mid2, Shape const & shape, int & counter, const int elements, const int numBandsOfType)
{
    // select hyperslab in HDF5 file
    hsize_t shapeHDF5[2];
    shapeHDF5[0] = 1;
    shapeHDF5[1] = elements;
    hsize_t startHDF5[2];
    startHDF5[0] = 0;
    startHDF5[1] = counter * numBandsOfType * shape[0]; // we have to reserve space for the pixel type channel(s)
    hsize_t strideHDF5[2];
    strideHDF5[0] = 1;
    strideHDF5[1] = 1;                        
    hsize_t countHDF5[2];
    countHDF5[0] = 1;
    countHDF5[1] = numBandsOfType * shape[0];
    hsize_t blockHDF5[2];
    blockHDF5[0] = 1;
    blockHDF5[1] = 1;
    mid1 = HDF5Handle(H5Screate_simple(2, shapeHDF5, NULL),
                      &H5Sclose, "unable to create hyperslabs."); 
    H5Sselect_hyperslab(mid1, H5S_SELECT_SET, startHDF5, strideHDF5, countHDF5, blockHDF5);
    // select hyperslab in input data object
    hsize_t shapeData[2];
    shapeData[0] = 1;
    shapeData[1] = numBandsOfType * shape[0];
    hsize_t startData[2];
    startData[0] = 0;
    startData[1] = 0;
    hsize_t strideData[2];
    strideData[0] = 1;
    strideData[1] = 1;
    hsize_t countData[2];
    countData[0] = 1;
    countData[1] = numBandsOfType * shape[0];
    hsize_t blockData[2];
    blockData[0] = 1;
    blockData[1] = 1;
    mid2 = HDF5Handle(H5Screate_simple(2, shapeData, NULL),
                      &H5Sclose, "unable to create hyperslabs."); 
    H5Sselect_hyperslab(mid2, H5S_SELECT_SET, startData, strideData, countData, blockData);
}

template <class DestIterator, class Shape, class T>
inline void
readHDF5Impl(DestIterator d, Shape const & shape, const hid_t dataset_id, const hid_t datatype, ArrayVector<T> & buffer, int & counter, const int elements, const int numBandsOfType, MetaInt<0>)
{
    HDF5Handle mid1, mid2;

    // select hyperslabs
    selectHyperslabs(mid1, mid2, shape, counter, elements, numBandsOfType);

    // read from hdf5
    H5Dread(dataset_id, datatype, mid2, mid1, H5P_DEFAULT, buffer.data());

    // increase counter
    counter++;


	//std::cout << "numBandsOfType: " << numBandsOfType << std::endl;
    DestIterator dend = d + shape[0];
    int k = 0;
    for(; d < dend; ++d, k++)
    {
        *d = buffer[k];
        //std::cout << buffer[k] << "| ";
    }

}

template <class DestIterator, class Shape, class T, int N>
void
readHDF5Impl(DestIterator d, Shape const & shape, const hid_t dataset_id, const hid_t datatype, ArrayVector<T> & buffer, int & counter, const int elements, const int numBandsOfType, MetaInt<N>)
{
    DestIterator dend = d + shape[N];
    for(; d < dend; ++d)
    {
        readHDF5Impl(d.begin(), shape, dataset_id, datatype, buffer, counter, elements, numBandsOfType, MetaInt<N-1>());
    }
}

} // namespace detail

// scalar and unstrided target multi array
template<unsigned int N, class T>
inline void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, T, UnstridedArrayTag> array) // scalar
{
	readHDF5(info, array, detail::getH5DataType<T>(), 1);
}

// non-scalar (TinyVector) and unstrided target multi array
template<unsigned int N, class T, int SIZE>
inline void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, TinyVector<T, SIZE>, UnstridedArrayTag> array)
{
	readHDF5(info, array, detail::getH5DataType<T>(), SIZE);
}

// non-scalar (RGBValue) and unstrided target multi array
template<unsigned int N, class T>
inline void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, RGBValue<T>, UnstridedArrayTag> array)
{
	readHDF5(info, array, detail::getH5DataType<T>(), 3);
}

// unstrided target multi array
template<unsigned int N, class T>
void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, T, UnstridedArrayTag> array, const hid_t datatype, const int numBandsOfType) 
{
	int offset = (numBandsOfType > 1);

	//std::cout << "offset: " << offset << ", N: " << N << ", dims: " << info.numDimensions() << std::endl;
	vigra_precondition(( (N + offset ) == info.numDimensions()), // the object in the HDF5 file may have one additional dimension which we then interpret as the pixel type bands
        "readHDF5(): Array dimension disagrees with HDF5ImportInfo.numDimensions().");

    typename MultiArrayShape<N>::type shape;
	for(int k=offset; k<info.numDimensions(); ++k) {
        shape[k-offset] = info.shapeOfDimension(k); 
	}

	vigra_precondition(shape == array.shape(), 
         "readHDF5(): Array shape disagrees with HDF5ImportInfo.");

	// simply read in the data as is
	H5Dread( info.getDatasetHandle(), datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data() ); // .data() possible since void pointer!
}

// scalar and strided target multi array
template<unsigned int N, class T>
inline void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, T, StridedArrayTag> array) // scalar
{
	readHDF5(info, array, detail::getH5DataType<T>(), 1);
}

// non-scalar (TinyVector) and strided target multi array
template<unsigned int N, class T, int SIZE>
inline void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, TinyVector<T, SIZE>, StridedArrayTag> array) 
{
	readHDF5(info, array, detail::getH5DataType<T>(), SIZE);
}

// non-scalar (RGBValue) and strided target multi array
template<unsigned int N, class T>
inline void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, RGBValue<T>, StridedArrayTag> array) 
{
	readHDF5(info, array, detail::getH5DataType<T>(), 3);
}

// strided target multi array
template<unsigned int N, class T>
void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, T, StridedArrayTag> array, const hid_t datatype, const int numBandsOfType)
{
	int offset = (numBandsOfType > 1);

	//std::cout << "offset: " << offset << ", N: " << N << ", dims: " << info.numDimensions() << std::endl;
	vigra_precondition(( (N + offset ) == info.numDimensions()), // the object in the HDF5 file may have one additional dimension which we then interpret as the pixel type bands
        "readHDF5(): Array dimension disagrees with HDF5ImportInfo.numDimensions().");

    typename MultiArrayShape<N>::type shape;
	for(int k=offset; k<info.numDimensions(); ++k) {
        shape[k-offset] = info.shapeOfDimension(k); 
	}

	vigra_precondition(shape == array.shape(), 
         "readHDF5(): Array shape disagrees with HDF5ImportInfo.");

    //Get the data
    int counter = 0;
    int elements = numBandsOfType;
    for(int i=0;i<N;++i)
        elements *= shape[i];
    ArrayVector<T> buffer(shape[0]);
    detail::readHDF5Impl(array.traverser_begin(), shape, info.getDatasetHandle(), datatype, buffer, counter, elements, numBandsOfType, vigra::MetaInt<N-1>());
}

inline hid_t openGroup(hid_t parent, std::string group_name)
{
    //std::cout << group_name << std::endl;
    size_t last_slash = group_name.find_last_of('/'); 
    if (last_slash == std::string::npos || last_slash != group_name.size() - 1)
        group_name = group_name + '/';
    std::string::size_type begin = 0, end = group_name.find('/');
    int ii =  0;
    while (end != std::string::npos)
    {
        std::string group(group_name.begin()+begin, group_name.begin()+end);
        hid_t prev_parent = parent; 
        parent = H5Gopen(prev_parent, group.c_str(), H5P_DEFAULT);

        if(ii != 0)     H5Gclose(prev_parent);
        if(parent < 0)  return parent;
        ++ii; 
        begin = end + 1;
        end = group_name.find('/', begin);
    }
    return parent; 
}

inline hid_t createGroup(hid_t parent, std::string group_name)
{
    if(group_name.size() == 0 ||*group_name.rbegin() != '/')
        group_name = group_name + '/';
    if(group_name == "/")
        return H5Gopen(parent, group_name.c_str(), H5P_DEFAULT);
    
    std::string::size_type begin = 0, end = group_name.find('/');
    int ii =  0;
    while (end != std::string::npos)
    {
        std::string group(group_name.begin()+begin, group_name.begin()+end);
        hid_t prev_parent = parent; 
        
        if(H5LTfind_dataset(parent, group.c_str()) == 0)
        {
            parent = H5Gcreate(prev_parent, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        } else {
            parent = H5Gopen(prev_parent, group.c_str(), H5P_DEFAULT);
        }

        if(ii != 0)     H5Gclose(prev_parent);
        if(parent < 0)  return parent;
        ++ii; 
        begin = end + 1;
        end = group_name.find('/', begin);
    }
    return parent; 
}

inline void deleteDataset(hid_t parent, std::string dataset_name)
{
    // delete existing data and create new dataset
    if(H5LTfind_dataset(parent, dataset_name.c_str()))
    {
        //std::cout << "dataset already exists" << std::endl;
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
		if(H5Gunlink(parent, dataset_name.c_str()) < 0)
        {
            vigra_postcondition(false, "writeToHDF5File(): Unable to delete existing data.");
        }
#else
		if(H5Ldelete(parent, dataset_name.c_str(), H5P_DEFAULT ) < 0)
        {
            vigra_postcondition(false, "createDataset(): Unable to delete existing data.");
        }
#endif
    } 
}

inline hid_t createFile(std::string filePath, bool append_ = true)
{
    FILE * pFile;
    pFile = fopen ( filePath.c_str(), "r" );
    hid_t file_id; 
    if ( pFile == NULL )
    {
        file_id = H5Fcreate(filePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    } 
    else if(append_)
    {
        fclose( pFile );
        file_id = H5Fopen(filePath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else
    {
        fclose(pFile);
        std::remove(filePath.c_str());
        file_id = H5Fcreate(filePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    return file_id; 
}

template<unsigned int N, class T, class Tag>
void createDataset(const char* filePath, const char* pathInFile, const MultiArrayView<N, T, Tag> & array, const hid_t datatype, const int numBandsOfType, HDF5Handle & file_handle, HDF5Handle & dataset_handle)
{
    std::string path_name(pathInFile), group_name, data_set_name, message;
    std::string::size_type delimiter = path_name.rfind('/');
    
    //create or open file
    file_handle = HDF5Handle(createFile(filePath), &H5Fclose, 
                       "createDataset(): unable to open output file.");

    // get the groupname and the filename
    if(delimiter == std::string::npos)
    {
        group_name    = "/";
        data_set_name = path_name;
    }
    else
    {
        group_name = std::string(path_name.begin(), path_name.begin()+delimiter);
        data_set_name = std::string(path_name.begin()+delimiter+1, path_name.end());
    }

    // create all groups
    HDF5Handle group(createGroup(file_handle, group_name), &H5Gclose, 
                     "createDataset(): Unable to create and open group. generic v");

    // delete the dataset if it already exists
    deleteDataset(group, data_set_name);

    // create dataspace
	// add an extra dimension in case that the data is non-scalar
	HDF5Handle dataspace_handle;
	if(numBandsOfType > 1) {
		// invert dimensions to guarantee c-order
		hsize_t shape_inv[N+1]; // one additional dimension for pixel type channel(s)
		for(unsigned int k=0; k<N; ++k) {
			shape_inv[N-1-k] = array.shape(k);  // the channels (eg of an RGB image) are represented by the first dimension (before inversion)
			//std::cout << shape_inv[N-k] << " (" << N << ")";
		}
		shape_inv[N] = numBandsOfType;

		// create dataspace
		dataspace_handle = HDF5Handle(H5Screate_simple(N+1, shape_inv, NULL),
									&H5Sclose, "createDataset(): unable to create dataspace for non-scalar data.");
	} else {
		// invert dimensions to guarantee c-order
		hsize_t shape_inv[N];
		for(unsigned int k=0; k<N; ++k)
			shape_inv[N-1-k] = array.shape(k);

		// create dataspace
		dataspace_handle = HDF5Handle(H5Screate_simple(N, shape_inv, NULL),
									&H5Sclose, "createDataset(): unable to create dataspace for scalar data.");
	}

    //alloc memory for dataset. 
    dataset_handle = HDF5Handle(H5Dcreate(group, 
                                        data_set_name.c_str(), 
                                        datatype, 
                                        dataspace_handle, 
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                              &H5Dclose, "createDataset(): unable to create dataset.");
}



namespace detail {

template <class DestIterator, class Shape, class T>
inline void
writeHDF5Impl(DestIterator d, Shape const & shape, const hid_t file_id, const hid_t dataset_id, const hid_t datatype, ArrayVector<T> & buffer, int & counter, const int elements, const int numBandsOfType, MetaInt<0>)
{
    DestIterator dend = d + (typename DestIterator::difference_type)shape[0];
    int k = 0;
	//std::cout << "new:" << std::endl;
	for(; d < dend; ++d, k++)
    {
        buffer[k] = *d; 
        //std::cout << buffer[k] << " ";
    }
	std::cout << std::endl;
    HDF5Handle mid1, mid2;

    // select hyperslabs
    selectHyperslabs(mid1, mid2, shape, counter, elements, numBandsOfType);

    // write to hdf5
    H5Dwrite(dataset_id, datatype, mid2, mid1, H5P_DEFAULT, buffer.data());
    // increase counter
    counter++;
}

template <class DestIterator, class Shape, class T, int N>
void
writeHDF5Impl(DestIterator d, Shape const & shape, const hid_t file_id, const hid_t dataset_id, const hid_t datatype, ArrayVector<T> & buffer, int & counter, const int elements, const int numBandsOfType, MetaInt<N>)
{
		DestIterator dend = d + (typename DestIterator::difference_type)shape[N];
		for(; d < dend; ++d)
		{
			writeHDF5Impl(d.begin(), shape, file_id, dataset_id, datatype, buffer, counter, elements, numBandsOfType, MetaInt<N-1>());
		}
}

} // namespace detail

/** write a MultiArrayView to hdf5 file
 */
// scalar and unstrided multi arrays
template<unsigned int N, class T>
inline void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, T, UnstridedArrayTag> & array) // scalar
{
	writeHDF5(filePath, pathInFile, array, detail::getH5DataType<T>(), 1);
}

// non-scalar (TinyVector) and unstrided multi arrays
template<unsigned int N, class T, int SIZE>
inline void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, TinyVector<T, SIZE>, UnstridedArrayTag> & array)
{
	writeHDF5(filePath, pathInFile, array, detail::getH5DataType<T>(), SIZE);
}

// non-scalar (RGBValue) and unstrided multi arrays
template<unsigned int N, class T>
inline void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, RGBValue<T>, UnstridedArrayTag> & array)
{
	writeHDF5(filePath, pathInFile, array, detail::getH5DataType<T>(), 3);
}

// unstrided multi arrays
template<unsigned int N, class T>
void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, T, UnstridedArrayTag> & array, const hid_t datatype, const int numBandsOfType)
{
	HDF5Handle file_handle;
	HDF5Handle dataset_handle;
	createDataset(filePath, pathInFile, array, datatype, numBandsOfType, file_handle, dataset_handle);
	
    // Write the data to the HDF5 dataset as is
	H5Dwrite( dataset_handle, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data()); // .data() possible since void pointer!

	H5Fflush(file_handle, H5F_SCOPE_GLOBAL);
}


// scalar and strided multi arrays
template<unsigned int N, class T>
inline void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, T, StridedArrayTag> & array) // scalar
{
	writeHDF5(filePath, pathInFile, array, detail::getH5DataType<T>(), 1);
}

// non-scalar (TinyVector) and strided multi arrays
template<unsigned int N, class T, int SIZE>
inline void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, TinyVector<T, SIZE>, StridedArrayTag> & array) 
{
	writeHDF5(filePath, pathInFile, array, detail::getH5DataType<T>(), SIZE);
}

// non-scalar (RGBValue) and strided multi arrays
template<unsigned int N, class T>
inline void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, RGBValue<T>, StridedArrayTag> & array) 
{
	writeHDF5(filePath, pathInFile, array, detail::getH5DataType<T>(), 3);
}

// strided multi arrays
template<unsigned int N, class T>
void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, T, StridedArrayTag> & array, const hid_t datatype, const int numBandsOfType)
{
	HDF5Handle file_handle;
	HDF5Handle dataset_handle;
	createDataset(filePath, pathInFile, array, datatype, numBandsOfType, file_handle, dataset_handle);
	
    vigra::TinyVector<int,N> shape;
    vigra::TinyVector<int,N> stride;
    int elements = numBandsOfType;
    for(unsigned int k=0; k<N; ++k)
    {
        shape[k] = array.shape(k);
		stride[k] = array.stride(k);
        elements *= (int)shape[k];
    }
    int counter = 0;

    ArrayVector<T> buffer((int)array.shape(0));
	detail::writeHDF5Impl(array.traverser_begin(), shape, file_handle, dataset_handle, datatype, buffer, counter, elements, numBandsOfType, vigra::MetaInt<N-1>());

	H5Fflush(file_handle, H5F_SCOPE_GLOBAL);

}

namespace detail
{
struct MaxSizeFnc
{
    size_t size;

    MaxSizeFnc()
    : size(0)
    {}

    void operator()(std::string const & in)
    {
        size = in.size() > size ? 
                    in.size() :
                    size;
    }
};
}


#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 8)
/** write a numeric MultiArray as a attribute of location identifier loc
 * with name name
 */
template<size_t N, class T, class C>
void writeHDF5Attr(hid_t loc, 
                   const char* name, 
                   MultiArrayView<N, T, C> const & array)
{
    if(H5Aexists(loc, name) > 0)
        H5Adelete(loc, name);
    
    ArrayVector<hsize_t> shape(array.shape().begin(), 
                               array.shape().end());
    HDF5Handle 
        dataspace_handle(H5Screate_simple(N, shape.data(), NULL),
                         &H5Sclose, 
                         "writeToHDF5File(): unable to create dataspace.");
    
    HDF5Handle attr(H5Acreate(loc, 
                              name, 
                              detail::getH5DataType<T>(), 
                              dataspace_handle,
                              H5P_DEFAULT ,H5P_DEFAULT ),
                    &H5Aclose,
                    "writeHDF5Attr: unable to create Attribute");

    //copy data - since attributes are small - who cares!
    ArrayVector<T> buffer;
    for(int ii = 0; ii < array.size(); ++ii)
        buffer.push_back(array[ii]);
    H5Awrite(attr, detail::getH5DataType<T>(), buffer.data());
}



/** write a String MultiArray as a attribute of location identifier
 *  loc with name name
 */
template<size_t N, class C>
void writeHDF5Attr(hid_t loc, 
                   const char* name, 
                   MultiArrayView<N, std::string, C> const & array)
{
    if(H5Aexists(loc, name) > 0)
        H5Adelete(loc, name);
    
    ArrayVector<hsize_t> shape(array.shape().begin(), 
                               array.shape().end());
    HDF5Handle 
        dataspace_handle(H5Screate_simple(N, shape.data(), NULL),
                         &H5Sclose, 
                         "writeToHDF5File(): unable to create dataspace.");
    
    HDF5Handle atype(H5Tcopy (H5T_C_S1), 
                     &H5Tclose, 
                     "writeToHDF5File(): unable to create type.");

    detail::MaxSizeFnc max_size;
    max_size = std::for_each(array.data(),array.data()+ array.size(), max_size);
    H5Tset_size (atype, max_size.size);
    
    HDF5Handle attr(H5Acreate(loc, 
                              name, 
                              atype, 
                              dataspace_handle,
                              H5P_DEFAULT ,H5P_DEFAULT ),
                    &H5Aclose,
                    "writeHDF5Attr: unable to create Attribute");
    
    std::string buf ="";
    for(int ii = 0; ii < array.size(); ++ii)
    {
        buf = buf + array[ii]
                  + std::string(max_size.size - array[ii].size(), ' ');
    }
    H5Awrite(attr, atype, buf.c_str());
}

/** write an ArrayVectorView as an attribute with name to a location identifier
 */
template<class T>
inline void writeHDF5Attr(  hid_t loc,
                            const char* name,
                            ArrayVectorView<T>  & array)
{
    writeHDF5Attr(loc, name, 
                  MultiArrayView<1, T>(MultiArrayShape<1>::type(array.size()),
                                       array.data()));
}

/** write an Attribute given a file and a path in the file.
 *  the path in the file should have the format 
 *  [attribute] or /[subgroups/]dataset.attribute or
 *  /[subgroups/]group.attribute.
 *  The attribute is written to the root group, a dataset or a subgroup
 *  respectively
 */
template<class Arr>
inline void writeHDF5Attr(  std::string filePath,
                            std::string pathInFile,
                            Arr  & ar)
{
    std::string path_name(pathInFile), group_name, data_set_name, message, attr_name;
    std::string::size_type delimiter = path_name.rfind('/');
    
    //create or open file
    HDF5Handle file_id(createFile(filePath), &H5Fclose, 
                       "writeToHDF5File(): unable to open output file.");

    // get the groupname and the filename
    if(delimiter == std::string::npos)
    {
        group_name    = "/";
        data_set_name = path_name;
    }

    else
    {
        group_name = std::string(path_name.begin(), path_name.begin()+delimiter);
        data_set_name = std::string(path_name.begin()+delimiter+1, path_name.end());
    }
    delimiter = data_set_name.rfind('.');
    if(delimiter == std::string::npos)
    {
        attr_name = path_name;
        data_set_name = "/";
    }
    else
    {
        attr_name = std::string(data_set_name.begin()+delimiter+1, data_set_name.end());
        data_set_name = std::string(data_set_name.begin(), data_set_name.begin()+delimiter);
    }
    
    HDF5Handle group(openGroup(file_id, group_name), &H5Gclose, 
                     "writeToHDF5File(): Unable to create and open group. attr ver");

    if(data_set_name != "/")
    {
        HDF5Handle dset(H5Dopen(group, data_set_name.c_str(), H5P_DEFAULT), &H5Dclose,
                        "writeHDF5Attr():unable to open dataset");
        writeHDF5Attr(hid_t(dset), attr_name.c_str(), ar);
    }
    else
    {
        writeHDF5Attr(hid_t(group), attr_name.c_str(), ar);
    }

}
#endif
} // namespace vigra

#endif // VIGRA_HDF5IMPEX_HXX
