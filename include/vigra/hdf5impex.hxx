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
#include <numeric>

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

/** Handle to automatically to a hid_t object. takes care of deallocation
 * once the hid_t is not needed anymore. 
 *
 * contains an implicit cast to hid_t - this means it can be used exactly like
 * a hid_t otherwise.
 * */

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
    
    /** create a HDF5 handle. 
     * \param hid_t returned by a H5X...() method (e.g. H5Acreate() or 
     *          H5Dcreate etc.
     * \param destructor function pointer to the right deallocation method
     *          (e.g. H5Ddelete, H5Adelete,  H5Tdelete etc...
     * \param error_message runtime exception string to be thrown if the 
     *        allocation or deallocation of memory pointed to by hid_t 
     */
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
    
    /** Implicit cast to hid_t
     */
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
VIGRA_H5_DATATYPE(std::string, H5T_C_S1)

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
void loadFromHDF5File(const HDF5ImportInfo &info, 
                      MultiArrayView<N, T, Tag> array, 
                      const bool rowMajorOrder = false) 
{


    typename MultiArrayShape<N>::type shape;
    for(unsigned int k=0; k<N; ++k)
        shape[k] = (MultiArrayIndex)info.shapeOfDimension(k);
    
    int elements = std::accumulate(info.shape().begin(), 
                                   info.shape().end(),
                                   1,
                                   std::multiplies<int>());
    
    vigra_precondition((array.size() == elements),
        "loadFromHDF5File(): Array size disagrees with number of elements in"
        "HDF5 Dataset");
    
    // commented out as we only have to ensure that the number of elements 
    // are the same. Remove comment block if this is ok.
    // vigra_precondition((N == info.numDimensions()),
    //    "loadFromHDF5File(): Array dimension disagrees with HDF5ImportInfo.numDimensions().");
    // vigra_precondition(shape == array.shape(), 
    //     "loadFromHDF5File(): Array shape disagrees with HDF5ImportInfo.");

    hid_t dataset = info.getDatasetHandle();
    hid_t cparms = H5Dget_create_plist (dataset);
    if (H5D_CHUNKED == H5Pget_layout (cparms))
    {    
        hsize_t      dimsr[N];
        for(int ii = 0; ii < N;++ii)
            dimsr[ii] = array.shape(ii);
        HDF5Handle filespace(H5Dget_space (dataset), H5Sclose, "space error");
        HDF5Handle memspace(H5Screate_simple (N,dimsr,NULL), H5Sclose, "second space error");

        if(rowMajorOrder)
        {
            H5Dread (dataset, detail::getH5DataType<T>(), memspace, filespace,
                              H5P_DEFAULT, array.data());
        }
        else
        {
            ArrayVector<T> buffer(elements);
            H5Dread (dataset, detail::getH5DataType<T>(), memspace, filespace,
                              H5P_DEFAULT, buffer.data());
            vigra::TinyVector<int,N> strideNew;
            vigra::TinyVector<int,N> shapeNew;
            for(unsigned int k=0; k<N; ++k)
            {
                strideNew[k] = array.stride(N-1-k);
                shapeNew[k] = array.shape(N-1-k);
            }
            MultiArrayView<N, T, StridedArrayTag> arrayNew (shapeNew, 
                                                            strideNew, 
                                                            array.data());
            for(int ii = 0; ii < elements; ++ii)
                arrayNew[ii] = buffer[ii];
        }
       return;
    }
    H5Pclose (cparms);


    //Get the data
    int counter = 0;
    if(rowMajorOrder)
    {
        ArrayVector<T> buffer(shape[0]);
        detail::readHDF5Impl(array.traverser_begin(), 
                             shape, 
                             info.getDatasetHandle(), 
                             buffer, 
                             counter, 
                             elements, 
                             vigra::MetaInt<N-1>());
    } else {
        vigra::TinyVector<int,N> strideNew;
        vigra::TinyVector<int,N> shapeNew;
        for(unsigned int k=0; k<N; ++k)
        {
            strideNew[k] = array.stride(N-1-k);
            shapeNew[k] = array.shape(N-1-k);
        }
        MultiArrayView<N, T, StridedArrayTag> arrayNew (shapeNew, 
                                                        strideNew, 
                                                        array.data());
        ArrayVector<T> buffer(arrayNew.shape(0));
        detail::readHDF5Impl(arrayNew.traverser_begin(), 
                             arrayNew.shape(), 
                             info.getDatasetHandle(), 
                             buffer, 
                             counter, 
                             elements, 
                             vigra::MetaInt<N-1>());
    }

}

/** Overloaded load method for loading string multiArrays
 */
template<unsigned int N, class Tag>
void loadFromHDF5File(const HDF5ImportInfo &info, 
                      MultiArrayView<N, std::string, Tag> array, 
                      const bool rowMajorOrder = false) 
{
    typename MultiArrayShape<N>::type shape;
    for(unsigned int k=0; k<N; ++k)
        shape[k] = (MultiArrayIndex)info.shapeOfDimension(k);

    int elements = std::accumulate(info.shape().begin(), 
                                   info.shape().end(),
                                   1,
                                   std::multiplies<int>());
    std::cerr << array.size() << " " << elements << std::endl;
    vigra_precondition((array.size() == elements),
        "loadFromHDF5File(): Array size disagrees with number of elements in"
        "HDF5 Dataset");
    
    // commented out as we only have to ensure that the number of elements 
    // are the same. Remove comment block if this is ok.
    // vigra_precondition((N == info.numDimensions()),
    //    "loadFromHDF5File(): Array dimension disagrees with HDF5ImportInfo.numDimensions().");
    // vigra_precondition(shape == array.shape(), 
    //     "loadFromHDF5File(): Array shape disagrees with HDF5ImportInfo.");

    HDF5Handle 
        dataspace_handle(H5Dget_space(info.getDatasetHandle()),
                     &H5Sclose, 
                     "loadToHDF5File(): unable to create dataspace.");

    HDF5Handle 
        datatype_handle(H5Dget_type(info.getDatasetHandle()),
                     &H5Tclose, 
                     "loadToHDF5File(): unable to create dataspace.");

    //TODO add precondition to check whether this is really a string type
    // not checking it will not make the method crash but will produce wierd
    // results when you try to open a numeric array.
    hssize_t type_size = H5Tget_size(datatype_handle);
    hssize_t size = H5Sget_simple_extent_npoints(dataspace_handle);

    char buffer[type_size * size];
    H5Dread(info.getDatasetHandle(), 
            datatype_handle, 
            H5S_ALL, 
            H5S_ALL, 
            H5P_DEFAULT, 
            buffer); 
    char* iter = buffer;
    if(rowMajorOrder)
    {
        for(int ii = 0; ii < array.size();++ii)
        {
            array[ii] = std::string(iter, type_size);
            while(array[ii].find('\0') != std::string::npos)
                array[ii].erase(array[ii].find('\0'), 1);
            iter += type_size;
        }
    } else {
       
        for(int ii = 0; ii < array.size();++ii)
        {
            array.transpose()[ii] = std::string(iter, type_size);
            while(array.transpose()[ii].find('\0') != std::string::npos)
                array.transpose()[ii].erase(array.transpose()[ii].find('\0'), 1);
            iter += type_size;
        }
    }

}

/** recursively open a group  given by a path string
 */
inline hid_t openGroup(hid_t parent, std::string group_name)
{
    //std::cout << group_name << std::endl;
    if(group_name == "/")
        return H5Gopen(parent, group_name.c_str(), H5P_DEFAULT);
    if(group_name[0] == '/')
        group_name.erase(0,1);
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

/** recursively create a path of groups given a path string
 */
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

/** delete a Dataset
 */
inline void deleteDataset(hid_t parent, std::string dataset_name)
{
    // delete existing data and create new dataset
    if(H5LTfind_dataset(parent, dataset_name.c_str()))
    {
        //std::cout << "dataset already exists" << std::endl;
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
		if(H5Gunlink(parent, dataset_name.c_str()) < 0)
        {
            vigra_postcondition(false, "deleteDataset(): Unable to delete existing data.");
        }
#else
		if(H5Ldelete(parent, dataset_name.c_str(), H5P_DEFAULT ) < 0)
        {
            vigra_postcondition(false, "deleteDataset(): Unable to delete existing data.");
        }
#endif
    } 
}

/** create a HDF5 file if it doesn't exist. 
 *  default behavior: append to existing file if it exists.
 *  if append_ = false: overwrite existing file
 */
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

/** write a MultiArrayView to hdf5 file
 */
template<unsigned int N, class T, class Tag>
void writeToHDF5File(const char* filePath, const char* pathInFile, const MultiArrayView<N, T, Tag> & array, const bool rowMajorOrder = false)
{

    std::string path_name(pathInFile), group_name, data_set_name, message;
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

    // create all groups
    HDF5Handle group(createGroup(file_id, group_name), &H5Gclose, 
                     "writeToHDF5File(): Unable to create and open group. generic v");

    // delete the dataset if it already exists
    deleteDataset(group, data_set_name);

    // create dataspace
    hsize_t shape[N];
    std::copy(array.shape().begin(), array.shape().end(), shape);
    HDF5Handle dataspace_handle(H5Screate_simple(N, shape, NULL),
                                &H5Sclose, "writeToHDF5File(): unable to create dataspace.");

    //alloc memory for dataset. 
    HDF5Handle dataset_handle(H5Dcreate(group, 
                                        data_set_name.c_str(), 
                                        detail::getH5DataType<T>(), 
                                        dataspace_handle, 
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                              &H5Dclose, "writeToHDF5File(): unable to create dataset.");
    
    // Write the data to the HDF5 dataset
    int elements = array.size();
    int counter = 0;
    
    if(rowMajorOrder)
    {
        ArrayVector<T> buffer(array.shape(0));
        detail::writeHDF5Impl(array.traverser_begin(), 
                              shape, 
                              file_id, 
                              dataset_handle, 
                              buffer, counter, 
                              elements, 
                              vigra::MetaInt<N-1>());
    } else {
        // for column major order we have to reverse the shape and strides before calling the write function
        vigra::TinyVector<int,N> strideNew, shapeNew;
        for(unsigned int k=0; k<N; ++k)
        {
            strideNew[k] = array.stride(N-1-k);
            shapeNew[k] = array.shape(N-1-k);
        }
        MultiArrayView<N, T, StridedArrayTag> arrayNew (shapeNew, strideNew, array.data());
        ArrayVector<T> buffer((int)arrayNew.shape(0));
        detail::writeHDF5Impl(arrayNew.traverser_begin(), 
                              arrayNew.shape(), 
                              file_id, 
                              dataset_handle, 
                              buffer, 
                              counter, 
                              elements, 
                              vigra::MetaInt<N-1>());
    }
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
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


/** write a string MultiArray array to pathInFile in the hdf5 file filePath
 */
template<unsigned int N, class Tag>
void writeToHDF5File(const char* filePath, 
                     const char* pathInFile, 
                     const MultiArrayView<N, std::string, Tag> & array, 
                     const bool rowMajorOrder = false)
{
    std::string path_name(pathInFile), group_name, data_set_name, message;
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
    // create all groups
    HDF5Handle group(createGroup(file_id, group_name), &H5Gclose, 
                     "writeToHDF5File(): Unable to create and open group. str ver");

    // delete the dataset if it already exists
    deleteDataset(group,data_set_name);

    // create dataspace
    hsize_t shape[N];
    std::copy(array.shape().begin(), array.shape().end(), shape);
    HDF5Handle dataspace_handle(H5Screate_simple(N, shape, NULL),
                                &H5Sclose, "writeToHDF5File(): unable to create dataspace.");
    
    HDF5Handle atype(H5Tcopy (H5T_C_S1), &H5Tclose, 
                    "writeToHDF5File(): unable to create type.");
    detail::MaxSizeFnc max_size;
    inspectMultiArray(srcMultiArrayRange(array), max_size);
    H5Tset_size (atype, max_size.size);

    //alloc memory for dataset. 
    HDF5Handle dataset_handle(H5Dcreate(group, 
                                        data_set_name.c_str(), 
                                        atype, 
                                        dataspace_handle, 
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                              &H5Dclose, "writeToHDF5File(): unable to create dataset.");
    std::string buf ="";
    for(int ii = 0; ii < array.size(); ++ii)
    {
        buf = buf + array[ii] + std::string(max_size.size - array[ii].size(), '\0');
    }
    H5Dwrite (dataset_handle, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.c_str());
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
}

// The creation and reading of HDF5 attributes is currently only supported 
// with HDF5 1.8 onwards. This is because I have no clue how to open an
// attribute by name with the earlier version. 
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 8)

/** Quick and dirty hack for loading arrays valued attributes without the 
 *  tedious usage of an info object. Down side is that you have to know the
 *  name and type of the attribute
 */ 
template<class T>
bool loadHDF5Attr(hid_t loc, const char* name, ArrayVector<T> & array)
{
        if(H5Aexists(loc, name) == 0)
            return false;

        HDF5Handle attr(H5Aopen_by_name(loc,".", name, H5P_DEFAULT, H5P_DEFAULT), 
                        &H5Aclose,
                        "loadHDF5Attr: unable to create Attribute");
        HDF5Handle 
            dataspace_handle(H5Aget_space(attr),
                         &H5Sclose, 
                         "writeToHDF5File(): unable to create dataspace.");
        HDF5Handle 
            datatype_handle(H5Aget_type(attr),
                         &H5Tclose, 
                         "writeToHDF5File(): unable to create dataspace.");
        if((H5Tget_class(datatype_handle) == H5T_STRING))
            return false;
        hssize_t size = H5Sget_simple_extent_npoints(dataspace_handle);
        
        array.resize(size);
        H5Aread(attr, detail::getH5DataType<T>(), array.data());
        return true;
}
/** loads a string valued array of attributes
 */
bool loadHDF5Attr(hid_t loc, const char* name, ArrayVector<std::string> & array)
{
        if(H5Aexists(loc, name) == 0)
            return false;

        HDF5Handle attr(H5Aopen_by_name(loc,".", name, H5P_DEFAULT, H5P_DEFAULT), 
                        &H5Aclose,
                        "loadHDF5Attr: unable to create Attribute");
        HDF5Handle 
            dataspace_handle(H5Aget_space(attr),
                         &H5Sclose, 
                         "writeToHDF5File(): unable to create dataspace.");

        HDF5Handle 
            datatype_handle(H5Aget_type(attr),
                         &H5Tclose, 
                         "writeToHDF5File(): unable to create dataspace.");
        if(!(H5Tget_class(datatype_handle) == H5T_STRING))
            return false;
        hssize_t type_size = H5Tget_size(datatype_handle);
        hssize_t size = H5Sget_simple_extent_npoints(dataspace_handle);

        char buffer[type_size * size];
        H5Aread(attr, datatype_handle, buffer);
        char* iter = buffer;
        for(int ii = 0; ii < size; ++ii)
        {
            array.push_back(std::string(iter, type_size));
            while(array.back().find('\0') != std::string::npos)
                array.back().erase(array[ii].find('\0'), 1);
            iter+= type_size; 
        }
        return true;
}

/** load given file and path in file. gorups are delimited by / and the 
 * attribute by .
 * eg /group/subgroup/dataset.attribute. Why for lords sake couldn't the guys
 * at HDF5 implemented it this way as well?
 */
template<class T>
bool loadHDF5Attr(const char* filePath, const char* pathInFile, ArrayVector<T> & array)
{
    std::string path_name(pathInFile), group_name, data_set_name, message, attr_name;
    std::string::size_type delimiter = path_name.rfind('/');
    
    //create or open file
    HDF5Handle file_id(createFile(filePath), &H5Fclose, 
                       "loadHDF5Attr(): unable to open output file.");

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
                     "loadHDF5Attr(): Unable to create and open group. attr ver");

    if(data_set_name != "/")
    {
        HDF5Handle dset(H5Dopen(group, data_set_name.c_str(), H5P_DEFAULT), &H5Dclose,
                        "loadHDF5Attr():unable to open dataset");
        return loadHDF5Attr(hid_t(dset), attr_name.c_str(), array);
    }
    else
    {
        return loadHDF5Attr(hid_t(group), attr_name.c_str(), array);
    }

}


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
inline void writeHDF5Attr(  const char* filePath,
                            const char* pathInFile,
                            Arr  & ar)
{
    std::string path_name(pathInFile), group_name, data_set_name, message, attr_name;
    std::string::size_type delimiter = path_name.rfind('/');
    
    //create or open file
    HDF5Handle file_id(createFile(filePath), &H5Fclose, 
                       "writeHDF5Attr(): unable to open output file.");

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
                     "writeHDF5Attr(): Unable to create and open group. attr ver");

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
