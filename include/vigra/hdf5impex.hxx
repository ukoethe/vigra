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
    enum PixelType { UINT8, INT16, UINT16, INT32, UINT32, FLOAT, DOUBLE };

	    /** Construct HDF5ImageImportInfo object.

            The dataset in the given HDF5 file is accessed and the properties 
			are set accordingly.
         **/
	VIGRA_EXPORT HDF5ImportInfo( const std::string &filename, const std::string &datasetname );

    VIGRA_EXPORT const std::string& getFileName() const;

    VIGRA_EXPORT const std::string& getDatasetName() const;

	VIGRA_EXPORT hid_t getH5FileHandle() const;

	VIGRA_EXPORT hid_t getDatasetHandle() const;

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
    HDF5Handle m_file, m_dataset;
    std::string m_filename, m_datasetname, m_pixeltype;
    int m_dimensions;
	ArrayVector<hsize_t> m_dims;
};

namespace detail {

template<class type>
inline hid_t getH5DataType()
{return 0;}

#define VIGRA_H5_DATATYPE(type, h5type) \
template<> \
inline hid_t getH5DataType<type>() \
{return h5type;}

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

template<unsigned int N, class T, class Tag>
void loadFromHDF5File(const HDF5ImportInfo &info, MultiArrayView<N, T, Tag> array) 
{
	vigra_precondition(N == info.numDimensions(), 
	     "loadFromHDF5File(): Array dimension disagrees with HDF5ImportInfo.numDimensions().");

    typename MultiArrayShape<N>::type shape;
    std::copy(info.shape().begin(), info.shape().end(), shape.begin());

	vigra_precondition(shape == array.shape(), 
	     "loadFromHDF5File(): Array shape disagrees with HDF5ImportInfo.");

	//Get the data
	//FIXME test if contiguous or use unstrided array tag
	vigra_postcondition(H5Dread(info.getDatasetHandle(), detail::getH5DataType<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data()) >= 0,
	       "loadFromHDF5File(): Unable to transfer data.");
}

template<unsigned int N, class T, class Alloc>
void loadFromHDF5File(HDF5ImportInfo const & info, MultiArray<N, T, Alloc> & array) 
{
    typedef typename MultiArray<N, T, Alloc>::view_type View;

    typename MultiArrayShape<N>::type shape;
    for(int k=0; k<N; ++k)
        shape[k] = info.shapeOfDimension(k);
        
	array.reshape(shape);

	return loadFromHDF5File(info, static_cast<View>(array));
}

inline hid_t createAllGroups(hid_t parent, std::string group_name)
{
	std::string::size_type begin = 0, end = group_name.find('/');
	while (end != std::string::npos)
	{
		std::string group(group_name.begin()+begin, group_name.begin()+end);
		// FIXME: handle relative and absolute paths correctly
		// FIXME: handle already existing groups
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
        parent = H5Gcreate(parent, group.c_str(), H5P_DEFAULT);
#else
        parent = H5Gcreate(parent, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
        if(parent < 0)
            return parent;
		begin = end + 1;
		end = group_name.find('/', begin);
	}
	std::string group(group_name.begin()+begin, group_name.end());
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
    return H5Gcreate(parent, group.c_str(), H5P_DEFAULT);
#else
    return H5Gcreate(parent, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
}



template<unsigned int N, class T, class Tag>
void writeToHDF5File(hid_t file, const char* path, const MultiArrayView<N, T, Tag> & array, const char* comment = "")
{
	std::string path_name(path), group_name, data_set_name, message;
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
	    group = HDF5Handle(createAllGroups(file, group_name), &H5Gclose, "writeToHDF5File(): Unable to open group.");
	}
	else
	{
	    group = HDF5Handle(file, 0, "");
	}
	
	hsize_t shape[N];              // dataset dimensions
	std::copy(array.shape().begin(), array.shape().end(), shape);

#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6)
    if(H5LTfind_dataset(group, data_set_name.c_str()))
    {
        if(H5Ldelete( group, data_set_name.c_str(), H5P_DEFAULT ) < 0)
        {
            vigra_postcondition(false, "writeToHDF5File(): Unable to delete existing data.");
        }
    }
#endif
    HDF5Handle dataset(H5LTmake_dataset(group, data_set_name.c_str(), N, shape, detail::getH5DataType<T>(), array.data()),
                       &H5Dclose, "writeToHDF5File(): Unable to write data set.");
}

} // namespace vigra

#endif // VIGRA_HDF5IMPEX_HXX
