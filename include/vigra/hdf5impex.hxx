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
#include <H5Cpp.h>
#ifndef H5std_string
#define H5std_string std::string
#endif

#include "impex.hxx"
#include "multi_array.hxx"
#include "multi_impex.hxx"

namespace vigra {

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
    VIGRA_EXPORT ~HDF5ImportInfo();

    VIGRA_EXPORT const std::string& getFileName() const;

    VIGRA_EXPORT const std::string& getDatasetName() const;

	VIGRA_EXPORT const H5::H5File& getH5FileHandle() const;

	VIGRA_EXPORT const H5::DataSet& getDatasetHandle() const;

	VIGRA_EXPORT MultiArrayIndex numDimensions() const;

	VIGRA_EXPORT MultiArrayIndex shapeOfDimension(const int dim) const;

	/** Get the file type of the image associated with this
            info object.

        /** Query the pixel type of the sataset.

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
    H5::H5File m_file;
	H5::DataSet m_dataset;
    std::string m_filename, m_datasetname, m_pixeltype;
    MultiArrayIndex m_dimensions;
	ArrayVector<hsize_t> m_dims;
};

template<class type>
inline H5::DataType GetH5DataType()
{return 0;}

#define VIGRA_H5_DATATYPE(type, h5type) \
template<> \
inline H5::DataType GetH5DataType<type>() \
{return H5::PredType::h5type;}

VIGRA_H5_DATATYPE(Int8, NATIVE_INT8)
VIGRA_H5_DATATYPE(Int16, NATIVE_INT16)
VIGRA_H5_DATATYPE(Int32, NATIVE_INT32)
VIGRA_H5_DATATYPE(Int64, NATIVE_INT64)
VIGRA_H5_DATATYPE(UInt8, NATIVE_UINT8)
VIGRA_H5_DATATYPE(UInt16, NATIVE_UINT16)
VIGRA_H5_DATATYPE(UInt32, NATIVE_UINT32)
VIGRA_H5_DATATYPE(UInt64, NATIVE_UINT64)
VIGRA_H5_DATATYPE(float, NATIVE_FLOAT)
VIGRA_H5_DATATYPE(double, NATIVE_DOUBLE)

#undef VIGRA_H5_DATATYPE

template<unsigned int N, class T>
bool loadFromHDF5File(H5::H5File &file, const char* data_set_name, MultiArray<N, T> & array)
{
   H5std_string comment("");
   return loadFromHDF5File(file, data_set_name, array, comment );
}


template<unsigned int N, class T>
bool loadFromHDF5File(H5::H5File &file, const char* data_set_name, MultiArray<N, T> & array, H5std_string& comment) 
{
	try {
        H5::Exception::dontPrint();

		H5::DataSet dset = file.openDataSet(data_set_name);

		//Check type
		if(dset.getTypeClass()!=GetH5DataType<T>().getClass())
		{
			vigra_precondition( false, "Data type in file does not match expected type." );
			return false;
		}
		//Check dimensions
		if(dset.getSpace().getSimpleExtentNdims()!=N)
		{
			vigra_precondition( false, "Number of dimensions in file does not match expected number." );
			return false;
		}
		//Get dimensions size
		hsize_t size[N];
		dset.getSpace().getSimpleExtentDims(size, NULL);

	    typename MultiArrayShape<N>::type shape;
	    std::copy(size, size+N, shape.begin());
		array.reshape(shape);

		//Get the data
		dset.read(array.data(), GetH5DataType<T>());
		comment = file.getComment(data_set_name);
		return true;
    }
    catch( GroupIException not_found_error )
    {
		vigra_precondition( false, "Dataset not found in HDF5 file." );
		return false;
    }
}


template<unsigned int N, class T, class Tag>
bool loadFromHDF5File(const HDF5ImportInfo &info, MultiArrayView<N, T, Tag> & array) 
{
   H5std_string comment("");
   return loadFromHDF5File(info, array, comment);
}

template<unsigned int N, class T, class Tag>
bool loadFromHDF5File(const HDF5ImportInfo &info, MultiArrayView<N, T, Tag> & array, H5std_string& comment) 
{
	H5::DataSet dset = info.getDatasetHandle();
	H5::H5File file = info.getH5FileHandle();

	//Get dimensions size
    typedef typename MultiArrayShape<N>::type shape_type;
	shape_type shape;
	for(int i=0; i<info.numDimensions(); ++i)
		shape[i] = info.shapeOfDimension(i);

	vigra_precondition(shape == array.shape(), "loadFromHDF5File(): array must be shaped according to HDF5ImportInfo.");

	//Get the data
	//FIXME test if contiguous or use unstrided array tag
	dset.read(array.data(), GetH5DataType<T>());
	comment = file.getComment(info.getDatasetName());
	return true;
}

inline bool createAllGroups(H5::H5File &file, const char* data_set_name)
{
	H5::Exception::dontPrint();

	// iterate over all subdirectories/groups:
	std::string name(data_set_name);
	std::string::size_type begin = 0, end = name.find('/');
	while (end != std::string::npos)
	{
		std::string group(name.begin()+begin, name.begin()+end);
		//std::cout << "extracted: [" << group << "]" << std::endl;
		try {
			file.createGroup( group );
		} catch (H5::Exception &) {
			// do nothing, happens if group already exists
			//FIXME what is the right exception to catch here?
			//std::cout << e.getDetailMsg() << std::endl;
		}
		begin = end + 1;
		end = name.find('/', begin);
	}
	return true;
}

template<unsigned int N, class T, class Tag>
bool writeToHDF5File(H5::H5File &file, const char* data_set_name, const MultiArrayView<N, T, Tag> & array, const char* comment = "")
{
	H5::Exception::dontPrint();

	/*
	* Define the size of the array and create the data space for fixed
	* size dataset.
	*/
	hsize_t shape[N];              // dataset dimensions
	for(int i=0;i<N;++i)
		shape[i] = array.shape(i);
	H5::DataSpace dataspace( N, shape );

	/*
	* Define datatype for the data in the file.
	*/
	H5::DataType datatype( GetH5DataType<T>() );
	//datatype.setOrder( H5T_ORDER_LE );

	/*
	* Make sure that all the groups in the file path exist
	*/
	createAllGroups(file, data_set_name);

	/*
	* Create a new dataset within the file using defined dataspace and
	* datatype and default dataset creation properties.
	*/
	// delete the dataset if it already exists
	// in HDF5 the data is not really created but unlinked (similar to file systems)
    try {  // attempt to unlink the dataset
       file.unlink( data_set_name );
    }
    catch( H5::FileIException & unlink_error )
    {
		// do nothing, apperently there was nothing to unlink!
    }
	H5::DataSet dataset = file.createDataSet( data_set_name, datatype, dataspace );

	/*
	* Write the data to the dataset using default memory space, file
	* space, and transfer properties.
	*/
	// FIXME: contiguous? or unstrided array tag only!
	dataset.write( array.data(), GetH5DataType<T>() );

	/*
	* Write the comment
	*/
	file.setComment( data_set_name, comment);

	return true;
}

} // namespace vigra

#endif // VIGRA_HDF5IMPEX_HXX
