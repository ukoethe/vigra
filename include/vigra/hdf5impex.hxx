#ifndef VIGRA_HDF5IMPEX_HXX
#define VIGRA_HDF5IMPEX_HXX

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include <H5Cpp.h>
#ifndef H5std_string
#define H5std_string std::string
#endif

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include "impex.hxx"
#include "multi_array.hxx"
#include "multi_impex.hxx"
#include <boost/xpressive/xpressive.hpp>

using namespace boost::xpressive;
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

	VIGRA_EXPORT const H5File& getH5FileHandle() const;

	VIGRA_EXPORT const DataSet& getDatasetHandle() const;

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
    H5File m_file;
	DataSet m_dataset;
    std::string m_filename, m_datasetname, m_pixeltype;
    MultiArrayIndex m_dimensions;
	ArrayVector<int> m_dims;
};

template<class type>
DataType GetH5DataType()
{return 0;}
template<>
DataType GetH5DataType<Int8>()
{return PredType::NATIVE_INT8;}
template<>
DataType GetH5DataType<Int16>()
{return PredType::NATIVE_INT16;}
template<>
DataType GetH5DataType<Int32>()
{return PredType::NATIVE_INT32;}
template<>
DataType GetH5DataType<Int64>()
{return PredType::NATIVE_INT64;}
template<>
DataType GetH5DataType<UInt8>()
{return PredType::NATIVE_UINT8;}
template<>
DataType GetH5DataType<UInt16>()
{return PredType::NATIVE_UINT16;}
template<>
DataType GetH5DataType<UInt32>()
{return PredType::NATIVE_UINT32;}
template<>
DataType GetH5DataType<UInt64>()
{return PredType::NATIVE_UINT64;}
template<>
DataType GetH5DataType<float>()
{return PredType::NATIVE_FLOAT;}
template<>
DataType GetH5DataType<double>()
{return PredType::NATIVE_DOUBLE;}


template<unsigned int N, class T>
bool loadFromHDF5File(H5File &file, const char* data_set_name, MultiArray<N, T> & array)
{
   H5std_string comment("");
   return loadFromHDF5File(file, data_set_name, array, comment );
}


template<unsigned int N, class T>
bool loadFromHDF5File(H5File &file, const char* data_set_name, MultiArray<N, T> & array, H5std_string& comment) 
{
	try {
        Exception::dontPrint();

		DataSet dset = file.openDataSet(data_set_name);

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
		typedef typename MultiArrayShape<N>::type shape_type;
		shape_type shape;
		hsize_t size[N];
		dset.getSpace().getSimpleExtentDims(size, NULL);
		for(int i=0; i<N; ++i)
			shape[i] = size[i];//N-i-1];
		//Create multi array
		array = MultiArray<N, T>(shape);
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
	DataSet dset = info.getDatasetHandle();
	H5File file = info.getH5FileHandle();

	//Get dimensions size
        typedef typename MultiArrayShape<N>::type shape_type;
	shape_type shape;
	for(int i=0; i<info.numDimensions(); ++i)
		shape[i] = info.shapeOfDimension(i);

	vigra_precondition(shape == array.shape(), "loadFromHDF5File(): array must be shaped according to HDF5ImportInfo.");

	//Get the data
	dset.read(array.data(), GetH5DataType<T>());
	comment = file.getComment(info.getDatasetName());
	return true;
}


bool createAllGroups(H5File &file, const char* data_set_name)
{
	sregex re = sregex::compile("[A-Za-z0-9_.,����\\-\\+]+[\\/]{1}"); // find groups (the last token ensures that a dataset is not considered a group)

	// iterate over all subdirectories/groups:
	std::string data_set_name_str(data_set_name);
	//std::cout << "input string is: " << data_set_name_str << std::endl;
	sregex_token_iterator begin( data_set_name_str.begin(), data_set_name_str.end(), re ), end;
	// create the group(s)
	std::string group("");
	while (begin != end)
	{
		group = group + *begin++; // subgroups
		//std::cout << "extracted: [" << group << "]" << std::endl;
		file.createGroup( group );
	}
	return true;
}

template<unsigned int N, class T, class Tag>
bool writeToHDF5File(H5File &file, const char* data_set_name, const MultiArrayView<N, T, Tag> & array, const char* comment = "")
{
	/*
	* Define the size of the array and create the data space for fixed
	* size dataset.
	*/
	hsize_t shape[N];              // dataset dimensions
	for(int i=0;i<N;++i)
		shape[i] = array.shape(i);
	DataSpace dataspace( N, shape );

	/*
	* Define datatype for the data in the file.
	*/
	DataType datatype( GetH5DataType<T>() );
	//datatype.setOrder( H5T_ORDER_LE );

	/*
	* Make sure that all the groups in the file path exist
	*/
	createAllGroups(file, data_set_name);

	/*
	* Create a new dataset within the file using defined dataspace and
	* datatype and default dataset creation properties.
	*/
	DataSet dataset = file.createDataSet( data_set_name, datatype, dataspace );

	/*
	* Write the data to the dataset using default memory space, file
	* space, and transfer properties.
	*/
	dataset.write( array.data(), GetH5DataType<T>() );

	/*
	* Write the comment
	*/
	file.setComment( data_set_name, comment);

	return true;
}

} // namespace vigra

#endif // VIGRA_HDF5IMPEX_HXX
