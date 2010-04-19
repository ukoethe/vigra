/************************************************************************/
/*                                                                      */
/*             Copyright 2009-2010 by Ullrich Koethe                    */
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

#ifdef HasHDF5

#include "vigra/hdf5impex.hxx"
#include "vigra/multi_array.hxx"
#include <iostream>
#include <cstring>
#include <cstdio>

namespace vigra {

HDF5ImportInfo::HDF5ImportInfo(const char* filePath, const char* pathInFile)
{
    m_file_handle = HDF5Handle(H5Fopen(filePath, H5F_ACC_RDONLY, H5P_DEFAULT),
                               &H5Fclose, "HDF5ImportInfo(): Unable to open file.");


    m_dataset_handle = HDF5Handle(H5Dopen(m_file_handle, pathInFile, H5P_DEFAULT),
                                  &H5Dclose, "HDF5ImportInfo(): Unable to open dataset.");


    //DataSet dset = m_file.openDataSet(datasetname);
    m_filename = filePath;
    m_path = pathInFile;
    HDF5Handle dataspace_handle(H5Dget_space(m_dataset_handle),
                                &H5Sclose, "HDF5ImportInfo(): could not access dataset dataspace.");
    m_dimensions = H5Sget_simple_extent_ndims(dataspace_handle);
    //m_dimensions = dset.getSpace().getSimpleExtentNdims();
    
    //why?
    //vigra_precondition( m_dimensions>=2, "HDF5ImportInfo(): Number of dimensions is lower than 2. Not an image!" );

	hid_t datatype = H5Dget_type(m_dataset_handle);
	H5T_class_t dataclass = H5Tget_class(datatype);
	size_t datasize  = H5Tget_size(datatype);
	H5T_sign_t datasign  = H5Tget_sign(datatype);

	if(dataclass == H5T_FLOAT)
	{
		if(datasize == 4)
			m_pixeltype = "FLOAT";
		else if(datasize == 8)
			m_pixeltype = "DOUBLE";
	}
	else if(dataclass == H5T_INTEGER)	
	{
		if(datasign == H5T_SGN_NONE)
		{
			if(datasize ==  1)
				m_pixeltype = "UINT8";
			else if(datasize == 2)
				m_pixeltype = "UINT16";
			else if(datasize == 4)
				m_pixeltype = "UINT32";
			else if(datasize == 8)
				m_pixeltype = "UINT64";
		}
		else
		{
			if(datasize ==  1)
				m_pixeltype = "INT8";
			else if(datasize == 2)
				m_pixeltype = "INT16";
			else if(datasize == 4)
				m_pixeltype = "INT32";
			else if(datasize == 8)
				m_pixeltype = "INT64";
		}
	}

    ArrayVector<hsize_t>::size_type ndims = ArrayVector<hsize_t>::size_type(m_dimensions);
    m_dims.resize(ndims);
    ArrayVector<hsize_t> size(ndims);
    ArrayVector<hsize_t> maxdims(ndims);
    H5Sget_simple_extent_dims(dataspace_handle, size.data(), maxdims.data());
    //dset.getSpace().getSimpleExtentDims(size, NULL);
	// invert the dimensions to guarantee c-order
	for(ArrayVector<hsize_t>::size_type i=0; i<ndims; i++) {
        m_dims[i] = size[ndims-1-i];
		//std::cout << "m_dims[" << i << "]=" << m_dims[i] << std::endl;
	}
}

HDF5ImportInfo::~HDF5ImportInfo()
{}


HDF5ImportInfo::PixelType HDF5ImportInfo::pixelType() const
{
   const std::string pixeltype=HDF5ImportInfo::getPixelType();
   if (pixeltype == "UINT8")
       return HDF5ImportInfo::UINT8;
   if (pixeltype == "UINT16")
     return HDF5ImportInfo::UINT16;
   if (pixeltype == "UINT32")
     return HDF5ImportInfo::UINT32;
   if (pixeltype == "UINT64")
     return HDF5ImportInfo::UINT64;
   if (pixeltype == "INT8")
       return HDF5ImportInfo::INT8;
   if (pixeltype == "INT16")
     return HDF5ImportInfo::INT16;
   if (pixeltype == "INT32")
     return HDF5ImportInfo::INT32;
   if (pixeltype == "INT64")
     return HDF5ImportInfo::INT64;
   if (pixeltype == "FLOAT")
     return HDF5ImportInfo::FLOAT;
   if (pixeltype == "DOUBLE")
     return HDF5ImportInfo::DOUBLE;
   vigra_fail( "internal error: unknown pixel type" );
   return HDF5ImportInfo::PixelType();
}

const char * HDF5ImportInfo::getPixelType() const
{
    return m_pixeltype.c_str();
}

MultiArrayIndex HDF5ImportInfo::shapeOfDimension(const int dim) const 
{ 
    return MultiArrayIndex(m_dims[dim]); 
}

MultiArrayIndex HDF5ImportInfo::numDimensions() const 
{ 
    return MultiArrayIndex(m_dimensions); 
}

const std::string & HDF5ImportInfo::getPathInFile() const 
{ 
    return m_path; 
}

const std::string & HDF5ImportInfo::getFilePath() const 
{ 
    return m_filename; 
}

hid_t HDF5ImportInfo::getH5FileHandle() const 
{ 
    return m_file_handle; 
}

hid_t HDF5ImportInfo::getDatasetHandle() const 
{ 
    return m_dataset_handle; 
}

} // namespace vigra

#endif // HasHDF5
