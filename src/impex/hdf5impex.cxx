#include "vigra/hdf5impex.hxx"
#include "vigra/numpy_array.hxx"
#include "vigra/numpy_array_converters.hxx"
#include <iostream>
#include <cstring>
#include <cstdio>

using namespace vigra;

HDF5ImportInfo::HDF5ImportInfo(const std::string &filename, const std::string &datasetname)
{
	try {
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        Exception::dontPrint();

		m_file = H5File( filename, H5F_ACC_RDONLY );

		DataSet dset = m_file.openDataSet(datasetname);

		m_filename = filename;
		m_datasetname = datasetname;
		m_dimensions = dset.getSpace().getSimpleExtentNdims();
		m_dataset = dset;

		vigra_precondition( m_dimensions>=2, "Number of dimensions is lower than 2. Not an image!" );

		if(dset.getTypeClass()==GetH5DataType<float>().getClass())
			m_pixeltype = "FLOAT";
		if(dset.getTypeClass()==GetH5DataType<UInt8>().getClass())
			m_pixeltype = "UINT8";
		if(dset.getTypeClass()==GetH5DataType<Int8>().getClass())
			m_pixeltype = "INT8";
		if(dset.getTypeClass()==GetH5DataType<UInt16>().getClass())
			m_pixeltype = "UINT16";
		if(dset.getTypeClass()==GetH5DataType<Int16>().getClass())
			m_pixeltype = "INT16";
		if(dset.getTypeClass()==GetH5DataType<UInt32>().getClass())
			m_pixeltype = "UINT32";
		if(dset.getTypeClass()==GetH5DataType<Int32>().getClass())
			m_pixeltype = "INT32";
		if(dset.getTypeClass()==GetH5DataType<double>().getClass())
			m_pixeltype = "DOUBLE";

		m_dims = ArrayVector<int>(m_dimensions);
		hsize_t* size = new hsize_t[m_dimensions];
		dset.getSpace().getSimpleExtentDims(size, NULL);
		for(int i=0; i<m_dimensions; ++i)
			m_dims[i] = size[i];
		delete size;

    }
    catch( GroupIException not_found_error )
    {
		vigra_precondition( false, "Dataset not found in HDF5 file." );
    }
}

HDF5ImportInfo::~HDF5ImportInfo()
{
}

HDF5ImportInfo::PixelType HDF5ImportInfo::pixelType() const
{
   const std::string pixeltype=HDF5ImportInfo::getPixelType();
   if (pixeltype == "UINT8")
	   return HDF5ImportInfo::UINT8;
   if (pixeltype == "INT16")
     return HDF5ImportInfo::INT16;
   if (pixeltype == "UINT16")
     return HDF5ImportInfo::UINT16;
   if (pixeltype == "INT32")
     return HDF5ImportInfo::INT32;
   if (pixeltype == "UINT32")
     return HDF5ImportInfo::UINT32;
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
MultiArrayIndex HDF5ImportInfo::shapeOfDimension(const int dim) const { return m_dims[dim]; };
MultiArrayIndex HDF5ImportInfo::numDimensions() const { return m_dimensions; }
const std::string & HDF5ImportInfo::getDatasetName() const { return m_datasetname; }
const std::string & HDF5ImportInfo::getFileName() const { return m_filename; }
const H5File& HDF5ImportInfo::getH5FileHandle() const { return m_file; }
const DataSet& HDF5ImportInfo::getDatasetHandle() const { return m_dataset; }


int main (void)
{
    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try
    {
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        //Exception::dontPrint();

        /*
         * Create a file if it does not exist and open for reading/writing
         */
		H5File out_file( "Z:\\Ilastik\\vigranumpy\\vigranumpy.new\\src\\testfile.hdf5", H5F_ACC_TRUNC );

        /*
        * Create property list for a dataset and set up fill values.
        */
        int fillvalue = 0;   /* Fill value for the dataset */
        DSetCreatPropList plist;
        plist.setFillValue(PredType::NATIVE_INT, &fillvalue);

		/*
		* Create some data and write it to the file
		*/
		MultiArray<3,double> out_data1(MultiArrayShape<3>::type(200, 100, 3));
		out_data1.init(5);
		bool b1 = writeToHDF5File(out_file, "my_data", out_data1);

		MultiArray<4,int> out_data2(MultiArrayShape<4>::type(500, 300, 5, 6));
		out_data2.init(3);
		bool b2 = writeToHDF5File(out_file, "/group/subgroup/my_data2", out_data2, "in subgroup");

		ImageImportInfo info("Z:\\Ilastik\\vigranumpy\\vigranumpy.new\\src\\giraffe.jpg");
		MultiArray<2,int> out_data3(MultiArrayShape<2>::type(info.width(), info.height()));
		importImage(info, destImage(out_data3));
		bool b3 = writeToHDF5File(out_file, "giraffe", out_data3);

		std::cout<<"Writing(old)... Success="<<b1*b2*b3<<std::endl;

		out_file.close();

		H5File out_file2( "Z:\\Ilastik\\vigranumpy\\vigranumpy.new\\src\\testfile.hdf5", H5F_ACC_RDWR );
		bool b4 = writeToHDF5File(out_file2, "giraffe3", out_data3);
		std::cout<<"Writing(new)... Success="<<b4<<std::endl;
		out_file2.close();

		/*
		* Read the data
		*/
		H5File in_file("Z:\\Ilastik\\vigranumpy\\vigranumpy.new\\src\\testfile.hdf5", H5F_ACC_RDONLY );
		H5std_string comment;

		MultiArray<3,double> in_data1;
		bool b5 = loadFromHDF5File(in_file, "my_data", in_data1, comment);
		std::cout<<"Size set 1: "<<in_data1.shape()<<std::endl;
		std::cout<<comment<<std::endl;

    	MultiArray<4,int> in_data2;
		bool b6 = loadFromHDF5File(in_file, "/group/subgroup/my_data2", in_data2, comment);
		std::cout<<"Size set 2: "<<in_data2.shape()<<std::endl;
		std::cout<<comment<<std::endl;

    	MultiArray<2,int> in_data3;
		bool b7 = loadFromHDF5File(in_file, "giraffe", in_data3, comment);
		std::cout<<"Size set 3: "<<in_data3.shape()<<std::endl;
		std::cout<<comment<<std::endl;

		std::cout<<"Reading(old)... Success="<<b5*b6*b7<<std::endl;

		in_file.close();
		
		HDF5ImportInfo infoHDF5("Z:\\Ilastik\\vigranumpy\\vigranumpy.new\\src\\testfile.hdf5", "giraffe");
        MultiArray<2,int> in_data4(MultiArrayShape<2>::type(infoHDF5.shapeOfDimension(0), infoHDF5.shapeOfDimension(1)));
        bool b8 = loadFromHDF5File(infoHDF5, in_data4);
		std::cout<<"Size set 4: "<<in_data4.shape()<<std::endl;

		std::cout<<"Reading(new1)... Success="<<b8<<std::endl;

		/*
	    HDF5ImportInfo infoHDF5b("Z:\\Ilastik\\vigranumpy\\vigranumpy.new\\src\\testfile.hdf5", "my_data");
		std::cout<<"test1"<<std::endl;
        MultiArray<2, RGBValue<double> > in_data5(MultiArrayShape<2>::type(infoHDF5b.shapeOfDimension(0), infoHDF5b.shapeOfDimension(1)));
		std::cout<<"test2"<<std::endl;
        bool b9 = loadFromHDF5File(infoHDF5b, in_data5);
		std::cout<<"Size set 5: "<<in_data5.shape()<<std::endl;

		std::cout<<"Reading(new2)... Success="<<b9<<std::endl;
		*/

		ImageExportInfo info2("Z:\\Ilastik\\vigranumpy\\vigranumpy.new\\src\\giraffentest.jpg");
		exportImage(srcImageRange(in_data3), info2);

   }  // end of try block

   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
        error.printError();
        return -1;
   }

   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
        error.printError();
        return -1;
   }

   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
        error.printError();
        return -1;
   }

   return 0;
}


