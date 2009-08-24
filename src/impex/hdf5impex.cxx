#include "vigra/hdf5impex.hxx"
#include "vigra/multi_array.hxx"
#include <iostream>
#include <cstring>
#include <cstdio>

namespace vigra {

HDF5ImportInfo::HDF5ImportInfo(const std::string &filename, const std::string &path)
: m_filename(filename),
  m_datasetname(path)
{
	std::string path_name(path), group_name, data_set_name, message;

    message = std::string("HDF5ImportInfo(): Unable to open file '") + filename + "'.";
    m_file = HDF5Handle(H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
                        &H5Fclose, message.c_str());
    
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
        message = std::string("HDF5ImportInfo(): Unable to open group '") + group_name + "'.";
#if H5_VERS_MINOR <= 6
	    group = HDF5Handle(H5Gopen(m_file, group_name.c_str()),
	                       &H5Gclose, message.c_str());
#else
	    group = HDF5Handle(H5Gopen(m_file, group_name.c_str(), H5P_DEFAULT),
	                       &H5Gclose, message.c_str());
#endif
	}
	else
	{
	    group = HDF5Handle(m_file, 0, "");
	}
	
	message = std::string("HDF5ImportInfo(): Unable to open data set '") + path + "'.";
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
    m_dataset = HDF5Handle(H5Dopen(group, data_set_name.c_str()),
	                       &H5Dclose, message.c_str());
#else
    m_dataset = HDF5Handle(H5Dopen(group, data_set_name.c_str(), H5P_DEFAULT),
	                       &H5Dclose, message.c_str());
#endif
    
    HDF5Handle dspace(H5Dget_space(m_dataset), &H5Sclose, "HDF5ImportInfo(): Unable to open data space.");
    
    m_dimensions = H5Sget_simple_extent_ndims(dspace);
    vigra_precondition(m_dimensions >= 2,
             "HDF5ImportInfo(): Number of dimensions is lower than 2. Not an image!");
    
    m_dims.resize(m_dimensions);
    H5Sget_simple_extent_dims(dspace, m_dims.begin(), 0);
    
    HDF5Handle dtype(H5Dget_type(m_dataset), &H5Tclose, "HDF5ImportInfo(): Unable to retrieve data type");

    HDF5Handle native_type(H5Tget_native_type(dtype, H5T_DIR_ASCEND), &H5Tclose, "HDF5ImportInfo(): Unable to retrieve data type");

    if(H5Tequal(native_type, H5T_NATIVE_LLONG) && sizeof(long long) == 4)
        m_pixeltype = "INT32";
    if(H5Tequal(native_type, H5T_NATIVE_LLONG) && sizeof(long long) == 8)
        m_pixeltype = "INT64";
    if(H5Tequal(native_type, H5T_NATIVE_LONG) && sizeof(long) == 4)
        m_pixeltype = "INT32";
    if(H5Tequal(native_type, H5T_NATIVE_LONG) && sizeof(long) == 8)
        m_pixeltype = "INT64";
    if(H5Tequal(native_type, H5T_NATIVE_INT) && sizeof(int) == 4)
        m_pixeltype = "INT32";
    if(H5Tequal(native_type, H5T_NATIVE_INT) && sizeof(int) == 8)
        m_pixeltype = "INT64";
    if(H5Tequal(native_type, H5T_NATIVE_SHORT))
        m_pixeltype = "INT16";
    if(H5Tequal(native_type, H5T_NATIVE_CHAR))
        m_pixeltype = "INT8";

    if(H5Tequal(native_type, H5T_NATIVE_ULLONG) && sizeof(unsigned long long) == 4)
        m_pixeltype = "UINT32";
    if(H5Tequal(native_type, H5T_NATIVE_ULLONG) && sizeof(unsigned long long) == 8)
        m_pixeltype = "UINT64";
    if(H5Tequal(native_type, H5T_NATIVE_ULONG) && sizeof(unsigned long) == 4)
        m_pixeltype = "UINT32";
    if(H5Tequal(native_type, H5T_NATIVE_ULONG) && sizeof(unsigned long) == 8)
        m_pixeltype = "UINT64";
    if(H5Tequal(native_type, H5T_NATIVE_UINT) && sizeof(unsigned int) == 4)
        m_pixeltype = "UINT32";
    if(H5Tequal(native_type, H5T_NATIVE_UINT) && sizeof(unsigned int) == 8)
        m_pixeltype = "UINT64";
    if(H5Tequal(native_type, H5T_NATIVE_USHORT))
        m_pixeltype = "UINT16";
    if(H5Tequal(native_type, H5T_NATIVE_UCHAR))
        m_pixeltype = "UINT8";

    if(H5Tequal(native_type, H5T_NATIVE_LDOUBLE))
        m_pixeltype = "LDOUBLE";
    if(H5Tequal(native_type, H5T_NATIVE_DOUBLE))
        m_pixeltype = "DOUBLE";
    if(H5Tequal(native_type, H5T_NATIVE_FLOAT))
        m_pixeltype = "FLOAT";
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

MultiArrayIndex HDF5ImportInfo::shapeOfDimension(const int dim) const { return (MultiArrayIndex)m_dims[dim]; };
MultiArrayIndex HDF5ImportInfo::numDimensions() const { return m_dimensions; }
const std::string & HDF5ImportInfo::getDatasetName() const { return m_datasetname; }
const std::string & HDF5ImportInfo::getFileName() const { return m_filename; }
hid_t HDF5ImportInfo::getH5FileHandle() const { return m_file.get(); }
hid_t HDF5ImportInfo::getDatasetHandle() const { return m_dataset.get(); }

} // namespace vigra

using namespace vigra;


int main (void)
{
    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try
    {
#if 0
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
#endif
    hid_t file = H5Fcreate( "test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    MultiArray<3, double> a(MultiArrayShape<3>::type(3,4,5));
    for(int k=0; k<a.size(); ++k)
        a[k] = k;
    writeToHDF5File(file, "array", a);
    MultiArray<3, int> aa(MultiArrayShape<3>::type(4,3,2));
    for(int k=0; k<aa.size(); ++k)
        aa[k] = k+100;
    writeToHDF5File(file, "l1/l2/array2", aa);
    H5Fclose(file);
    
    {
        MultiArray<3, double> b;
        HDF5ImportInfo info("test.h5", "array");
        loadFromHDF5File(info, b);
        std::cerr << "data type in group 'array': " << info.getPixelType() << "\n";
        std::cerr << "shape in group 'array': " << info.shape() << "\n";
        std::cerr << b(1,0,0) << "\n";
    }
    {
        MultiArray<3, float> b;
        HDF5ImportInfo info("test.h5", "l1/l2/array2");
        loadFromHDF5File(info, b);
        std::cerr << "data type in group 'l1/l2/array2': " << info.getPixelType() << "\n";
        std::cerr << "shape in group 'l1/l2/array2': " << info.shape() << "\n";
        std::cerr << b(1,0,0) << "\n";
    }    
   }  // end of try block

   // catch failure caused by the H5File operations
   catch( std::exception & e )
   {
        std::cerr << e.what() << "\n";
        return -1;
   }

   return 0;
}


