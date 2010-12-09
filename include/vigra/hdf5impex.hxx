/************************************************************************/
/*                                                                      */
/*       Copyright 2009 by Michael Hanselmann and Ullrich Koethe        */
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

#ifndef VIGRA_HDF5IMPEX_HXX
#define VIGRA_HDF5IMPEX_HXX

#include <string>

#define H5Gcreate_vers 2
#define H5Gopen_vers 2
#define H5Dopen_vers 2
#define H5Dcreate_vers 2
#define H5Acreate_vers 2

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

/** \addtogroup VigraHDF5Impex Import/Export of Images and Arrays in HDF5 Format

    Supports arrays with arbitrary element types and arbitrary many dimensions.
    See the <a href="http://www.hdfgroup.org/HDF5/">HDF5 Website</a> for more
    information on the HDF5 file format.
*/
//@{

    /** \brief Wrapper for hid_t objects.

    Newly created or opened HDF5 handles are usually stored as objects of type 'hid_t'. When the handle
    is no longer needed, the appropriate close function must be called. However, if a function is 
    aborted by an exception, this is difficult to ensure. Class HDF5Handle is a smart pointer that 
    solves this problem by calling the close function in the destructor (This is analogous to how 
    std::auto_ptr calls 'delete' on the contained pointer). A pointer to the close function must be 
    passed to the constructor, along with an error message that is raised when creation/opening fails. 
    
    Since HDF5Handle objects are convertible to hid_t, they can be used in the code in place 
    of the latter.

    <b>Usage:</b>

    \code
    HDF5Handle file_id(H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT), 
                       &H5Fclose, 
                       "Error message.");
                       
    ... // use file_id in the same way as a plain hid_t object
    \endcode

    <b>\#include</b> \<vigra/hdf5impex.hxx\><br>
    Namespace: vigra
    */
class HDF5Handle
{
public:
    typedef herr_t (*Destructor)(hid_t);
    
private:
    hid_t handle_;
    Destructor destructor_;
    
public:

        /** \brief Default constuctor.
            Creates a NULL handle.
        **/
    HDF5Handle()
    : handle_( 0 ),
      destructor_(0)
    {}

        /** \brief Create a wrapper for a hid_t object.

        The hid_t object \a h is assumed to be the return value of an open or create function.
        It will be closed with the given close function \a destructor as soon as this 
        HDF5Handle is destructed, except when \a destructor is a NULL pointer (in which
        case nothing happens at destruction time). If \a h has a value that indicates
        failed opening or creation (by HDF5 convention, this means if it is a negative number),
        an exception is raised by calling <tt>vigra_fail(error_message)</tt>.

        <b>Usage:</b>

        \code
        HDF5Handle file_id(H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT), 
                           &H5Fclose, 
                           "Error message.");
                           
        ... // use file_id in the same way
        \endcode
        */
    HDF5Handle(hid_t h, Destructor destructor, const char * error_message)
    : handle_( h ),
      destructor_(destructor)
    {
        if(handle_ < 0)
            vigra_fail(error_message);
    }

        /** \brief Copy constructor.
            Hands over ownership of the RHS handle (analogous to std::auto_ptr).
        */
    HDF5Handle(HDF5Handle const & h)
    : handle_( h.handle_ ),
      destructor_(h.destructor_)
    {
        const_cast<HDF5Handle &>(h).handle_ = 0;
    }
    
        /** \brief Assignment.
            Calls close() for the LHS handle and hands over ownership of the 
            RHS handle (analogous to std::auto_ptr).
        */
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

        /** \brief Destreuctor.
            Calls close() for the contained handle.
        */
    ~HDF5Handle()
    {
        close();
    }
    
        /** \brief Explicitly call the stored function (if one has been stored within
             this object) for the contained handle and set the handle to NULL.
        */
    herr_t close()
    {
        herr_t res = 1;
        if(handle_ && destructor_)
            res = (*destructor_)(handle_);
        handle_ = 0;
        return res;
    }

        /** \brief Get a temporary hid_t object for the contained handle.
            Do not call a close function on the return value - a crash will be likely
            otherwise.
        */
    hid_t get() const
    {
        return handle_;
    }

        /** \brief Convert to a plain hid_t object.

        This function ensures that hid_t objects can be transparently replaced with 
        HDF5Handle objects in user code. Do not call a close function on the return 
        value - a crash will be likely otherwise.
        */
    operator hid_t() const
    {
        return handle_;
    }

        /** \brief Equality comparison of the contained handle.
        */
    bool operator==(HDF5Handle const & h) const
    {
        return handle_ == h.handle_;
    }

        /** \brief Equality comparison of the contained handle.
        */
    bool operator==(hid_t h) const
    {
        return handle_ == h;
    }

        /** \brief Unequality comparison of the contained handle.
        */
    bool operator!=(HDF5Handle const & h) const
    {
        return handle_ != h.handle_;
    }

        /** \brief Unequality comparison of the contained handle.
        */
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

/** \brief Argument object for the function readHDF5().

See \ref readHDF5() for a usage example. This object must be
used to read an image or array from a HDF5 file 
and enquire about its properties.

<b>\#include</b> \<vigra/hdf5impex.hxx\><br>
Namespace: vigra
*/
class HDF5ImportInfo
{
  public:
    enum PixelType { UINT8, UINT16, UINT32, UINT64, 
	   				 INT8, INT16, INT32, INT64,
					 FLOAT, DOUBLE };

        /** Construct HDF5ImportInfo object.

            The dataset \a pathInFile in the HDF5 file \a filename is accessed to 
            read its properties. \a pathInFile may contain '/'-separated group
            names, but must end with the name of the desired dataset:
            
            \code
            HDF5ImportInfo info(filename, "/group1/group2/my_dataset");
            \endcode
         */
    VIGRA_EXPORT HDF5ImportInfo( const char* filePath, const char* pathInFile );

    VIGRA_EXPORT ~HDF5ImportInfo();

        /** Get the filename of this HDF5 object.
         */
    VIGRA_EXPORT const std::string& getFilePath() const;

        /** Get the dataset's full name in the HDF5 file.
         */
    VIGRA_EXPORT const std::string& getPathInFile() const;

        /** Get a handle to the file represented by this info object.
         */
    VIGRA_EXPORT hid_t getH5FileHandle() const;

        /** Get a handle to the dataset represented by this info object.
         */
    VIGRA_EXPORT hid_t getDatasetHandle() const;

        /** Get the number of dimensions of the dataset represented by this info object.
         */
    VIGRA_EXPORT MultiArrayIndex numDimensions() const;

        /** Get the shape of the dataset represented by this info object.
         */
    VIGRA_EXPORT ArrayVector<hsize_t> const & shape() const
    {
        return m_dims;
    }

        /** Get the shape (length) of the dataset along dimension \a dim.
         */
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
         */
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
         */
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




/********************************************************/
/*                                                      */
/*                     HDF5File                         */
/*                                                      */
/********************************************************/


/** \brief Access to HDF5 files

HDF5File provides a convenient way of accessing data in HDF5 files. vigra::MultiArray
structures of any dimension can be stored to / loaded from HDF5 files. Typical
HDF5 features like subvolume access, chunks and data compression are available,
string attributes can be attached to any dataset or group. Group- or dataset-handles
are encapsulated in the class and managed automatically. The internal file-system like
structure can be accessed by functions like "cd()" or "mkdir()".


<b>Example:</b>
Write the MultiArray out_multi_array to file. Change the current directory to
"/group" and read in the same MultiArray as in_multi_array.
\code
HDF5File file("/path/to/file",HDF5File::New);
file.mkdir("group");
file.write("/group/dataset", out_multi_array);

file.cd("/group");
file.read("dataset", in_multi_array);

\endcode

<b>\#include</b> \<vigra/hdf5impex.hxx\><br>
Namespace: vigra
*/
class HDF5File
{
  private:
    HDF5Handle fileHandle_;

    // current group handle
    HDF5Handle cGroupHandle_;


    // datastructure to hold a list of dataset and group names
    struct lsOpData
    {
        std::vector<std::string> objects;
    };

    // operator function, used by H5Literate
    static herr_t opFunc (hid_t loc_id, const char *name, const H5L_info_t *info, void *operator_data)
    {
        // get information about object
        H5O_info_t      infobuf;
        H5Oget_info_by_name (loc_id, name, &infobuf, H5P_DEFAULT);

        // add name to list, if object is a dataset or a group
        if(infobuf.type == H5O_TYPE_GROUP)
        {
            (*(struct lsOpData *) operator_data).objects.push_back(std::string(name)+"/");
        }
        if(infobuf.type == H5O_TYPE_DATASET)
        {
            (*(struct lsOpData *) operator_data).objects.push_back(std::string(name));
        }

        return 0;
    }

  public:
    /** \brief Set how a file is opened.
      OpenMode::New creates a new file. If the file already exists, overwrite it.

      OpenMode::Open opens a file for reading/writing. The file will be created,
      if nescessary.
    */
    enum OpenMode {
        New,           // Create new empty file (existing file will be deleted).
        Open           // Open file. Create if not existing.
    };


    /** \brief Create a HDF5File object.

    Creates or opens HDF5 file at position filename. The current group is set
    to "/".
    */
    HDF5File(std::string filename, OpenMode mode)
    {
        std::string errorMessage = "HDF5File: Could not create file '" + filename + "'.";
        fileHandle_ = HDF5Handle(createFile_(filename, mode), &H5Fclose, errorMessage.c_str());
        cGroupHandle_ = HDF5Handle(openCreateGroup_("/"), &H5Gclose, "HDF5File(): Failed to open root group.");
    }


    /** \brief Destructor to make sure that all data is flushed before closing the file.
     */
    ~HDF5File()
    {
        //Write everything to disk before closing
        H5Fflush(fileHandle_, H5F_SCOPE_GLOBAL);
    }


    /** \brief Change current group to "/".
     */
    void root()
    {
        std::string message = "HDF5File::root(): Could not open group '/'.";
        cGroupHandle_ = HDF5Handle(H5Gopen(fileHandle_, "/", H5P_DEFAULT),&H5Gclose,message.c_str());
    }


    /** \brief Change the current group.
      If the first character is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    void cd(std::string groupName)
    {
        std::string message = "HDF5File::cd(): Could not open group '" + groupName + "'.\n";
        if(groupName == "/")
        {
            cGroupHandle_ = HDF5Handle(openCreateGroup_("/"),&H5Gclose,message.c_str());
            return;
        }
        else if(groupName =="..")
        {
            cd_up();
        }
        else{
            if(relativePath_(groupName))
            {
                if (H5Lexists(cGroupHandle_, groupName.c_str(), H5P_DEFAULT) == 0)
                {
                    std::cerr << message;
                    return;
                }
                cGroupHandle_ = HDF5Handle(openCreateGroup_(groupName),&H5Gclose,message.c_str());
            }
            else
            {
                if (H5Lexists(fileHandle_, groupName.c_str(), H5P_DEFAULT) == 0)
                {
                    std::cerr << message;
                    return;
                }
                cGroupHandle_ = HDF5Handle(openCreateGroup_(groupName),&H5Gclose,message.c_str());
            }
        }
    }

    /** \brief Change the current group to its parent group.
      returns true if successful, false otherwise.
     */
    bool cd_up()
    {
        std::string groupName = currentGroupName_();

        //do not try to move up if we already in "/"
        if(groupName == "/"){
            return false;
        }

        size_t lastSlash = groupName.find_last_of('/');

        std::string parentGroup (groupName.begin(), groupName.begin()+lastSlash+1);

        cd(parentGroup);

        return true;
    }
    void cd_up(int levels)
    {
        for(int i = 0; i<levels; i++)
            cd_up();
    }


    /** \brief Create a new group.
      If the first character is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    void mkdir(std::string groupName)
    {
        hid_t handle = openCreateGroup_(groupName.c_str());
        if (handle != cGroupHandle_){
            H5Gclose(handle);
        }
    }

    /** \brief Change the current group; create it if nescessary.
      If the first character is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    void cd_mk(std::string groupName)
    {
        std::string  message = "HDF5File::cd_mk(): Could not create group '" + groupName + "'.";
        cGroupHandle_ = HDF5Handle(openCreateGroup_(groupName.c_str()),&H5Gclose,message.c_str());
    }



    /** \brief List the content of the current group.
      The function returns a vector of strings holding the entries of the current
      group. Only datasets and groups are listed, other objects (e.g. datatypes)
      are ignored. Group names always have an ending "/".
     */
    std::vector<std::string> ls()
    {
        lsOpData data;
        H5Literate(cGroupHandle_,H5_INDEX_NAME,H5_ITER_NATIVE,NULL, &opFunc, (void *) &data);

        return data.objects;
    }


    /** \brief Get the path of the current group.
     */
    std::string pwd()
    {
        return currentGroupName_();
    }

    /** \brief Get the name of the associated file.
     */
    std::string filename()
    {
        return fileName_();
    }

    /** \brief Get the number of dimensions of a certain dataset
      If the first character is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    hssize_t getDatasetDimensions(std::string datasetName)
    {
        //Open dataset and dataspace
        std::string errorMessage = "HDF5File::getDatasetDimensions(): Unable to open dataset '" + datasetName + "'.";
        HDF5Handle datasetHandle = HDF5Handle(getDatasetHandle_(datasetName), &H5Dclose, errorMessage.c_str());

        errorMessage = "HDF5File::getDatasetDimensions(): Unable to access dataspace.";
        HDF5Handle dataspaceHandle(H5Dget_space(datasetHandle), &H5Sclose, errorMessage.c_str());

        //return dimension information
        return H5Sget_simple_extent_ndims(dataspaceHandle);
    }

    /** \brief Get the shape of each dimension of a certain dataset.
      Normally, this function is called after determining the dimension of the
      dataset using \ref getDatasetDimensions().
      If the first character is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    ArrayVector<hsize_t> getDatasetShape(std::string datasetName)
    {
        //Open dataset and dataspace
        std::string errorMessage = "HDF5File::getDatasetShape(): Unable to open dataset '" + datasetName + "'.";
        HDF5Handle datasetHandle = HDF5Handle(getDatasetHandle_(datasetName), &H5Dclose, errorMessage.c_str());

        errorMessage = "HDF5File::getDatasetShape(): Unable to access dataspace.";
        HDF5Handle dataspaceHandle(H5Dget_space(datasetHandle), &H5Sclose, errorMessage.c_str());

        //get dimension information
		ArrayVector<hsize_t>::size_type dimensions = H5Sget_simple_extent_ndims(dataspaceHandle);

        ArrayVector<hsize_t> shape(dimensions);
        ArrayVector<hsize_t> maxdims(dimensions);
        H5Sget_simple_extent_dims(dataspaceHandle, shape.data(), maxdims.data());


        // invert the dimensions to guarantee c-order
        ArrayVector<hsize_t> shape_inv(dimensions);
        for(ArrayVector<hsize_t>::size_type i=0; i<shape.size(); i++) {
            shape_inv[i] = shape[dimensions-1-i];
        }

        return shape_inv;
    }

    /** \brief Attach a string attribute to an existing object.
      The attribute can be attached to datasets and groups. The string may have arbitrary length.
     */
    void setAttribute(std::string datasetName, std::string attributeName, std::string text)
    {
        std::string groupname = SplitString(datasetName).first();
        std::string setname = SplitString(datasetName).last();

        std::string errorMessage ("HDF5File::setAttribute(): Unable to open group '" + groupname + "'.");
        HDF5Handle groupHandle (openCreateGroup_(groupname), &H5Gclose, errorMessage.c_str());

        H5LTset_attribute_string(groupHandle,setname.c_str(), attributeName.c_str(),text.c_str());
    }


    /** \brief Get an attribute string of an object.
     */
    std::string getAttribute(std::string datasetName, std::string attributeName)
    {
        typedef ArrayVector<char>::size_type SizeType;
		if (H5Lexists(fileHandle_, datasetName.c_str(), H5P_DEFAULT) == 0)
        {
            std::cerr << "HDF5File::getAttribute(): Dataset '" << datasetName << "' does not exist.\n";
            return "error - dataset not found";
        }

        std::string groupname = SplitString(datasetName).first();
        std::string setname = SplitString(datasetName).last();

        std::string errorMessage ("HDF5File::setAttribute(): Unable to open group '" + groupname + "'.");
        HDF5Handle groupHandle (openCreateGroup_(groupname), &H5Gclose, errorMessage.c_str());

        // get the size of the attribute
        HDF5Handle AttrHandle (H5Aopen_by_name(groupHandle,setname.c_str(),attributeName.c_str(),H5P_DEFAULT, H5P_DEFAULT),&H5Aclose, "HDF5File::getAttribute(): Unable to open attribute.");
		SizeType len = (SizeType)H5Aget_storage_size(AttrHandle);

        //read the attribute
        ArrayVector<char> text (len+1,0);
        H5LTget_attribute_string(groupHandle, setname.c_str(), attributeName.c_str(), text.begin());

        return std::string(text.begin());
    }


    // Writing data

    /** \brief Write multi arrays.
      
      Chunks can be activated by setting 
      \code iChunkSize = size; //size \> 0 
      \endcode .
      The chunks will be hypercubes with edge length size.

      Compression can be activated by setting 
      \code compression = parameter; // 0 \< parameter \<= 9 
      \endcode
      where 0 stands for no compression and 9 for maximum compression.

      If the first character of datasetName is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    template<unsigned int N, class T>
    inline void write(std::string datasetName, const MultiArrayView<N, T, UnstridedArrayTag> & array, int iChunkSize = 0, int compression = 0)
    {
        typename MultiArrayShape<N>::type chunkSize;
        for(int i = 0; i < N; i++){
            chunkSize[i] = iChunkSize;
        }
        write_(datasetName, array, detail::getH5DataType<T>(), 1, chunkSize, compression);
    }

    /** \brief Write multi arrays.
      Chunks can be activated by providing a MultiArrayShape as chunkSize.
      chunkSize must have equal dimension as array.

      Compression can be activated by setting 
      \code compression = parameter; // 0 \< parameter \<= 9 
      \endcode
      where 0 stands for no compression and 9 for maximum compression.

      If the first character of datasetName is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    template<unsigned int N, class T>
    inline void write(std::string datasetName, const MultiArrayView<N, T, UnstridedArrayTag> & array, typename MultiArrayShape<N>::type chunkSize, int compression = 0)
    {
        write_(datasetName, array, detail::getH5DataType<T>(), 1, chunkSize, compression);
    }

    /** \brief Write a multi array into a larger volume.
      blockOffset determines the position, where array is written.

      Chunks can be activated by providing a MultiArrayShape as chunkSize.
      chunkSize must have equal dimension as array.

      Compression can be activated by setting 
      \code compression = parameter; // 0 \< parameter \<= 9 
      \endcode
      where 0 stands for no compression and 9 for maximum compression.

      If the first character of datasetName is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    template<unsigned int N, class T>
    inline void writeBlock(std::string datasetName, typename MultiArrayShape<N>::type blockOffset, const MultiArrayView<N, T, UnstridedArrayTag> & array)
    {
        writeBlock_(datasetName, blockOffset, array, detail::getH5DataType<T>(), 1);
    }

    // non-scalar (TinyVector) and unstrided multi arrays
    template<unsigned int N, class T, int SIZE>
    inline void write(std::string datasetName, const MultiArrayView<N, TinyVector<T, SIZE>, UnstridedArrayTag> & array, int iChunkSize = 0, int compression = 0)
    {
        typename MultiArrayShape<N>::type chunkSize;
        for(int i = 0; i < N; i++){
            chunkSize[i] = iChunkSize;
        }
        write_(datasetName, array, detail::getH5DataType<T>(), SIZE, chunkSize, compression);
    }

    template<unsigned int N, class T, int SIZE>
    inline void write(std::string datasetName, const MultiArrayView<N, TinyVector<T, SIZE>, UnstridedArrayTag> & array, typename MultiArrayShape<N>::type chunkSize, int compression = 0)
    {
        write_(datasetName, array, detail::getH5DataType<T>(), SIZE, chunkSize, compression);
    }

    template<unsigned int N, class T, int SIZE>
    inline void writeBlock(std::string datasetName, typename MultiArrayShape<N>::type blockOffset, const MultiArrayView<N, TinyVector<T, SIZE>, UnstridedArrayTag> & array)
    {
        writeBlock_(datasetName, blockOffset, array, detail::getH5DataType<T>(), SIZE);
    }

    // non-scalar (RGBValue) and unstrided multi arrays
    template<unsigned int N, class T>
    inline void write(std::string datasetName, const MultiArrayView<N, RGBValue<T>, UnstridedArrayTag> & array, int iChunkSize = 0, int compression = 0)
    {
        typename MultiArrayShape<N>::type chunkSize;
        for(int i = 0; i < N; i++){
            chunkSize[i] = iChunkSize;
        }
        write_(datasetName, array, detail::getH5DataType<T>(), 3, chunkSize, compression);
    }

    template<unsigned int N, class T>
    inline void write(std::string datasetName, const MultiArrayView<N, RGBValue<T>, UnstridedArrayTag> & array, typename MultiArrayShape<N>::type chunkSize, int compression = 0)
    {
        write_(datasetName, array, detail::getH5DataType<T>(), 3, chunkSize, compression);
    }

    template<unsigned int N, class T>
    inline void writeBlock(std::string datasetName, typename MultiArrayShape<N>::type blockOffset, const MultiArrayView<N, RGBValue<T>, UnstridedArrayTag> & array)
    {
        writeBlock_(datasetName, blockOffset, array, detail::getH5DataType<T>(), 3);
    }

    /** \brief Write single value as dataset.
      This functions allows to write data of atomic datatypes (int, long, double)
      as a dataset in the HDF5 file. So it is not nescessary to create a MultiArray
      of size 1 to write a single number.

      If the first character of datasetName is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
   template<class T>
    inline void writeAtomic(std::string datasetName, const T data)
    {
        typename MultiArrayShape<1>::type chunkSize;
        chunkSize[0] = 0;
        MultiArray<1,T> array(MultiArrayShape<1>::type(1));
        array[0] = data;
        write_(datasetName, array, detail::getH5DataType<T>(), 1, chunkSize,0);
    }


    /* Reading data. */

    /** \brief Read data into a multi array.
      If the first character of datasetName is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    template<unsigned int N, class T>
    inline void read(std::string datasetName, MultiArrayView<N, T, UnstridedArrayTag> array)
    {
        read_(datasetName, array, detail::getH5DataType<T>(), 1);
    }

    /** \brief Read a block of data into s multi array.
      This function allows to read a small block out of a larger volume stored
      in a HDF5 dataset.

      blockOffset determines the position of the block.
      blockSize determines the size in each dimension of the block.

      If the first character of datasetName is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    template<unsigned int N, class T>
    inline void readBlock(std::string datasetName, typename MultiArrayShape<N>::type blockOffset, typename MultiArrayShape<N>::type blockShape, MultiArrayView<N, T, UnstridedArrayTag> array)
    {
        readBlock_(datasetName, blockOffset, blockShape, array, detail::getH5DataType<T>(), 1);
    }

    // non-scalar (TinyVector) and unstrided target multi array
    template<unsigned int N, class T, int SIZE>
    inline void read(std::string datasetName, MultiArrayView<N, TinyVector<T, SIZE>, UnstridedArrayTag> array)
    {
        read_(datasetName, array, detail::getH5DataType<T>(), SIZE);
    }

    template<unsigned int N, class T, int SIZE>
    inline void readBlock(std::string datasetName, typename MultiArrayShape<N>::type blockOffset, typename MultiArrayShape<N>::type blockShape, MultiArrayView<N, TinyVector<T, SIZE>, UnstridedArrayTag> array)
    {
        readBlock_(datasetName, blockOffset, blockShape, array, detail::getH5DataType<T>(), SIZE);
    }

    // non-scalar (RGBValue) and unstrided target multi array
    template<unsigned int N, class T>
    inline void read(std::string datasetName, MultiArrayView<N, RGBValue<T>, UnstridedArrayTag> array)
    {
        read_(datasetName, array, detail::getH5DataType<T>(), 3);
    }

    template<unsigned int N, class T>
    inline void readBlock(std::string datasetName, typename MultiArrayShape<N>::type blockOffset, typename MultiArrayShape<N>::type blockShape, MultiArrayView<N, RGBValue<T>, UnstridedArrayTag> array)
    {
        readBlock_(datasetName, blockOffset, blockShape, array, detail::getH5DataType<T>(), 3);
    }

    /** \brief Read a single value.
      This functions allows to read a single datum of atomic datatype (int, long, double)
      from the HDF5 file. So it is not nescessary to create a MultiArray
      of size 1 to read a single number.

      If the first character of datasetName is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    template<class T>
    inline void readAtomic(std::string datasetName, T & data)
    {
        MultiArray<1,T> array(MultiArrayShape<1>::type(1));
        read_(datasetName, array, detail::getH5DataType<T>(), 1);
        data = array[0];
    }


    /** \brief Create a new dataset.
      This function can be used to create a dataset filled with a default value,
      for example before writing data into it using \ref writeBlock().
      Attention: only atomic datatypes are provided. For spectral data, add an
      dimension (case RGB: add one dimension of size 3).

      shape determines the dimension and the size of the dataset.

      Chunks can be activated by providing a MultiArrayShape as chunkSize.
      chunkSize must have equal dimension as array.

      Compression can be activated by setting 
      \code compression = parameter; // 0 \< parameter \<= 9 
      \endcode
      where 0 stands for no compression and 9 for maximum compression.

      If the first character of datasetName is a "/", the path will be interpreted as absolute path,
      otherwise it will be interpreted as path relative to the current group.
     */
    template<unsigned int N, class T>
    inline void createDataset(std::string datasetName, typename MultiArrayShape<N>::type shape, T init, int iChunkSize = 0, int compressionParameter = 0)
    {
        typename MultiArrayShape<N>::type chunkSize;
        for(int i = 0; i < N; i++){
            chunkSize[i] = iChunkSize;
        }
        createDataset<N,T>(datasetName, shape, init, chunkSize, compressionParameter);
    }

    template<unsigned int N, class T>
    inline void createDataset(std::string datasetName, typename MultiArrayShape<N>::type shape, T init, typename MultiArrayShape<N>::type chunkSize, int compressionParameter = 0)
    {
        std::string groupname = SplitString(datasetName).first();
        std::string setname = SplitString(datasetName).last();

        hid_t parent = openCreateGroup_(groupname);

        // delete the dataset if it already exists
        deleteDataset_(parent, setname);

        // create dataspace
        // add an extra dimension in case that the data is non-scalar
        HDF5Handle dataspaceHandle;

        // invert dimensions to guarantee c-order
        hsize_t shape_inv[N];
        for(unsigned int k=0; k<N; ++k)
            shape_inv[N-1-k] = shape[k];

        // create dataspace
        dataspaceHandle = HDF5Handle(H5Screate_simple(N, shape_inv, NULL),
                                    &H5Sclose, "HDF5File::createDataset(): unable to create dataspace for scalar data.");

        // set fill value
        HDF5Handle plist ( H5Pcreate(H5P_DATASET_CREATE), &H5Pclose, "HDF5File::createDataset(): unable to create property list." );
        H5Pset_fill_value(plist,detail::getH5DataType<T>(), &init);

        // enable chunks
        if(chunkSize[0] > 0)
        {
            hsize_t cSize [N];
            for(int i = 0; i<N; i++)
            {
                cSize[i] = chunkSize[N-1-i];
            }
            H5Pset_chunk (plist, N, cSize);
        }

        // enable compression
        if(compressionParameter > 0)
        {
            H5Pset_deflate(plist, compressionParameter);
        }

        //create the dataset.
        HDF5Handle datasetHandle ( H5Dcreate(parent, setname.c_str(), detail::getH5DataType<T>(), dataspaceHandle, H5P_DEFAULT, plist, H5P_DEFAULT),
                                  &H5Dclose, "HDF5File::createDataset(): unable to create dataset.");
        if(parent != cGroupHandle_)
            H5Gclose(parent);
    }

    /** \brief Immediately write all data to disk
     */
    void flushToDisk()
    {
        H5Fflush(fileHandle_, H5F_SCOPE_GLOBAL);
    }


  private:

    /* Simple extension of std::string for splitting into two parts
     *
     *  Strings (in particular: file/dataset paths) will be split into two
     *  parts. The split is made at the last occurance of the delimiter.
     *
     *  For example, "/path/to/some/file" will be split (delimiter = "/") into
     *  first() = "/path/to/some" and last() = "file".
     */
    class SplitString: public std::string {
    public:
        SplitString(std::string &sstring): std::string(sstring) {};

        // return the part of the string before the delimiter
        std::string first(char delimiter = '/')
        {
            size_t last = find_last_of(delimiter);
            if(last == std::string::npos) // delimiter not found --> no first
                return "";

            return std::string(begin(), begin()+last+1);
        }

        // return the part of the string after the delimiter
        std::string last(char delimiter = '/')
        {
            size_t last = find_last_of(delimiter);
            if(last == std::string::npos) // delimiter not found --> only last
                return std::string(*this);
            return std::string(begin()+last+1, end());
        }
    };

    inline bool relativePath_(std::string & path)
    {
        std::string::size_type pos = path.find('/') ;
        if(pos == 0)
            return false;

        return true;
    }


    std::string currentGroupName_()
    {
        int len = H5Iget_name(cGroupHandle_,NULL,1000);
        ArrayVector<char> name (len+1,0);
        H5Iget_name(cGroupHandle_,name.begin(),len+1);

        return std::string(name.begin());
    }

    std::string fileName_()
    {
        int len = H5Fget_name(fileHandle_,NULL,1000);
        ArrayVector<char> name (len+1,0);
        H5Fget_name(fileHandle_,name.begin(),len+1);

        return std::string(name.begin());
    }


    inline hid_t createFile_(std::string filePath, OpenMode mode = Open)
    {
        // try to open file
        FILE * pFile;
        pFile = fopen ( filePath.c_str(), "r" );
        hid_t fileId;

        // check if opening was successful (= file exists)
        if ( pFile == NULL )
        {
            fileId = H5Fcreate(filePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
        else if(mode == Open)
        {
            fclose( pFile );
            fileId = H5Fopen(filePath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        }
        else
        {
            fclose(pFile);
            std::remove(filePath.c_str());
            fileId = H5Fcreate(filePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
        return fileId;
    }


    // open a group and subgroups. Create if nescessary.
    inline hid_t openCreateGroup_(std::string groupName)
    {
        // check if groupName is absolute or relative
        hid_t parent;
        if(relativePath_(groupName))
        {
            parent = cGroupHandle_;
        }
        else
        {
            // open root group
            parent = H5Gopen(fileHandle_, "/", H5P_DEFAULT);
            if(groupName == "/")
            {
                return parent;
            }

            // remove leading /
            groupName = std::string(groupName.begin()+1, groupName.end());
        }


        // check if the groupName has the right format
        if( groupName.size() != 0 && *groupName.rbegin() != '/')
        {
            groupName = groupName + '/';
        }


        //create subgroups one by one
        std::string::size_type begin = 0, end = groupName.find('/');
        int ii =  0;
        while (end != std::string::npos)
        {
            std::string group(groupName.begin()+begin, groupName.begin()+end);
            hid_t prevParent = parent;

            if(H5LTfind_dataset(parent, group.c_str()) == 0)
            {
                parent = H5Gcreate(prevParent, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            } else {
                parent = H5Gopen(prevParent, group.c_str(), H5P_DEFAULT);
            }

            if(ii != 0)
            {
                H5Gclose(prevParent);
            }
            if(parent < 0)
            {
                return parent;
            }
            ++ii;
            begin = end + 1;
            end = groupName.find('/', begin);
        }


        return parent;
    }


    inline void deleteDataset_(hid_t parent, std::string datasetName)
    {
        // delete existing data and create new dataset
        if(H5LTfind_dataset(parent, datasetName.c_str()))
        {
            //std::cout << "dataset already exists" << std::endl;
    #if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
            if(H5Gunlink(parent, datasetName.c_str()) < 0)
            {
                vigra_postcondition(false, "HDF5File::deleteDataset_(): Unable to delete existing data.");
            }
    #else
            if(H5Ldelete(parent, datasetName.c_str(), H5P_DEFAULT ) < 0)
            {
                vigra_postcondition(false, "HDF5File::deleteDataset_(): Unable to delete existing data.");
            }
    #endif
        }
    }

    hid_t getDatasetHandle_(std::string datasetName)
    {
        std::string groupname = SplitString(datasetName).first();
        std::string setname = SplitString(datasetName).last();

        if(relativePath_(datasetName))
        {
            if (H5Lexists(cGroupHandle_, datasetName.c_str(), H5P_DEFAULT) != 1)
            {
                std::cerr << "HDF5File::getDatasetHandle_(): Dataset '" << datasetName << "' does not exist.\n";
                return -1;
            }
            //Open parent group
            hid_t groupHandle = openCreateGroup_(groupname);

            hid_t datasetHandle = H5Dopen(groupHandle, setname.c_str(), H5P_DEFAULT);

            if(groupHandle != cGroupHandle_){
                H5Gclose(groupHandle);
            }

            //return dataset handle
            return datasetHandle;
        }
        else
        {
            if (H5Lexists(fileHandle_, datasetName.c_str(), H5P_DEFAULT) != 1)
            {
                std::cerr << "HDF5File::getDatasetHandle_(): Dataset '" << datasetName << "' does not exist.\n";
                return -1;
            }

            //Open parent group
            hid_t groupHandle = openCreateGroup_(groupname);

            hid_t datasetHandle = H5Dopen(groupHandle, setname.c_str(), H5P_DEFAULT);

            if(groupHandle != cGroupHandle_){
                H5Gclose(groupHandle);
            }

            //return dataset handle
            return datasetHandle;
        }
    }


    // unstrided multi arrays
    template<unsigned int N, class T>
    void write_(std::string &datasetName, const MultiArrayView<N, T, UnstridedArrayTag> & array, const hid_t datatype, const int numBandsOfType, typename MultiArrayShape<N>::type &chunkSize, int compressionParameter = 0)
    {
        std::string groupname = SplitString(datasetName).first();
        std::string setname = SplitString(datasetName).last();

        // shape of the array. Add one dimension, if array contains non-scalars.
        ArrayVector<hsize_t> shape(N + (numBandsOfType > 1),0);
        for(int i = 0; i < N; i++){
            shape[N-1-i] = array.shape(i); // reverse order
        }

        if(numBandsOfType > 1)
            shape[N] = numBandsOfType;

        HDF5Handle dataspace ( H5Screate_simple(N + (numBandsOfType > 1), shape.begin(), NULL), &H5Sclose, "HDF5File::write(): Can not create dataspace.");

        // create and open group:
        std::string errorMessage ("HDF5File::write(): can not create group '" + groupname + "'.");
        hid_t groupHandle = openCreateGroup_(groupname);
        if(groupHandle <= 0)
        {
            std::cerr << errorMessage << "\n";
        }

        // delete dataset, if it already exists
        deleteDataset_(groupHandle, setname.c_str());

        // set up properties list
        HDF5Handle plist ( H5Pcreate(H5P_DATASET_CREATE), &H5Pclose, "HDF5File::write(): unable to create property list." );

        // enable chunks
        if(chunkSize[0] > 0)
        {
            ArrayVector<hsize_t> cSize(N + (numBandsOfType > 1),0);
            for(int i = 0; i<N; i++)
            {
                cSize[i] = chunkSize[N-1-i];
            }
            if(numBandsOfType > 1)
                cSize[N] = numBandsOfType;
            
            H5Pset_chunk (plist, N + (numBandsOfType > 1), cSize.begin());
        }

        // enable compression
        if(compressionParameter > 0)
        {
            H5Pset_deflate(plist, compressionParameter);
        }

        // create dataset
        HDF5Handle datasetHandle (H5Dcreate(groupHandle, setname.c_str(), datatype, dataspace,H5P_DEFAULT, plist, H5P_DEFAULT), &H5Dclose, "HDF5File::write(): Can not create dataset.");

        // Write the data to the HDF5 dataset as is
        H5Dwrite( datasetHandle, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data());

        if(groupHandle != cGroupHandle_)
        {
            H5Gclose(groupHandle);
        }
    }

    // unstrided target multi array
    template<unsigned int N, class T>
    void read_(std::string datasetName, MultiArrayView<N, T, UnstridedArrayTag> array, const hid_t datatype, const int numBandsOfType)
    {
        //Prepare to read without using HDF5ImportInfo
        ArrayVector<hsize_t> dimshape = getDatasetShape(datasetName);
        hssize_t dimensions = getDatasetDimensions(datasetName);

        std::string errorMessage ("HDF5File::read(): Unable to open dataset '" + datasetName + "'.");
        HDF5Handle datasetHandle (getDatasetHandle_(datasetName), &H5Dclose, errorMessage.c_str());

        int offset = (numBandsOfType > 1);

        vigra_precondition(( (N + offset ) ==  MultiArrayIndex(dimensions)), // the object in the HDF5 file may have one additional dimension which we then interpret as the pixel type bands
            "HDF5File::read(): Array dimension disagrees with dataset dimension.");

        typename MultiArrayShape<N>::type shape;
        for(int k=offset; k< MultiArrayIndex(dimensions); ++k) {
            shape[k-offset] = MultiArrayIndex(dimshape[k]);
        }

        vigra_precondition(shape == array.shape(),
                           "HDF5File::read(): Array shape disagrees with dataset shape.");

        // simply read in the data as is
        H5Dread( datasetHandle, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data() ); // .data() possible since void pointer!
    }

    template<unsigned int N, class T>
    void writeBlock_(std::string datasetName, typename MultiArrayShape<N>::type &blockOffset, const MultiArrayView<N, T, UnstridedArrayTag> & array, const hid_t datatype, const int numBandsOfType)
    {
        // open dataset if it exists
        std::string errorMessage = "HDF5File::writeBlock(): Error opening dataset '" + datasetName + "'.";
        HDF5Handle datasetHandle (getDatasetHandle_(datasetName), &H5Dclose, errorMessage.c_str());


        // hyperslab parameters for position, size, ...
        hsize_t boffset [N];
        hsize_t bshape [N];
        hsize_t bones [N];

        for(int i = 0; i < N; i++){
            boffset[i] = blockOffset[N-1-i];
            bshape[i] = array.size(N-1-i);
            bones[i] = 1;
        }

        // create a target dataspace in memory with the shape of the desired block
        HDF5Handle memspace_handle (H5Screate_simple(N,bshape,NULL),&H5Sclose,"Unable to get origin dataspace");

        // get file dataspace and select the desired block
        HDF5Handle dataspaceHandle (H5Dget_space(datasetHandle),&H5Sclose,"Unable to create target dataspace");
        H5Sselect_hyperslab(dataspaceHandle, H5S_SELECT_SET, boffset, bones, bones, bshape);

        // Write the data to the HDF5 dataset as is
        H5Dwrite( datasetHandle, datatype, memspace_handle, dataspaceHandle, H5P_DEFAULT, array.data()); // .data() possible since void pointer!
    }



    // unstrided target multi array
    // read a block of a HDF5 dataset into a MultiArray
    template<unsigned int N, class T>
    void readBlock_(std::string datasetName, typename MultiArrayShape<N>::type &blockOffset, typename MultiArrayShape<N>::type &blockShape, MultiArrayView<N, T, UnstridedArrayTag> &array, const hid_t datatype, const int numBandsOfType)
    {
        //Prepare to read without using HDF5ImportInfo
        //ArrayVector<hsize_t> dimshape = getDatasetShape(datasetName) ;
        hssize_t dimensions = getDatasetDimensions(datasetName);

        std::string errorMessage ("HDF5File::readBlock(): Unable to open dataset '" + datasetName + "'.");
        HDF5Handle datasetHandle (getDatasetHandle_(datasetName), &H5Dclose, errorMessage.c_str());


        int offset = (numBandsOfType > 1);

        vigra_precondition(( (N + offset ) ==  MultiArrayIndex(dimensions)), // the object in the HDF5 file may have one additional dimension which we then interpret as the pixel type bands
            "readHDF5_block(): Array dimension disagrees with data dimension.");

        vigra_precondition(blockShape == array.shape(),
             "readHDF5_block(): Array shape disagrees with block size.");

        // hyperslab parameters for position, size, ...
        hsize_t boffset [N];
        hsize_t bshape [N];
        hsize_t bones [N];

        for(int i = 0; i < N; i++){
            // virgra and hdf5 use different indexing
            boffset[i] = blockOffset[N-1-i];
            //bshape[i] = blockShape[i];
            bshape[i] = blockShape[N-1-i];
            //boffset[i] = blockOffset[N-1-i];
            bones[i] = 1;
        }

        // create a target dataspace in memory with the shape of the desired block
        HDF5Handle memspace_handle (H5Screate_simple(N,bshape,NULL),&H5Sclose,"Unable to create target dataspace");

        // get file dataspace and select the desired block
        HDF5Handle dataspaceHandle (H5Dget_space(datasetHandle),&H5Sclose,"Unable to get dataspace");
        H5Sselect_hyperslab(dataspaceHandle, H5S_SELECT_SET, boffset, bones, bones, bshape);

        // now read the data
        H5Dread( datasetHandle, datatype, memspace_handle, dataspaceHandle, H5P_DEFAULT, array.data() ); // .data() possible since void pointer!
    }

};  /* class HDF5File */







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

/** \brief Read the data specified by the given \ref vigra::HDF5ImportInfo object
                and write the into the given 'array'.
                
    The array must have the correct number of dimensions and shape for the dataset 
    represented by 'info'. When the element type of 'array' differs from the stored element
    type, HDF5 will convert the type on the fly (except when the HDF5 version is 1.6 or below,
    in which case an error will result). Multi-channel element types (i.e. \ref vigra::RGBValue
    and \ref vigra::TinyVector) are recognized and handled correctly.
    
    <b> Declaration:</b>
    
    \code
    namespace vigra {
        template<unsigned int N, class T, class StrideTag>
        void 
        readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, T, StrideTag> array);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/hdf5impex.hxx\><br>
    Namespace: vigra
    
    \code
    
    HDF5ImportInfo info(filename, dataset_name);
    vigra_precondition(info.numDimensions() == 3, "Dataset must be 3-dimensional.");
    
    MultiArrayShape<3>::type shape(info.shape().begin());
    MultiArray<3, int> array(shape);
    
    readHDF5(info, array);
    \endcode
*/
doxygen_overloaded_function(template <...> void readHDF5)

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
    for(unsigned int i=0;i<N;++i)
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
writeHDF5Impl(DestIterator d, Shape const & shape, const hid_t dataset_id, const hid_t datatype, ArrayVector<T> & buffer, int & counter, const int elements, const int numBandsOfType, MetaInt<0>)
{
    DestIterator dend = d + (typename DestIterator::difference_type)shape[0];
    int k = 0;
	//std::cout << "new:" << std::endl;
	for(; d < dend; ++d, k++)
    {
        buffer[k] = *d; 
        //std::cout << buffer[k] << " ";
    }
	//std::cout << std::endl;
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
writeHDF5Impl(DestIterator d, Shape const & shape, const hid_t dataset_id, const hid_t datatype, ArrayVector<T> & buffer, int & counter, const int elements, const int numBandsOfType, MetaInt<N>)
{
		DestIterator dend = d + (typename DestIterator::difference_type)shape[N];
		for(; d < dend; ++d)
		{
			writeHDF5Impl(d.begin(), shape, dataset_id, datatype, buffer, counter, elements, numBandsOfType, MetaInt<N-1>());
		}
}

} // namespace detail

/** \brief Store array data in an HDF5 file.
                
    The number of dimensions, shape and element type of the stored dataset is automatically 
    determined from the properties of the given \a array. Strided arrays are stored in an
    unstrided way, i.e. in contiguous scan-order. Multi-channel element types 
    (i.e. \ref vigra::RGBValue and \ref vigra::TinyVector) are recognized and handled correctly
    (in particular, the will form the innermost dimension of the stored dataset).
    \a pathInFile may contain '/'-separated group names, but must end with the name 
    of the dataset to be created.
    
    <b> Declaration:</b>
    
    \code
    namespace vigra {
        template<unsigned int N, class T, class StrideTag>
        void 
        writeHDF5(const char* filePath, const char* pathInFile, 
                  MultiArrayView<N, T, StrideTag>const  & array);
    }
    \endcode
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/hdf5impex.hxx\><br>
    Namespace: vigra
    
    \code
    MultiArrayShape<3>::type shape(100, 200, 20);
    MultiArray<3, int> array(shape);
    ... // fill array with data
    
    writeHDF5("mydata.h5", "/group1/my_dataset", array);
    \endcode
*/
doxygen_overloaded_function(template <...> void writeHDF5)

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
	detail::writeHDF5Impl(array.traverser_begin(), shape, dataset_handle, datatype, buffer, counter, elements, numBandsOfType, vigra::MetaInt<N-1>());

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


#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 8) || DOXYGEN
/** Write a numeric MultiArray as an attribute with name \a name 
    of the dataset specified by the handle \a loc. 

    <b>\#include</b> \<vigra/hdf5impex.hxx\><br>
    Namespace: vigra
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

/** Write a string MultiArray as an attribute with name \a name 
    of the dataset specified by the handle \a loc. 

    <b>\#include</b> \<vigra/hdf5impex.hxx\><br>
    Namespace: vigra
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

/** Write a numeric ArrayVectorView as an attribute with name \a name 
    of the dataset specified by the handle \a loc. 

    <b>\#include</b> \<vigra/hdf5impex.hxx\><br>
    Namespace: vigra
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
    the path in the file should have the format 
    [attribute] or /[subgroups/]dataset.attribute or
    /[subgroups/]group.attribute.
    The attribute is written to the root group, a dataset or a subgroup
    respectively
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

//@}

} // namespace vigra

#endif // VIGRA_HDF5IMPEX_HXX
