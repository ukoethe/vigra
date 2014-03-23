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
#include <algorithm>
#include <utility>

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
# ifndef H5Pset_obj_track_times
#  define H5Pset_obj_track_times(a, b) do {} while (0)
# endif
# include <H5LT.h>
#else
# include <hdf5_hl.h>
#endif

#include "impex.hxx"
#include "multi_array.hxx"
#include "multi_iterator_coupled.hxx"
#include "multi_impex.hxx"
#include "utilities.hxx"
#include "error.hxx"

#if defined(_MSC_VER)
#  include <io.h>
#else
#  include <unistd.h>
#endif

namespace vigra {

/** \addtogroup VigraHDF5Impex Import/Export of Images and Arrays in HDF5 Format

    Supports arrays with arbitrary element types and arbitrary many dimensions.
    See the <a href="http://www.hdfgroup.org/HDF5/">HDF5 Website</a> for more
    information on the HDF5 file format.
*/
//@{

    /** \brief Check if given filename refers to a HDF5 file.
    */
inline bool isHDF5(char const * filename)
{
#ifdef _MSC_VER
    return _access(filename, 0) != -1 && H5Fis_hdf5(filename);
#else
    return access(filename, F_OK) == 0 && H5Fis_hdf5(filename);
#endif
}

    /** \brief Wrapper for unique hid_t objects.

    This class offers the functionality of <tt>std::unique_ptr</tt> for HDF5 handles 
    (type <tt>hid_t</tt>). Unfortunately, <tt>std::unique_ptr</tt> cannot be used directly 
    for this purpose because it only works with pointers, whereas <tt>hid_t</tt> is an integer type.
    
    Newly created or opened HDF5 handles are stored as objects of type <tt>hid_t</tt>. When the handle
    is no longer needed, the appropriate close function must be called. However, if a function is 
    aborted by an exception, this is difficult to ensure. Class HDF5Handle is a smart pointer that 
    solves this problem by calling the close function in the destructor (This is analogous to how 
    <tt>std::unique_ptr</tt> calls 'delete' on the contained pointer). A pointer to the close function 
    must be passed to the constructor, along with an error message that is raised when 
    creation/opening fails. 
    
    When a <tt>HDF5Handle</tt> is created or assigned from another one, ownership passes on to the  
    left-hand-side handle object, and the right-hand-side objects is resest to a NULL handle.
    
    Since <tt>HDF5Handle</tt> objects are convertible to <tt>hid_t</tt>, they can be used in the code 
    in place of the latter.

    <b>Usage:</b>

    \code
    HDF5Handle file_id(H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT), 
                       &H5Fclose, 
                       "Error message when H5Fopen() fails.");
                       
    ... // use file_id in the same way as a plain hid_t object
    
    // pass ownership to a new handle object
    HDF5Handle new_handle(file_id);
    
    assert(file_id.get() == 0);
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

        /** \brief Default constructor.
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
        failed opening or creation (by HDF5 convention, this means that \a h is negative),
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
        
            Hands over ownership of the RHS handle (analogous to <tt>std::unique_pt</tt>).
        */
    HDF5Handle(HDF5Handle const & h)
    : handle_( h.handle_ ),
      destructor_(h.destructor_)
    {
        const_cast<HDF5Handle &>(h).handle_ = 0;
    }
    
        /** \brief Assignment.
            Calls close() for the LHS handle and hands over ownership of the 
            RHS handle (analogous to <tt>std::unique_pt</tt>).
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

        /** \brief Destructor.
            Calls close() for the contained handle.
        */
    ~HDF5Handle()
    {
        close();
    }
    
        /** \brief Explicitly call the stored destructor (if one has been registered in the
             constructor) for the contained handle and set the wrapper to NULL. Returns
             a negative value when the destructor call for the handle fails, and
             a non-negative value otherwise.
        */
    herr_t close()
    {
        herr_t res = 1;
        if(handle_ && destructor_)
            res = (*destructor_)(handle_);
        handle_ = 0;
        destructor_ = 0;
        return res;
    }
    
        /** \brief Return the contained handle and set the wrapper to NULL
            without calling <tt>close()</tt>.
        */
    hid_t release()
    {
        hid_t res = handle_;
        handle_ = 0;
        destructor_ = 0;
        return res;
    }
    
        /** \brief Reset the wrapper to a new handle.
        
             Equivalent to <tt>handle = HDF5Handle(h, destructor, error_message)</tt>.
        */
    void reset(hid_t h, Destructor destructor, const char * error_message)
    {
        if(h < 0)
            vigra_fail(error_message);
        if(h != handle_)
        {
            close();
            handle_ = h;
            destructor_ = destructor;
        }
    }
    
        /** \brief Swap the contents of two handle wrappers.
        
            Also available as <tt>std::swap(handle1, handle2)</tt>.
        */    
    void swap(HDF5Handle & h)
    {
        std::swap(handle_, h.handle_);
        std::swap(destructor_, h.destructor_);
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

        /** \brief Inequality comparison of the contained handle.
        */
    bool operator!=(HDF5Handle const & h) const
    {
        return handle_ != h.handle_;
    }

        /** \brief Inequality comparison of the contained handle.
        */
    bool operator!=(hid_t h) const
    {
        return handle_ != h;
    }
};


    /** \brief Wrapper for shared hid_t objects.

    This class offers the functionality of <tt>std::shared_ptr</tt> for HDF5 handles 
    (type <tt>hid_t</tt>). Unfortunately, <tt>std::shared_ptr</tt> cannot be used directly 
    for this purpose because it only works with pointers, whereas <tt>hid_t</tt> is an integer type.
    
    Newly created or opened HDF5 handles are stored as objects of type <tt>hid_t</tt>. When the handle
    is no longer needed, the appropriate close function must be called. However, if a function is 
    aborted by an exception, this is difficult to ensure. Class HDF5HandleShared is a smart pointer 
    that solves this problem by calling the close function in the destructor of the handle's last 
    owner (This is analogous to how <tt>std::shared_ptr</tt> calls 'delete' on the contained
    pointer). A pointer to the close function must be passed to the constructor, along with an error 
    message that is raised when creation/opening fails. 
    
    When a <tt>HDF5HandleShared</tt> is created or assigned from another one, ownership is shared
    between the two handles, and the value returned by <tt>use_count()</tt> increases by one.
    
    Since <tt>HDF5HandleShared</tt> objects are convertible to <tt>hid_t</tt>, they can be used in the code 
    in place of the latter.

    <b>Usage:</b>

    \code
    HDF5HandleShared file_id(H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT), 
                             &H5Fclose, 
                             "Error message when H5Fopen() fails.");
                       
    ... // use file_id in the same way as a plain hid_t object
    
    // share ownership between same_id and file_id
    HDF5HandleShared same_id(file_id);
    assert(same_id.use_count() == 2);
    assert(same_id.get() == file_id.get());
    \endcode

    <b>\#include</b> \<vigra/hdf5impex.hxx\><br>
    Namespace: vigra
    */
class HDF5HandleShared
{
public:
    typedef herr_t (*Destructor)(hid_t);
    
private:
    hid_t handle_;
    Destructor destructor_;
    size_t * refcount_;
    
public:

        /** \brief Default constructor.
            Creates a NULL handle.
        **/
    HDF5HandleShared()
    : handle_( 0 ),
      destructor_(0),
      refcount_(0)
    {}

        /** \brief Create a shared wrapper for a plain hid_t object.

        The hid_t object \a h is assumed to be the return value of an open or create function.
        It will be closed with the given close function \a destructor as soon as this 
        HDF5HandleShared is destructed, except when \a destructor is a NULL pointer (in which
        case nothing happens at destruction time). If \a h has a value that indicates
        failed opening or creation (by HDF5 convention, this means that \a h is negative),
        an exception is raised by calling <tt>vigra_fail(error_message)</tt>.

        <b>Usage:</b>

        \code
        HDF5HandleShared file_id(H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT), 
                                 &H5Fclose, 
                                 "Error message.");
                           
        ... // use file_id in the same way
        \endcode
        */
    HDF5HandleShared(hid_t h, Destructor destructor, const char * error_message)
    : handle_( h ),
      destructor_(destructor),
      refcount_(0)
    {
        if(handle_ < 0)
            vigra_fail(error_message);
        if(handle_ > 0)
            refcount_ = new size_t(1);
    }

        /** \brief Copy constructor.
            Shares ownership with the RHS handle (analogous to <tt>std::shared_ptr</tt>).
        */
    HDF5HandleShared(HDF5HandleShared const & h)
    : handle_( h.handle_ ),
      destructor_(h.destructor_),
      refcount_(h.refcount_)
    {
        if(refcount_)
            ++(*refcount_);
    }
    
        /** \brief Assignment.
            Call close() for the present LHS handle and share ownership with the 
            RHS handle (analogous to <tt>std::shared_ptr</tt>).
        */
    HDF5HandleShared & operator=(HDF5HandleShared const & h)
    {
        if(h.handle_ != handle_)
        {
            close();
            handle_ = h.handle_;
            destructor_ = h.destructor_;
            refcount_ = h.refcount_;
            if(refcount_)
                ++(*refcount_);
        }
        return *this;
    }

        /** \brief Destructor (calls close())
        */
    ~HDF5HandleShared()
    {
        close();
    }
    
        /** \brief Close the handle if this is the unique (i.e. last) owner. 
        
             Decrements the reference counter and calls the destructor function of 
             the handle (if one has been registered in the constructor) when the counter
             reaches zero. Sets this wrapper to NULL in any case. Returns
             a negative value when the destructor call for the handle fails, and
             a non-negative value otherwise.
        */
    herr_t close()
    {
        herr_t res = 1;
        if(refcount_)
        {
            --(*refcount_);
            if(*refcount_ == 0) 
            {
                if(destructor_)
                    res = (*destructor_)(handle_);
                delete refcount_;
            }
        }
        handle_ = 0;
        destructor_ = 0;
        refcount_ = 0;
        return res;
    }
    
        /** \brief Reset the handle to a new value. 
        
             Equivalent to <tt>handle = HDF5HandleShared(h, destructor, error_message)</tt>.
        */
    void reset(hid_t h, Destructor destructor, const char * error_message)
    {
        if(h < 0)
            vigra_fail(error_message);
        if(h != handle_)
        {
            close();
            handle_ = h;
            destructor_ = destructor;
            if(handle_ > 0)
                refcount_ = new size_t(1);
        }
    }
    
        /** \brief Get the number of owners of the contained handle.
        */
    size_t use_count() const
    {
        return refcount_
                 ? *refcount_
                 : 0;
    }
    
        /** \brief Check if this is the unique owner of the conatined handle.
        
            Equivalent to <tt>handle.use_count() == 1</tt>.
        */
    bool unique() const
    {
        return use_count() == 1;
    }
    
        /** \brief Swap the contents of two handle wrappers.
        
            Also available as <tt>std::swap(handle1, handle2)</tt>.
        */    
    void swap(HDF5HandleShared & h)
    {
        std::swap(handle_, h.handle_);
        std::swap(destructor_, h.destructor_);
        std::swap(refcount_, h.refcount_);
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
        HDF5HandleShared objects in user code. Do not call a close function on the return 
        value - a crash will be likely otherwise.
        */
    operator hid_t() const
    {
        return handle_;
    }

        /** \brief Equality comparison of the contained handle.
        */
    bool operator==(HDF5HandleShared const & h) const
    {
        return handle_ == h.handle_;
    }

        /** \brief Equality comparison of the contained handle.
        */
    bool operator==(hid_t h) const
    {
        return handle_ == h;
    }

        /** \brief Inequality comparison of the contained handle.
        */
    bool operator!=(HDF5HandleShared const & h) const
    {
        return handle_ != h.handle_;
    }

        /** \brief Inequality comparison of the contained handle.
        */
    bool operator!=(hid_t h) const
    {
        return handle_ != h;
    }
};

} // namespace vigra

namespace std {

inline void swap(vigra::HDF5Handle & l, vigra::HDF5Handle & r)
{
    l.swap(r);
}

inline void swap(vigra::HDF5HandleShared & l, vigra::HDF5HandleShared & r)
{
    l.swap(r);
}

} // namespace std

namespace vigra {

/********************************************************/
/*                                                      */
/*                   HDF5ImportInfo                     */
/*                                                      */
/********************************************************/

/** \brief Argument object for the function readHDF5().

See \ref readHDF5() for a usage example. This object must be
used to read an image or array from an HDF5 file 
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
            
            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This function therefore reverses the axis
            order relative to the file contents. That is, when the axes in the file are 
            ordered as 'z', 'y', 'x', this function will return the shape in the order
            'x', 'y', 'z'.
         */
    VIGRA_EXPORT ArrayVector<hsize_t> const & shape() const
    {
        return m_dims;
    }

        /** Get the shape (length) of the dataset along dimension \a dim.
            
            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This function therefore reverses the axis
            order relative to the file contents. That is, when the axes in the file are 
            ordered as 'z', 'y', 'x', this function will return the shape in the order
            'x', 'y', 'z'.
         */
    VIGRA_EXPORT MultiArrayIndex shapeOfDimension(const int dim) const;

        /** Query the pixel type of the dataset.

            Possible values are:
            <DL>
            <DT>"INT8"<DD> 8-bit signed integer (unsigned char)
            <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"INT16"<DD> 16-bit signed integer (short)
            <DT>"UINT16"<DD> 16-bit unsigned integer (unsigned short)
            <DT>"INT32"<DD> 32-bit signed integer (long)
            <DT>"UINT32"<DD> 32-bit unsigned integer (unsigned long)
            <DT>"INT64"<DD> 64-bit signed integer (long long)
            <DT>"UINT64"<DD> 64-bit unsigned integer (unsigned long long)
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
    HDF5HandleShared m_file_handle, m_dataset_handle;
    std::string m_filename, m_path, m_pixeltype;
    hssize_t m_dimensions;
    ArrayVector<hsize_t> m_dims;
};


namespace detail {

template <class T>
struct HDF5TypeTraits
{
    static hid_t getH5DataType()
    {
        std::runtime_error("getH5DataType(): invalid type");
        return 0;
    }
    
    static int numberOfBands()
    {
        std::runtime_error("numberOfBands(): invalid type");
        return 0;
    }
};

template<class T>
inline hid_t getH5DataType()
{
    return HDF5TypeTraits<T>::getH5DataType();
}

#define VIGRA_H5_DATATYPE(type, h5type) \
template <> \
struct HDF5TypeTraits<type> \
{ \
    static hid_t getH5DataType() \
    { \
        return h5type; \
    } \
    static int numberOfBands() \
    { \
        return 1; \
    } \
    typedef type value_type; \
}; \
template <int M> \
struct HDF5TypeTraits<TinyVector<type, M> > \
{ \
    static hid_t getH5DataType() \
    { \
        return h5type; \
    } \
    static int numberOfBands() \
    { \
        return M; \
    } \
    typedef type value_type; \
}; \
template <> \
struct HDF5TypeTraits<RGBValue<type> > \
{ \
    static hid_t getH5DataType() \
    { \
        return h5type; \
    } \
    static int numberOfBands() \
    { \
        return 3; \
    } \
    typedef type value_type; \
}; \

VIGRA_H5_DATATYPE(char, H5T_NATIVE_CHAR)
VIGRA_H5_DATATYPE(signed char, H5T_NATIVE_SCHAR)
VIGRA_H5_DATATYPE(unsigned char, H5T_NATIVE_UCHAR)
VIGRA_H5_DATATYPE(signed short, H5T_NATIVE_SHORT)
VIGRA_H5_DATATYPE(unsigned short, H5T_NATIVE_USHORT)
VIGRA_H5_DATATYPE(signed int, H5T_NATIVE_INT)
VIGRA_H5_DATATYPE(unsigned int, H5T_NATIVE_UINT)
VIGRA_H5_DATATYPE(signed long, H5T_NATIVE_LONG)
VIGRA_H5_DATATYPE(unsigned long, H5T_NATIVE_ULONG)
VIGRA_H5_DATATYPE(signed long long, H5T_NATIVE_LLONG)
VIGRA_H5_DATATYPE(unsigned long long, H5T_NATIVE_ULLONG)
VIGRA_H5_DATATYPE(float, H5T_NATIVE_FLOAT)
VIGRA_H5_DATATYPE(double, H5T_NATIVE_DOUBLE)
VIGRA_H5_DATATYPE(long double, H5T_NATIVE_LDOUBLE)

// char arrays with flexible length require 'handcrafted' H5 datatype
template <>
struct HDF5TypeTraits<char*>
{
    static hid_t getH5DataType()
    {
        hid_t stringtype = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringtype, H5T_VARIABLE);
        return stringtype;
    }
    
    static int numberOfBands()
    {
        return 1;
    }
};

template <>
struct HDF5TypeTraits<const char*>
{
    static hid_t getH5DataType()
    {
        hid_t stringtype = H5Tcopy (H5T_C_S1);
        H5Tset_size(stringtype, H5T_VARIABLE);
        return stringtype;
    }
    
    static int numberOfBands()
    {
        return 1;
    }
};

#undef VIGRA_H5_DATATYPE

#if 0
template<>
inline hid_t getH5DataType<FFTWComplex<float> >()
{
    hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof (FFTWComplex<float>));
    H5Tinsert (complex_id, "real", 0, H5T_NATIVE_FLOAT);
    H5Tinsert (complex_id, "imaginary", sizeof(float), H5T_NATIVE_FLOAT);
    return complex_id;
}

template<>
inline hid_t getH5DataType<FFTWComplex<double> >()
{
    hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof (FFTWComplex<double>));
    H5Tinsert (complex_id, "real", 0, H5T_NATIVE_DOUBLE);
    H5Tinsert (complex_id, "imaginary", sizeof(double), H5T_NATIVE_DOUBLE);
    return complex_id;
}
#endif


} // namespace detail

// helper friend function for callback HDF5_ls_inserter_callback()
void HDF5_ls_insert(void*, const std::string &);
// callback function for ls(), called via HDF5File::H5Literate()
// see http://www.parashift.com/c++-faq-lite/pointers-to-members.html#faq-33.2
// for as to why.

VIGRA_EXPORT H5O_type_t HDF5_get_type(hid_t, const char*);
extern "C" VIGRA_EXPORT herr_t HDF5_ls_inserter_callback(hid_t, const char*, const H5L_info_t*, void*);

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

Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
Fortran-order, while HDF5 uses C-order. This means that a VIGRA MultiArray,
whose indices represent the 'x'-, 'y'-, and 'z'-axis in that order, is reversed
upon writing to an HDF5 file, i.e. in the file the axis order is 'z', 'y', 'x'. 
Likewise, the order is reversed upon reading.

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
  protected:
    HDF5HandleShared fileHandle_;

    // current group handle
    HDF5Handle cGroupHandle_;
    
  private:
    // time tagging of datasets, turned off (= 0) by default.
    int track_time;
    
    bool read_only_;

    // helper class for ls()
    struct ls_closure
    {
        virtual void insert(const std::string &) = 0;
        virtual ~ls_closure() {}
    };
    // datastructure to hold a list of dataset and group names
    struct lsOpData : public ls_closure
    {
        std::vector<std::string> & objects;
        lsOpData(std::vector<std::string> & o) : objects(o) {}
        void insert(const std::string & x)
        {
            objects.push_back(x);
        }
    };
    // (associative-)container closure
    template<class Container>
    struct ls_container_data : public ls_closure
    {
        Container & objects;
        ls_container_data(Container & o) : objects(o) {}
        void insert(const std::string & x)
        {
            objects.insert(std::string(x));
        }
    };

  public:

        // helper for callback HDF5_ls_inserter_callback(), used by ls()
    friend void HDF5_ls_insert(void*, const std::string &);

        /** \brief Set how a file is opened.

            OpenMode::New creates a new file. If the file already exists, overwrite it.

            OpenMode::Open opens a file for reading/writing. The file will be created,
                           if necessary.
        */
    enum OpenMode {
        New,           // Create new empty file (existing file will be deleted).
        Open,          // Open file. Create if not existing.
        OpenReadOnly   // Open file in read-only mode.
    };

        /** \brief Default constructor.

        A file can later be opened via the open() function.
        
        If \a track_creation_times is non-zero, time tagging of datasets will be enabled (it is disabled
        by default).
        */
    HDF5File(int track_creation_times = 0)
    : track_time(track_creation_times),
      read_only_(true)
    {}

        /** \brief Open or create an HDF5File object.

        Creates or opens HDF5 file with given filename. 
        The current group is set to "/".
        
        Note that the HDF5File class is not copyable (the copy constructor is 
        private to enforce this).
        */
    HDF5File(std::string filename, OpenMode mode, int track_creation_times = 0)
    : track_time(track_creation_times)
    {
        open(filename, mode);
    }

        /** \brief Copy a HDF5File object.

            The new object will refer to the same file and group as \a other.
        */
    HDF5File(HDF5File const & other)
    : fileHandle_(other.fileHandle_),
      track_time(other.track_time),
      read_only_(other.read_only_)
    {
        cGroupHandle_ = HDF5Handle(openCreateGroup_(other.currentGroupName_()), &H5Gclose, 
                                   "HDF5File(HDF5File const &): Failed to open group.");
    }

        /** \brief The destructor flushes and closes the file.
         */
    ~HDF5File()
    {
        // The members fileHandle_ and cGroupHandle_ are automatically closed
        // as they are of type HDF5Handle and are properly initialised.
        // The closing of fileHandle_ implies flushing the file to
        // the operating system, see
        // http://www.hdfgroup.org/HDF5/doc/RM/RM_H5F.html#File-Close .
    }
    
        /** \brief Assign a HDF5File object.

            Calls close() on the present file and The new object will refer to the same file and group as \a other.
        */
    HDF5File & operator=(HDF5File const & other)
    {
        if(this != &other)
        {
            close();
            fileHandle_ = other.fileHandle_;
            cGroupHandle_ = HDF5Handle(openCreateGroup_(other.currentGroupName_()), &H5Gclose, 
                                       "HDF5File::operator=(): Failed to open group.");
            track_time = other.track_time;
            read_only_ = other.read_only_;
        }
        return *this;
    }

    int file_use_count() const
    {
        return fileHandle_.use_count();
    }
    
    bool isReadOnly() const
    {
        return read_only_;
    }
    
    void setReadOnly(bool stat=true)
    {
        read_only_ = stat;
    }
  
        /** \brief Open or create the given file in the given mode and set the group to "/".
            If another file is currently open, it is first closed.
         */
    void open(std::string filename, OpenMode mode)
    {
        close();
        
        std::string errorMessage = "HDF5File.open(): Could not open or create file '" + filename + "'.";
        fileHandle_ = HDF5HandleShared(createFile_(filename, mode), &H5Fclose, errorMessage.c_str());
        cGroupHandle_ = HDF5Handle(openCreateGroup_("/"), &H5Gclose, "HDF5File.open(): Failed to open root group.");
        setReadOnly(mode == OpenReadOnly);
    }

        /** \brief Close the current file.
        
            Calls close() on the present file and then assigns itself to the same file and group as \a other.
        */
    void close()
    {
        bool success = cGroupHandle_.close() >= 0 && fileHandle_.close() >= 0;
        vigra_postcondition(success, "HDF5File.close() failed.");
    }

        /** \brief Change current group to "/".
         */
    inline void root()
    {
        std::string message = "HDF5File::root(): Could not open group '/'.";
        cGroupHandle_ = HDF5Handle(H5Gopen(fileHandle_, "/", H5P_DEFAULT),&H5Gclose,message.c_str());
    }

        /** \brief Change the current group.
            Both absolute and relative group names are allowed.
         */
    inline void cd(std::string groupName)
    {
        cGroupHandle_ = getGroupHandle(groupName, "HDF5File::cd()");
    }

        /** \brief Change the current group to its parent group.
            Returns true if successful, false otherwise. If unsuccessful,
            the group will not change.
         */
    inline bool cd_up()
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
    
        /** \brief Change the current group to its parent group.
            Returns true if successful, false otherwise. If unsuccessful,
            the group will not change.
         */
    inline bool cd_up(int levels)
    {
        std::string groupName = currentGroupName_();
        
        for(int i = 0; i<levels; i++)
        {
            if(!cd_up())
            {
                // restore old group if neccessary
                if(groupName != currentGroupName_())
                    cd(groupName);
                return false;
            }
        }
        return true;
    }

        /** \brief Create a new group.
             If the first character is a "/", the path will be interpreted as absolute path,
             otherwise it will be interpreted as path relative to the current group.
        */
    inline void mkdir(std::string groupName)
    {
        vigra_precondition(!isReadOnly(),
            "HDF5File::mkdir(): file is read-only.");
        
        std::string message = "HDF5File::mkdir(): Could not create group '" + groupName + "'.\n";

        // make groupName clean
        groupName = get_absolute_path(groupName);
        
        HDF5Handle(openCreateGroup_(groupName.c_str()),&H5Gclose,message.c_str());
    }

        /** \brief Change the current group; create it if necessary.
             If the first character is a "/", the path will be interpreted as absolute path,
             otherwise it will be interpreted as path relative to the current group.
        */
    inline void cd_mk(std::string groupName)
    {
        vigra_precondition(!isReadOnly(),
            "HDF5File::cd_mk(): file is read-only.");
        
        std::string  message = "HDF5File::cd_mk(): Could not create group '" + groupName + "'.";

        // make groupName clean
        groupName = get_absolute_path(groupName);

        cGroupHandle_ = HDF5Handle(openCreateGroup_(groupName.c_str()),&H5Gclose,message.c_str());
    }

        // helper function for the various ls() variants.
    void ls_H5Literate(ls_closure & data) const
    {
        H5Literate(cGroupHandle_, H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
                   HDF5_ls_inserter_callback, static_cast<void*>(&data));
    }

        /** \brief List the contents of the current group.
            The function returns a vector of strings holding the entries of the
            current group. Only datasets and groups are listed, other objects
            (e.g. datatypes) are ignored. Group names always have a trailing "/".
        */
    inline std::vector<std::string> ls() const
    {
        std::vector<std::string> list;
        lsOpData data(list);
        ls_H5Literate(data);
        return list;
    }

        /** \brief List the contents of the current group into a container-like
                   object via insert(). 
                   
            Only datasets and groups are inserted, other objects (e.g., datatypes) are ignored. 
            Group names always have a trailing "/".

            The argument cont is presumably an associative container, however,
            only its member function <tt>cont.insert(std::string)</tt> will be
            called.
            \param cont      reference to a container supplying a member function
                             <tt>insert(const i_type &)</tt>, where <tt>i_type</tt>
                             is convertible to <tt>std::string</tt>.
        */
    template<class Container>
    void ls(Container & cont) const
    {
        ls_container_data<Container> data(cont);
        ls_H5Literate(data);
    }

        /** \brief Get the path of the current group.
        */
    inline std::string pwd() const
    {
        return currentGroupName_();
    }

        /** \brief Get the name of the associated file.
        */
    inline std::string filename() const
    {
        return fileName_();
    }

        /** \brief Check if given datasetName exists.
        */
    inline bool existsDataset(std::string datasetName) const
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);
        return (H5Lexists(fileHandle_, datasetName.c_str(), H5P_DEFAULT) > 0);
    }

        /** \brief Get the number of dimensions of a certain dataset
             If the first character is a "/", the path will be interpreted as absolute path,
             otherwise it will be interpreted as path relative to the current group.
        */
    hssize_t getDatasetDimensions(std::string datasetName) const
    {
        HDF5Handle datasetHandle = getDatasetHandle(datasetName);

        return getDatasetDimensions_(datasetHandle);
    }

    hssize_t getDatasetDimensions_(hid_t dataset) const
    {
        std::string errorMessage = "HDF5File::getDatasetDimensions(): Unable to access dataspace.";
        HDF5Handle dataspaceHandle(H5Dget_space(dataset), &H5Sclose, errorMessage.c_str());

        //return dimension information
        return H5Sget_simple_extent_ndims(dataspaceHandle);
    }

        /** \brief Get the shape of each dimension of a certain dataset.
            
           Normally, this function is called after determining the dimension of the
            dataset using \ref getDatasetDimensions().
            If the first character is a "/", the path will be interpreted as absolute path,
            otherwise it will be interpreted as path relative to the current group.
            
            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This function therefore reverses the axis
            order relative to the file contents. That is, when the axes in the file are 
            ordered as 'z', 'y', 'x', this function will return the shape in the order
            'x', 'y', 'z'.
        */
    ArrayVector<hsize_t> getDatasetShape(std::string datasetName) const
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

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

        // invert the dimensions to guarantee VIGRA-compatible order.
        std::reverse(shape.begin(), shape.end());
        return shape;
    }

        /** Query the pixel type of the dataset.

            Possible values are:
            <DL>
            <DT>"INT8"<DD> 8-bit signed integer (unsigned char)
            <DT>"UINT8"<DD> 8-bit unsigned integer (unsigned char)
            <DT>"INT16"<DD> 16-bit signed integer (short)
            <DT>"UINT16"<DD> 16-bit unsigned integer (unsigned short)
            <DT>"INT32"<DD> 32-bit signed integer (long)
            <DT>"UINT32"<DD> 32-bit unsigned integer (unsigned long)
            <DT>"INT64"<DD> 64-bit signed integer (long long)
            <DT>"UINT64"<DD> 64-bit unsigned integer (unsigned long long)
            <DT>"FLOAT"<DD> 32-bit floating point (float)
            <DT>"DOUBLE"<DD> 64-bit floating point (double)
            <DT>"UNKNOWN"<DD> any other type
            </DL>
         */
    std::string getDatasetType(std::string const & datasetName) const
    {
        HDF5Handle datasetHandle = getDatasetHandle(datasetName);

        hid_t datatype = H5Dget_type(datasetHandle);
        H5T_class_t dataclass = H5Tget_class(datatype);
        size_t datasize  = H5Tget_size(datatype);
        H5T_sign_t datasign  = H5Tget_sign(datatype);

        if(dataclass == H5T_FLOAT)
        {
            if(datasize == 4)
                return "FLOAT";
            else if(datasize == 8)
                return "DOUBLE";
        }
        else if(dataclass == H5T_INTEGER)   
        {
            if(datasign == H5T_SGN_NONE)
            {
                if(datasize ==  1)
                    return "UINT8";
                else if(datasize == 2)
                    return "UINT16";
                else if(datasize == 4)
                    return "UINT32";
                else if(datasize == 8)
                    return "UINT64";
            }
            else
            {
                if(datasize ==  1)
                    return "INT8";
                else if(datasize == 2)
                    return "INT16";
                else if(datasize == 4)
                    return "INT32";
                else if(datasize == 8)
                    return "INT64";
            }
        }
        return "UNKNOWN";
    }
        
        /** \brief Obtain the HDF5 handle of a dataset.
        */
    HDF5Handle getDatasetHandle(std::string const & datasetName) const
    {
        std::string errorMessage = "HDF5File::getDatasetHandle(): Unable to open dataset '" + datasetName + "'.";
        return HDF5Handle(getDatasetHandle_(get_absolute_path(datasetName)), &H5Dclose, errorMessage.c_str());
    }
        
        /** \brief Obtain a shared HDF5 handle of a dataset.
        */
    HDF5HandleShared getDatasetHandleShared(std::string const & datasetName) const
    {
        std::string errorMessage = "HDF5File::getDatasetHandle(): Unable to open dataset '" + datasetName + "'.";
        return HDF5HandleShared(getDatasetHandle_(get_absolute_path(datasetName)), &H5Dclose, errorMessage.c_str());
    }

        /** \brief Obtain the HDF5 handle of a group (create the group if it doesn't exist).
         */
    HDF5Handle getGroupHandle(std::string group_name, 
                              std::string function_name = "HDF5File::getGroupHandle()")
    {
        std::string errorMessage = function_name + ": Group '" + group_name + "' not found.";

        // make group_name clean
        group_name = get_absolute_path(group_name);

        // group must exist
        vigra_precondition(group_name == "/" || H5Lexists(fileHandle_, group_name.c_str(), H5P_DEFAULT) != 0, 
                           errorMessage.c_str());

        // open group and return group handle
        return HDF5Handle(openCreateGroup_(group_name), &H5Gclose, "Internal error");
    }

        /** \brief Obtain the HDF5 handle of a attribute.
         */
    HDF5Handle getAttributeHandle(std::string dataset_name, std::string attribute_name) const
    {
        std::string message = "HDF5File::getAttributeHandle(): Attribute '" + attribute_name + "' not found.";
        return HDF5Handle(H5Aopen(getDatasetHandle(dataset_name), attribute_name.c_str(), H5P_DEFAULT),
                          &H5Aclose, message.c_str());
    }

    /* Writing Attributes */

        /** \brief Write MultiArray Attributes.
          * In contrast to datasets, subarray access, chunks and compression are not available.
          */
    template<unsigned int N, class T, class Stride>
    inline void writeAttribute(std::string object_name, 
                               std::string attribute_name, 
                               const MultiArrayView<N, T, Stride> & array)
    {
        // make object_name clean
        object_name = get_absolute_path(object_name);

        write_attribute_(object_name, attribute_name, array, detail::getH5DataType<T>(), 1);
    }

    template<unsigned int N, class T, int SIZE, class Stride>
    inline void writeAttribute(std::string datasetName, 
                               std::string attributeName, 
                               const MultiArrayView<N, TinyVector<T, SIZE>, Stride> & array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        write_attribute_(datasetName, attributeName, array, detail::getH5DataType<T>(), SIZE);
    }

    template<unsigned int N, class T, class Stride>
    inline void writeAttribute(std::string datasetName, 
                               std::string attributeName, 
                               const MultiArrayView<N, RGBValue<T>, Stride> & array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        write_attribute_(datasetName, attributeName, array, detail::getH5DataType<T>(), 3);
    }

        /** \brief Write a single value.
          Specialization of the write function for simple datatypes
         */
    inline void writeAttribute(std::string object_name, std::string attribute_name, char data) 
        { writeAtomicAttribute(object_name,attribute_name,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, signed char data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, signed short data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, signed int data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, signed long data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, signed long long data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, unsigned char data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, unsigned short data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, unsigned int data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, unsigned long data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, unsigned long long data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, float data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, double data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, long double data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, const char* data) 
        { writeAtomicAttribute(datasetName,attributeName,data); }
    inline void writeAttribute(std::string datasetName, std::string attributeName, std::string const & data) 
        { writeAtomicAttribute(datasetName,attributeName,data.c_str()); }

        /** \brief Test if attribute exists.
        */
    bool existsAttribute(std::string object_name, std::string attribute_name)
    {
        std::string obj_path = get_absolute_path(object_name);
        htri_t exists = H5Aexists_by_name(fileHandle_, obj_path.c_str(),
                                          attribute_name.c_str(), H5P_DEFAULT);
        vigra_precondition(exists >= 0, "HDF5File::existsAttribute(): "
                                        "object '" + object_name + "' "
                                        "not found.");
        return exists != 0;
    }

    // Reading Attributes

        /** \brief Read MultiArray Attributes.
          * In contrast to datasets, subarray access is not available.
          */
    template<unsigned int N, class T, class Stride>
    inline void readAttribute(std::string object_name, 
                              std::string attribute_name, 
                              MultiArrayView<N, T, Stride> array)
    {
        // make object_name clean
        object_name = get_absolute_path(object_name);

        read_attribute_(object_name, attribute_name, array, detail::getH5DataType<T>(), 1);
    }

    template<unsigned int N, class T, int SIZE, class Stride>
    inline void readAttribute(std::string datasetName, 
                              std::string attributeName, 
                              MultiArrayView<N, TinyVector<T, SIZE>, Stride> array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        read_attribute_(datasetName, attributeName, array, detail::getH5DataType<T>(), SIZE);
    }

    template<unsigned int N, class T, class Stride>
    inline void readAttribute(std::string datasetName, 
                              std::string attributeName, 
                              MultiArrayView<N, RGBValue<T>, Stride> array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        read_attribute_(datasetName, attributeName, array, detail::getH5DataType<T>(), 3);
    }

        /** \brief Read a single value.
          Specialization of the read function for simple datatypes
         */
    inline void readAttribute(std::string object_name, std::string attribute_name, char &data)       
        { readAtomicAttribute(object_name,attribute_name,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, signed char &data)        
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, signed short &data)       
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, signed int &data)       
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, signed long &data)       
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, signed long long &data)       
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, unsigned char &data)       
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, unsigned short &data)      
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, unsigned int &data)      
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, unsigned long &data)      
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, unsigned long long &data)      
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, float &data)       
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, double &data)      
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, long double &data) 
        { readAtomicAttribute(datasetName,attributeName,data); }
    inline void readAttribute(std::string datasetName, std::string attributeName, std::string &data) 
        { readAtomicAttribute(datasetName,attributeName,data); }

    // Writing data

        /** \brief Write multi arrays.
          
            Chunks can be activated by setting 
            \code iChunkSize = size; //size \> 0 
            \endcode .
            The chunks will be hypercubes with edge length size. When <tt>iChunkSize == 0</tt>
            (default), the behavior depends on the <tt>compression</tt> setting: If no
            compression is requested, the data is written without chunking. Otherwise,
            chuning is required, and the chunk size is automatically selected such that
            each chunk contains about 300k pixels.

            Compression can be activated by setting 
            \code compression = parameter; // 0 \< parameter \<= 9 
            \endcode
            where 0 stands for no compression and 9 for maximum compression.

            If the first character of datasetName is a "/", the path will be interpreted as absolute path,
            otherwise it will be interpreted as path relative to the current group.

            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This means that a VIGRA MultiArray,
            whose indices represent the 'x'-, 'y'-, and 'z'-axis in that order, is reversed
            upon writing to an HDF5 file, i.e. in the file the axis order is 'z', 'y', 'x'. 
        */
    template<unsigned int N, class T, class Stride>
    inline void write(std::string datasetName, 
                      const MultiArrayView<N, T, Stride> & array, 
                      int iChunkSize = 0, int compression = 0)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        typename MultiArrayShape<N>::type chunkSize;
        for(unsigned int i = 0; i < N; i++){
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

            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This means that a VIGRA MultiArray,
            whose indices represent the 'x'-, 'y'-, and 'z'-axis in that order, is reversed
            upon writing to an HDF5 file, i.e. in the file the axis order is 'z', 'y', 'x'. 
        */
    template<unsigned int N, class T, class Stride>
    inline void write(std::string datasetName, 
                      const MultiArrayView<N, T, Stride> & array, 
                      typename MultiArrayShape<N>::type chunkSize, int compression = 0)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        write_(datasetName, array, detail::getH5DataType<T>(), 1, chunkSize, compression);
    }

        /** \brief Write a multi array into a larger volume.
            blockOffset determines the position, where array is written.

            If the first character of datasetName is a "/", the path will be interpreted as absolute path,
            otherwise it will be interpreted as path relative to the current group.

            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This means that a VIGRA MultiArray,
            whose indices represent the 'x'-, 'y'-, and 'z'-axis in that order, is reversed
            upon writing to an HDF5 file, i.e. in the file the axis order is 'z', 'y', 'x'. 
        */
    template<unsigned int N, class T, class Stride>
    inline void writeBlock(std::string datasetName, 
                           typename MultiArrayShape<N>::type blockOffset, 
                           const MultiArrayView<N, T, Stride> & array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);
        typedef detail::HDF5TypeTraits<T> TypeTraits;
        writeBlock_(datasetName, blockOffset, array, 
                    TypeTraits::getH5DataType(), TypeTraits::numberOfBands());
    }

    template<unsigned int N, class T, class Stride>
    inline herr_t writeBlock(HDF5HandleShared dataset, 
                             typename MultiArrayShape<N>::type blockOffset, 
                             const MultiArrayView<N, T, Stride> & array)
    {
        typedef detail::HDF5TypeTraits<T> TypeTraits;
        return writeBlock_(dataset, blockOffset, array, 
                           TypeTraits::getH5DataType(), TypeTraits::numberOfBands());
    }

    // non-scalar (TinyVector) and unstrided multi arrays
    template<unsigned int N, class T, int SIZE, class Stride>
    inline void write(std::string datasetName, 
                      const MultiArrayView<N, TinyVector<T, SIZE>, Stride> & array, 
                      int iChunkSize = 0, int compression = 0)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        typename MultiArrayShape<N>::type chunkSize;
        for(int i = 0; i < N; i++){
            chunkSize[i] = iChunkSize;
        }
        write_(datasetName, array, detail::getH5DataType<T>(), SIZE, chunkSize, compression);
    }

    template<unsigned int N, class T, int SIZE, class Stride>
    inline void write(std::string datasetName, 
                      const MultiArrayView<N, TinyVector<T, SIZE>, Stride> & array, 
                      typename MultiArrayShape<N>::type chunkSize, int compression = 0)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        write_(datasetName, array, detail::getH5DataType<T>(), SIZE, chunkSize, compression);
    }

        /** \brief Write array vectors.
          
            Compression can be activated by setting 
            \code compression = parameter; // 0 \< parameter \<= 9 
            \endcode
            where 0 stands for no compression and 9 for maximum compression.

            If the first character of datasetName is a "/", the path will be interpreted as absolute path,
            otherwise it will be interpreted as path relative to the current group.
        */
    template<class T>
    void write(const std::string & datasetName,
                      const ArrayVectorView<T> & array,
                      int compression = 0)
    {
        // convert to a (trivial) MultiArrayView and forward.
        MultiArrayShape<1>::type shape(array.size());
        const MultiArrayView<1, T> m_array(shape, const_cast<T*>(array.data()));
        write(datasetName, m_array, compression);
    }

    // non-scalar (RGBValue) and unstrided multi arrays
    template<unsigned int N, class T, class Stride>
    inline void write(std::string datasetName, 
                      const MultiArrayView<N, RGBValue<T>, Stride> & array, 
                      int iChunkSize = 0, int compression = 0)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        typename MultiArrayShape<N>::type chunkSize;
        for(int i = 0; i < N; i++){
            chunkSize[i] = iChunkSize;
        }
        write_(datasetName, array, detail::getH5DataType<T>(), 3, chunkSize, compression);
    }

    template<unsigned int N, class T, class Stride>
    inline void write(std::string datasetName, 
                      const MultiArrayView<N, RGBValue<T>, Stride> & array, 
                      typename MultiArrayShape<N>::type chunkSize, int compression = 0)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        write_(datasetName, array, detail::getH5DataType<T>(), 3, chunkSize, compression);
    }

         /** \brief Write a single value.
            Specialization of the write function for simple datatypes
         */
    inline void write(std::string datasetName, char data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, signed char data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, signed short data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, signed int data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, signed long data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, signed long long data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, unsigned char data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, unsigned short data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, unsigned int data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, unsigned long data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, unsigned long long data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, float data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, double data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, long double data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, const char* data) { writeAtomic(datasetName,data); }
    inline void write(std::string datasetName, std::string const & data) { writeAtomic(datasetName,data.c_str()); }

    // Reading data
    
        /** \brief Read data into a multi array.
            If the first character of datasetName is a "/", the path will be interpreted as absolute path,
            otherwise it will be interpreted as path relative to the current group.

            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This means that a HDF5 dataset,
            whose indices represent the 'z'-, 'y'-, and 'x'-axis in that order, is reversed
            upon reading into a MultiArrayView, i.e. in the array axis order must be 'x', 'y', 'z'. 
        */
    template<unsigned int N, class T, class Stride>
    inline void read(std::string datasetName, MultiArrayView<N, T, Stride> array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        read_(datasetName, array, detail::getH5DataType<T>(), 1);
    }

        /** \brief Read data into a MultiArray. Resize MultiArray to the correct size.
            If the first character of datasetName is a "/", the path will be interpreted as absolute path,
            otherwise it will be interpreted as path relative to the current group.

            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This means that a HDF5 dataset,
            whose indices represent the 'z'-, 'y'-, and 'x'-axis in that order, is reversed
            upon reading into a MultiArray, i.e. in the array axis order will be 'x', 'y', 'z'. 
        */
    template<unsigned int N, class T, class Alloc>
    inline void readAndResize(std::string datasetName, MultiArray<N, T, Alloc> & array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        // get dataset dimension
        ArrayVector<hsize_t> dimshape = getDatasetShape(datasetName);

        // check if dimensions are correct
        vigra_precondition(N == MultiArrayIndex(dimshape.size()), // the object in the HDF5 file may have one additional dimension which we then interpret as the pixel type bands
            "HDF5File::readAndResize(): Array dimension disagrees with dataset dimension.");

        // reshape target MultiArray
        typename MultiArrayShape<N>::type shape;
        for(int k=0; k < (int)dimshape.size(); ++k)
            shape[k] = (MultiArrayIndex)dimshape[k];
        array.reshape(shape);

        read_(datasetName, array, detail::getH5DataType<T>(), 1);
    }

        /** \brief Read data into an array vector.
          If the first character of datasetName is a "/", the path will be interpreted as absolute path,
          otherwise it will be interpreted as path relative to the current group.
        */
    template<class T>
    inline void read(const std::string & datasetName, ArrayVectorView<T> array)
    {
        // convert to a (trivial) MultiArrayView and forward.
        MultiArrayShape<1>::type shape(array.size());
        MultiArrayView<1, T> m_array(shape, (array.data()));
        read(datasetName, m_array);
    }

        /** \brief Read data into an array vector. Resize the array vector to the correct size.
            If the first character of datasetName is a "/", the path will be interpreted as absolute path,
            otherwise it will be interpreted as path relative to the current group.
        */
    template<class T>
    inline void readAndResize(std::string datasetName,
                              ArrayVector<T> & array)
    {
        // make dataset name clean
        datasetName = get_absolute_path(datasetName);

        // get dataset dimension
        ArrayVector<hsize_t> dimshape = getDatasetShape(datasetName);

        // check if dimensions are correct
        vigra_precondition(1 == MultiArrayIndex(dimshape.size()),
            "HDF5File::readAndResize(): Array dimension disagrees with Dataset dimension must equal one for vigra::ArrayVector.");

        // resize target array vector
        array.resize((typename ArrayVector<T>::size_type)dimshape[0]);
        // convert to a (trivial) MultiArrayView and forward.
        MultiArrayShape<1>::type shape(array.size());
        MultiArrayView<1, T> m_array(shape, (array.data()));

        read_(datasetName, m_array, detail::getH5DataType<T>(), 1);
    }

        /** \brief Read a block of data into a multi array.
            This function allows to read a small block out of a larger volume stored
            in an HDF5 dataset.

            blockOffset determines the position of the block.
            blockSize determines the size in each dimension of the block.

            If the first character of datasetName is a "/", the path will be interpreted as absolute path,
            otherwise it will be interpreted as path relative to the current group.

            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This means that a HDF5 dataset,
            whose indices represent the 'z'-, 'y'-, and 'x'-axis in that order, is reversed
            upon reading into a MultiArray, i.e. in the array axis order will be 'x', 'y', 'z'. 
        */
    template<unsigned int N, class T, class Stride>
    inline void readBlock(std::string datasetName, 
                          typename MultiArrayShape<N>::type blockOffset, 
                          typename MultiArrayShape<N>::type blockShape, 
                          MultiArrayView<N, T, Stride> array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);
        typedef detail::HDF5TypeTraits<T> TypeTraits;
        readBlock_(datasetName, blockOffset, blockShape, array, 
                   TypeTraits::getH5DataType(), TypeTraits::numberOfBands());
    }

    template<unsigned int N, class T, class Stride>
    inline herr_t readBlock(HDF5HandleShared dataset, 
                          typename MultiArrayShape<N>::type blockOffset, 
                          typename MultiArrayShape<N>::type blockShape, 
                          MultiArrayView<N, T, Stride> array)
    {
        typedef detail::HDF5TypeTraits<T> TypeTraits;
        return readBlock_(dataset, blockOffset, blockShape, array, 
                          TypeTraits::getH5DataType(), TypeTraits::numberOfBands());
    }

    // non-scalar (TinyVector) and unstrided target MultiArrayView
    template<unsigned int N, class T, int SIZE, class Stride>
    inline void read(std::string datasetName, MultiArrayView<N, TinyVector<T, SIZE>, Stride> array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        read_(datasetName, array, detail::getH5DataType<T>(), SIZE);
    }

    // non-scalar (TinyVector) MultiArray
    template<unsigned int N, class T, int SIZE, class Alloc>
    inline void readAndResize(std::string datasetName, MultiArray<N, TinyVector<T, SIZE>, Alloc> & array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        // get dataset dimension
        ArrayVector<hsize_t> dimshape = getDatasetShape(datasetName);

        // check if dimensions are correct
        vigra_precondition((N+1) ==  MultiArrayIndex(dimshape.size()) &&
                           SIZE == dimshape[0], // the object in the HDF5 file must have one additional dimension which we interpret as the pixel type bands
            "HDF5File::readAndResize(): Array dimension disagrees with dataset dimension.");
        
        // reshape target MultiArray
        typename MultiArrayShape<N>::type shape;
        for(int k=1; k < (int)dimshape.size(); ++k)
            shape[k-1] = (MultiArrayIndex)dimshape[k];
        array.reshape(shape);

        read_(datasetName, array, detail::getH5DataType<T>(), SIZE);
    }

    // non-scalar (RGBValue) and unstrided target MultiArrayView
    template<unsigned int N, class T, class Stride>
    inline void read(std::string datasetName, MultiArrayView<N, RGBValue<T>, Stride> array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        read_(datasetName, array, detail::getH5DataType<T>(), 3);
    }

    // non-scalar (RGBValue) MultiArray
    template<unsigned int N, class T, class Alloc>
    inline void readAndResize(std::string datasetName, MultiArray<N, RGBValue<T>, Alloc> & array)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        // get dataset dimension
        ArrayVector<hsize_t> dimshape = getDatasetShape(datasetName);

        // check if dimensions are correct
        vigra_precondition((N+1) ==  MultiArrayIndex(dimshape.size()) &&
                           3 == dimshape[0], // the object in the HDF5 file must have one additional dimension which we interpret as the pixel type bands
            "HDF5File::readAndResize(): Array dimension disagrees with dataset dimension.");

        // reshape target MultiArray
        typename MultiArrayShape<N>::type shape;
        for(int k=1; k < (int)dimshape.size(); ++k)
            shape[k-1] = (MultiArrayIndex)dimshape[k];
        array.reshape(shape);

        read_(datasetName, array, detail::getH5DataType<T>(), 3);
    }

        /** \brief Read a single value.
            Specialization of the read function for simple datatypes
         */
    inline void read(std::string datasetName, char &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, signed char &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, signed short &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, signed int &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, signed long &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, signed long long &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, unsigned char &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, unsigned short &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, unsigned int &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, unsigned long &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, unsigned long long &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, float &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, double &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, long double &data) { readAtomic(datasetName,data); }
    inline void read(std::string datasetName, std::string &data) { readAtomic(datasetName,data); }

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

            Note that the memory order between VIGRA and HDF5 files differs: VIGRA uses 
            Fortran-order, while HDF5 uses C-order. This means that a VIGRA MultiArray,
            whose indices represent the 'x'-, 'y'-, and 'z'-axis in that order, is reversed
            upon writing to an HDF5 file, i.e. in the file the axis order is 'z', 'y', 'x'. 
        */
    template<unsigned int N, class T>
    HDF5HandleShared createDataset(std::string datasetName, 
                                   TinyVector<MultiArrayIndex, N> const & shape, 
                                   T init = T(), 
                                   int iChunkSize = 0, 
                                   int compressionParameter = 0)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        typename MultiArrayShape<N>::type chunkSize;
        for(int i = 0; i < N; i++){
            chunkSize[i] = iChunkSize;
        }
        return createDataset<N,T>(datasetName, shape, init, chunkSize, compressionParameter);
    }

        // FIXME: implement this in terms of createDatasetImpl()
    template<unsigned int N, class T>
    HDF5HandleShared createDataset(std::string datasetName, 
                                   TinyVector<MultiArrayIndex, N> const & shape, 
                                   T init, 
                                   TinyVector<MultiArrayIndex, N> const & chunkSize, 
                                   int compressionParameter = 0)
    {
        vigra_precondition(!isReadOnly(),
            "HDF5File::createDataset(): file is read-only.");
        
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

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

        // turn off time tagging of datasets by default.
        H5Pset_obj_track_times(plist, track_time);

        // enable chunks
        ArrayVector<hsize_t> chunks(defineChunks(chunkSize, shape, 1, compressionParameter));
        if(chunks.size() > 0)
        {
            std::reverse(chunks.begin(), chunks.end());
            H5Pset_chunk (plist, chunks.size(), chunks.begin());
        }

        // enable compression
        if(compressionParameter > 0)
        {
            H5Pset_deflate(plist, compressionParameter);
        }

        //create the dataset.
        HDF5HandleShared datasetHandle(H5Dcreate(parent, setname.c_str(), detail::getH5DataType<T>(), 
                                                 dataspaceHandle, H5P_DEFAULT, plist, H5P_DEFAULT),
                                       &H5Dclose, 
                                       "HDF5File::createDataset(): unable to create dataset.");
        if(parent != cGroupHandle_)
            H5Gclose(parent);
            
        return datasetHandle;
    }
    
    template<class T, unsigned int N>
    HDF5HandleShared 
    createDatasetImpl(std::string datasetName, 
                      TinyVector<MultiArrayIndex, N> const & shape, 
                      TinyVector<MultiArrayIndex, N> const & chunkSize,
                      int compressionParameter = 0)
    {
        vigra_precondition(!isReadOnly(),
            "HDF5File::createDataset(): file is read-only.");
        
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        std::string groupname = SplitString(datasetName).first();
        std::string setname = SplitString(datasetName).last();

        hid_t parent = openCreateGroup_(groupname);

        // delete the dataset if it already exists
        deleteDataset_(parent, setname);

        // invert dimensions to guarantee c-order
        // add an extra dimension in case that the data is non-scalar
        typedef detail::HDF5TypeTraits<T> TypeTraits;
        ArrayVector<hsize_t> shape_inv;
        if(TypeTraits::numberOfBands() > 1)
        {
            shape_inv.resize(N+1);
            shape_inv[N] = TypeTraits::numberOfBands();
        }
        else
        {
            shape_inv.resize(N);
        }
        for(unsigned int k=0; k<N; ++k)
            shape_inv[N-1-k] = shape[k];

        // create dataspace
        HDF5Handle 
        dataspaceHandle = HDF5Handle(H5Screate_simple(shape_inv.size(), shape_inv.data(), NULL),
                                    &H5Sclose, "HDF5File::createDataset(): unable to create dataspace for scalar data.");

        // set fill value
        HDF5Handle plist ( H5Pcreate(H5P_DATASET_CREATE), &H5Pclose, "HDF5File::createDataset(): unable to create property list." );
        typename TypeTraits::value_type init = TypeTraits::value_type();
        H5Pset_fill_value(plist, TypeTraits::getH5DataType(), &init);

        // turn off time tagging of datasets by default.
        H5Pset_obj_track_times(plist, track_time);

        // enable chunks
        ArrayVector<hsize_t> chunks(defineChunks(chunkSize, shape, TypeTraits::numberOfBands(), compressionParameter));
        if(chunks.size() > 0)
        {
            std::reverse(chunks.begin(), chunks.end());
            H5Pset_chunk (plist, chunks.size(), chunks.begin());
        }

        // enable compression
        if(compressionParameter > 0)
        {
            H5Pset_deflate(plist, compressionParameter);
        }

        //create the dataset.
        HDF5HandleShared datasetHandle(H5Dcreate(parent, setname.c_str(), 
                                                 TypeTraits::getH5DataType(), 
                                                 dataspaceHandle, H5P_DEFAULT, plist, H5P_DEFAULT),
                                       &H5Dclose, 
                                       "HDF5File::createDataset(): unable to create dataset.");
        if(parent != cGroupHandle_)
            H5Gclose(parent);
            
        return datasetHandle;
    }

        /** \brief Immediately write all data to disk
        */
    inline void flushToDisk()
    {
        H5Fflush(fileHandle_, H5F_SCOPE_GLOBAL);
    }

  private:

        /* Simple extension of std::string for splitting into two parts
         *
         *  Strings (in particular: file/dataset paths) will be split into two
         *  parts. The split is made at the last occurrence of the delimiter.
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
    
    template <class Shape>
    ArrayVector<hsize_t> 
    defineChunks(Shape const & chunks, Shape const & shape, int numBands, int compression = 0)
    {
        if(chunks[0] > 0)
        {
            ArrayVector<hsize_t> res(chunks.begin(), chunks.end());
            if(numBands > 1)
                res.insert(res.begin(), numBands);
            return res;
        }
        else if(compression > 0)
        {
            // set default chunks to enable compression 
            // (arbitrarily include about 300k pixels into each chunk, but make sure
            //  that the chunk size doesn't exceed the shape)
            ArrayVector<hsize_t> res(shape.begin(), shape.end());
            hsize_t chunk_length = (hsize_t)std::pow(300000.0, 1.0 / shape.size());
            for(unsigned int k=0; k < shape.size(); ++k)
                if(res[k] > chunk_length)
                    res[k] = chunk_length;
            if(numBands > 1)
                res.insert(res.begin(), numBands);
            return res;
        }
        else
        {
            return ArrayVector<hsize_t>();
        }
    }

  public:

        /** \brief takes any path and converts it into an absolute path
             in the current file.
           
             Elements like "." and ".." are treated as expected.
             Links are not supported or resolved.
        */
    inline std::string get_absolute_path(std::string path) const {
        // check for empty input or "." and return the current folder
        if(path.length() == 0 || path == "."){
            return currentGroupName_();
        }

        std::string str;
        // convert to absolute path
        if(relativePath_(path)){
            std::string cname = currentGroupName_();
            if (cname == "/")
                str = currentGroupName_()+path;
            else
                str = currentGroupName_()+"/"+path;
        }else{
            str = path;
        }

        // cut out "./"
        std::string::size_type startpos = 0;
        while(str.find(std::string("./"), startpos) != std::string::npos){
            std::string::size_type pos = str.find(std::string("./"), startpos);
            startpos = pos+1;
            // only cut if "./" is not part of "../" (see below)
            if(str.substr(pos-1,3) != "../"){
                // cut out part of the string
                str = str.substr(0,pos) + str.substr(pos+2,str.length()-pos-2);
                startpos = pos;
            }
        }

        // cut out pairs of "bla/../"
        while(str.find(std::string("..")) != std::string::npos){
            std::string::size_type pos = str.find(std::string(".."));

            // find first slash after ".."
            std::string::size_type end = str.find("/",pos);
            if(end != std::string::npos){
                // also include slash
                end++;
            }else{
                // no "/" after ".." --> this is a group, add a "/"
                str = str + "/";
                end = str.length();
            }

            // find first slash before ".."
            std::string::size_type prev_slash = str.rfind("/",pos);
            // if the root slash is the first before ".." --> Error
            vigra_invariant(prev_slash != 0 && prev_slash != std::string::npos,
                            "Error parsing path: "+str);
            // find second slash before ".."
            std::string::size_type begin = str.rfind("/",prev_slash-1);

            // cut out part of the string
            str = str.substr(0,begin+1) + str.substr(end,str.length()-end);
        }

        return str;
    }
    
  protected:

        /* checks if the given path is a relative path.
         */
    inline bool relativePath_(std::string & path) const
    {
        std::string::size_type pos = path.find('/') ;
        if(pos == 0)
            return false;

        return true;
    }

        /* return the name of the current group
         */
    inline std::string currentGroupName_() const
    {
        int len = H5Iget_name(cGroupHandle_,NULL,1000);
        ArrayVector<char> name (len+1,0);
        H5Iget_name(cGroupHandle_,name.begin(),len+1);

        return std::string(name.begin());
    }

        /* return the name of the current file
         */
    inline std::string fileName_() const
    {
        int len = H5Fget_name(fileHandle_,NULL,1000);
        ArrayVector<char> name (len+1,0);
        H5Fget_name(fileHandle_,name.begin(),len+1);

        return std::string(name.begin());
    }

        /* create an empty file and open is
         */
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
        else if(mode == OpenReadOnly) {
            fclose( pFile );
            fileId = H5Fopen(filePath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        }
        else
        {
            fclose(pFile);
            std::remove(filePath.c_str());
            fileId = H5Fcreate(filePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
        return fileId;
    }

        /* \brief Open a group. 
        
           A negative value is returned when the group does not exist or when opening 
           fails for other reasons.
         */
    hid_t openGroup_(std::string groupName) const
    {
        return const_cast<HDF5File *>(this)->openCreateGroup_(groupName, false);
    }

        /* \brief Open or create a group. 
        
           If \a create is <tt>true</tt> and the group does not exist, it will be created, 
           including all necessary parent groups. If group creation fails, a negative
           value is returned. Likewise, a negative value is returned when \a create
           is <tt>false</tt> and the group does not exist or when opening of the group 
           fails for other reasons.
         */
    hid_t openCreateGroup_(std::string groupName, bool create = true)
    {
        // make groupName clean
        groupName = get_absolute_path(groupName);

        // open root group
        hid_t parent = H5Gopen(fileHandle_, "/", H5P_DEFAULT);
        if(groupName == "/")
        {
            return parent;
        }

        // remove leading /
        groupName = std::string(groupName.begin()+1, groupName.end());

        // check if the groupName has finishing slash
        if( groupName.size() != 0 && *groupName.rbegin() != '/')
        {
            groupName = groupName + '/';
        }

        // open or create subgroups one by one
        std::string::size_type begin = 0, end = groupName.find('/');
        while (end != std::string::npos)
        {
            std::string group(groupName.begin()+begin, groupName.begin()+end);
            hid_t prevParent = parent;

            if(H5LTfind_dataset(parent, group.c_str()) == 0)
            {
                if(create)
                    parent = H5Gcreate(prevParent, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                else
                    parent = -1;
            } 
            else 
            {
                parent = H5Gopen(prevParent, group.c_str(), H5P_DEFAULT);
            }
            H5Gclose(prevParent);

            if(parent < 0)
            {
                return parent;
            }
            begin = end + 1;
            end = groupName.find('/', begin);
        }

        return parent;
    }

        /* delete a dataset by unlinking it from the file structure. This does not
           delete the data!
         */
    inline void deleteDataset_(hid_t parent, std::string datasetName)
    {
        // delete existing data and create new dataset
        if(H5LTfind_dataset(parent, datasetName.c_str()))
        {

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

        /* get the handle of a dataset specified by a string
         */
    hid_t getDatasetHandle_(std::string datasetName) const
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        std::string groupname = SplitString(datasetName).first();
        std::string setname = SplitString(datasetName).last();

        if(H5Lexists(fileHandle_, datasetName.c_str(), H5P_DEFAULT) <= 0)
        {
            std::cerr << "HDF5File::getDatasetHandle_(): Dataset '" << datasetName << "' does not exist.\n";
            return -1;
        }

        // Open parent group
        HDF5Handle groupHandle(openGroup_(groupname), &H5Gclose, "HDF5File::getDatasetHandle_(): Internal error");

        return H5Dopen(groupHandle, setname.c_str(), H5P_DEFAULT);
    }

        /* get the type of an object specified by a string
         */
    H5O_type_t get_object_type_(std::string name) const
    {
        name = get_absolute_path(name);
        std::string group_name = SplitString(name).first();
        std::string object_name = SplitString(name).last();
        if (!object_name.size())
            return H5O_TYPE_GROUP;

        htri_t exists = H5Lexists(fileHandle_, name.c_str(), H5P_DEFAULT);
        vigra_precondition(exists > 0,  "HDF5File::get_object_type_(): "
                                        "object \"" + name + "\" "
                                        "not found.");
        // open parent group
        HDF5Handle group_handle(openGroup_(group_name), &H5Gclose, "Internal error");
        return HDF5_get_type(group_handle, name.c_str());
    }

        /* low-level write function to write vigra MultiArray data as an attribute
         */
    template<unsigned int N, class T, class Stride>
    void write_attribute_(std::string name, 
                          const std::string & attribute_name,
                          const MultiArrayView<N, T, Stride> & array,
                          const hid_t datatype, 
                          const int numBandsOfType);

        /* Write single value attribute
           This function allows to write data of atomic datatypes (int, long, double)
           as an attribute in the HDF5 file. So it is not necessary to create a MultiArray
           of size 1 to write a single number.
        */
    template<class T>
    inline void writeAtomicAttribute(std::string datasetName, std::string attributeName, const T data)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        typename MultiArrayShape<1>::type chunkSize;
        chunkSize[0] = 0;
        MultiArray<1,T> array(MultiArrayShape<1>::type(1));
        array[0] = data;
        write_attribute_(datasetName, attributeName, array, detail::getH5DataType<T>(), 1);
    }

        /* low-level read function to read vigra MultiArray data from attributes
         */
    template<unsigned int N, class T, class Stride>
    void read_attribute_(std::string datasetName, 
                         std::string attributeName, 
                         MultiArrayView<N, T, Stride> array, 
                         const hid_t datatype, const int numBandsOfType);

        /* Read a single value attribute.
           This functions allows to read a single value attribute of atomic datatype (int, long, double)
           from the HDF5 file. So it is not necessary to create a MultiArray
           of size 1 to read a single number.
        */
    template<class T>
    inline void readAtomicAttribute(std::string datasetName, std::string attributeName, T & data)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        MultiArray<1,T> array(MultiArrayShape<1>::type(1));
        read_attribute_(datasetName, attributeName, array, detail::getH5DataType<T>(), 1);
        data = array[0];
    }

    inline void readAtomicAttribute(std::string datasetName, std::string attributeName, std::string & data)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        MultiArray<1,const char *> array(MultiArrayShape<1>::type(1));
        read_attribute_(datasetName, attributeName, array, detail::getH5DataType<const char *>(), 1);
        data = std::string(array[0]);
    }

        /* low-level write function to write vigra unstrided MultiArray data
        */
    template<unsigned int N, class T, class Stride>
    void write_(std::string &datasetName, 
                       const MultiArrayView<N, T, Stride> & array, 
                       const hid_t datatype, 
                       const int numBandsOfType, 
                       typename MultiArrayShape<N>::type &chunkSize, 
                       int compressionParameter = 0);

        /* Write single value as dataset.
           This functions allows to write data of atomic datatypes (int, long, double)
           as a dataset in the HDF5 file. So it is not necessary to create a MultiArray
           of size 1 to write a single number.

           If the first character of datasetName is a "/", the path will be interpreted as absolute path,
           otherwise it will be interpreted as path relative to the current group.
        */
    template<class T>
    inline void writeAtomic(std::string datasetName, const T data)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        typename MultiArrayShape<1>::type chunkSize;
        chunkSize[0] = 0;
        MultiArray<1,T> array(MultiArrayShape<1>::type(1));
        array[0] = data;
        write_(datasetName, array, detail::getH5DataType<T>(), 1, chunkSize,0);
    }

        /* low-level read function to read vigra unstrided MultiArray data
         */
    template<unsigned int N, class T, class Stride>
    void read_(std::string datasetName, 
                      MultiArrayView<N, T, Stride> array, 
                      const hid_t datatype, const int numBandsOfType);

        /* Read a single value.
           This functions allows to read a single datum of atomic datatype (int, long, double)
           from the HDF5 file. So it is not necessary to create a MultiArray
           of size 1 to read a single number.

           If the first character of datasetName is a "/", the path will be interpreted as absolute path,
           otherwise it will be interpreted as path relative to the current group.
        */
    template<class T>
    inline void readAtomic(std::string datasetName, T & data)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        MultiArray<1,T> array(MultiArrayShape<1>::type(1));
        read_(datasetName, array, detail::getH5DataType<T>(), 1);
        data = array[0];
    }

    inline void readAtomic(std::string datasetName, std::string & data)
    {
        // make datasetName clean
        datasetName = get_absolute_path(datasetName);

        MultiArray<1,const char *> array(MultiArrayShape<1>::type(1));
        read_(datasetName, array, detail::getH5DataType<const char *>(), 1);
        data = std::string(array[0]);
    }

       /* low-level write function to write vigra unstrided MultiArray data into a sub-block of a dataset
       */
    template<unsigned int N, class T, class Stride>
    void writeBlock_(std::string datasetName, 
                     typename MultiArrayShape<N>::type &blockOffset, 
                     const MultiArrayView<N, T, Stride> & array, 
                     const hid_t datatype, 
                     const int numBandsOfType)
    {
        // open dataset if it exists
        std::string errorMessage = "HDF5File::writeBlock(): Error opening dataset '" + datasetName + "'.";
        HDF5HandleShared dataset(getDatasetHandle_(datasetName), &H5Dclose, errorMessage.c_str());
        herr_t status = writeBlock_(dataset, blockOffset, array, datatype, numBandsOfType);
        vigra_postcondition(status >= 0,
            "HDF5File::writeBlock(): write to dataset '" + datasetName + "' via H5Dwrite() failed.");
    }

       /* low-level write function to write vigra unstrided MultiArray data into a 
          sub-block of a dataset.  Returns the result of the internal call
           to <tt>H5Dwrite()</tt>.
       */
    template<unsigned int N, class T, class Stride>
    herr_t writeBlock_(HDF5HandleShared dataset, 
                       typename MultiArrayShape<N>::type &blockOffset, 
                       const MultiArrayView<N, T, Stride> & array, 
                       const hid_t datatype, 
                       const int numBandsOfType);

        /* low-level read function to read vigra unstrided MultiArray data from a sub-block of a dataset.
        
           The array must have the same shape as the block.
        */
    template<unsigned int N, class T, class Stride>
    void readBlock_(std::string datasetName, 
                    typename MultiArrayShape<N>::type &blockOffset, 
                    typename MultiArrayShape<N>::type &blockShape, 
                    MultiArrayView<N, T, Stride> &array, 
                    const hid_t datatype, const int numBandsOfType)
    {
        std::string errorMessage ("HDF5File::readBlock(): Unable to open dataset '" + datasetName + "'.");
        HDF5HandleShared dataset(getDatasetHandle_(datasetName), &H5Dclose, errorMessage.c_str());
        herr_t status = readBlock_(dataset, blockOffset, blockShape, array, datatype, numBandsOfType);
        vigra_postcondition(status >= 0,
            "HDF5File::readBlock(): read from dataset '" + datasetName + "' via H5Dread() failed.");
    }

        /* low-level read function to read vigra unstrided MultiArray data from a sub-block of a dataset.
        
           The array must have the same shape as the block. Returns the result of the internal call
           to <tt>H5Dread()</tt>.
        */
    template<unsigned int N, class T, class Stride>
    herr_t readBlock_(HDF5HandleShared dataset, 
                      typename MultiArrayShape<N>::type &blockOffset, 
                      typename MultiArrayShape<N>::type &blockShape, 
                      MultiArrayView<N, T, Stride> &array, 
                      const hid_t datatype, const int numBandsOfType);
};  /* class HDF5File */

/********************************************************************/

template<unsigned int N, class T, class Stride>
void HDF5File::write_(std::string &datasetName, 
                      const MultiArrayView<N, T, Stride> & array, 
                      const hid_t datatype, 
                      const int numBandsOfType, 
                      typename MultiArrayShape<N>::type &chunkSize, 
                      int compressionParameter)
{
    vigra_precondition(!isReadOnly(),
        "HDF5File::write(): file is read-only.");
        
    std::string groupname = SplitString(datasetName).first();
    std::string setname = SplitString(datasetName).last();

    // shape of the array. Add one dimension, if array contains non-scalars.
    ArrayVector<hsize_t> shape(array.shape().begin(), array.shape().end());
    std::reverse(shape.begin(), shape.end());

    if(numBandsOfType > 1)
        shape.push_back(numBandsOfType);

    HDF5Handle dataspace(H5Screate_simple(shape.size(), shape.begin(), NULL), &H5Sclose, 
                         "HDF5File::write(): Can not create dataspace.");

    // create and open group:
    std::string errorMessage ("HDF5File::write(): can not create group '" + groupname + "'.");
    HDF5Handle groupHandle(openCreateGroup_(groupname), &H5Gclose, errorMessage.c_str());

    // delete dataset, if it already exists
    deleteDataset_(groupHandle, setname.c_str());

    // set up properties list
    HDF5Handle plist(H5Pcreate(H5P_DATASET_CREATE), &H5Pclose, 
                     "HDF5File::write(): unable to create property list." );

    // turn off time tagging of datasets by default.
    H5Pset_obj_track_times(plist, track_time);

    // enable chunks
    ArrayVector<hsize_t> chunks(defineChunks(chunkSize, array.shape(), numBandsOfType, compressionParameter));
    if(chunks.size() > 0)
    {
        std::reverse(chunks.begin(), chunks.end());
        H5Pset_chunk (plist, chunks.size(), chunks.begin());
    }

    // enable compression
    if(compressionParameter > 0)
    {
        H5Pset_deflate(plist, compressionParameter);
    }

    // create dataset
    HDF5Handle datasetHandle(H5Dcreate(groupHandle, setname.c_str(), datatype, dataspace,H5P_DEFAULT, plist, H5P_DEFAULT), 
                             &H5Dclose, "HDF5File::write(): Can not create dataset.");

    herr_t status = 0;
    if(array.isUnstrided())
    {
        // Write the data directly from the array data buffer
        status = H5Dwrite(datasetHandle, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data());
    }
    else
    {
        // otherwise, we need an intermediate buffer
        // FIXME: right now, the buffer has the same size as the array to be read
        //        incomplete code for better solutions is below
        // MultiArray<N, T> buffer(array);
        // status = H5Dwrite(datasetHandle, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data());
        
        int offset = numBandsOfType > 1 ? 1 : 0;
        std::reverse(shape.begin(), shape.end());
        if(chunks.size() > 0)
        {
            // if the file is chunked, we use a buffer that matches the chunk size.
            std::reverse(chunks.begin(), chunks.end());
        }
        else
        {
            // otherwise, we compute a suitable chunk size.
            ArrayVector<hsize_t>(shape.size(), 1).swap(chunks);
            chunks[0] = numBandsOfType;
            MultiArrayIndex prod = 1;
            for(unsigned int k=0; k<N;  ++k)
            {
                chunks[k+offset] = array.shape(k);
                prod *= array.shape(k);
                if(prod > 300000)
                    break;
            }
        }

        ArrayVector<hsize_t> null(shape.size(), 0),
                             start(shape.size(), 0),
                             count(shape.size(), 1);
        
        count[N-1-offset] = numBandsOfType;
        
        typedef typename MultiArrayShape<N>::type Shape;
        Shape chunkCount, chunkMaxShape;
        for(unsigned int k=offset; k<chunks.size(); ++k)
        {
            chunkMaxShape[k-offset] = chunks[k];
            chunkCount[k-offset] = (MultiArrayIndex)std::ceil(double(shape[k]) / chunks[k]);
        }
        
        typename CoupledIteratorType<N>::type chunkIter = createCoupledIterator(chunkCount),
                                              chunkEnd  = chunkIter.getEndIterator();
        for(; chunkIter != chunkEnd; ++chunkIter)
        {
            Shape chunkStart(chunkIter.point() * chunkMaxShape),
                  chunkStop(min(chunkStart + chunkMaxShape, array.shape()));
            MultiArray<N, T> buffer(array.subarray(chunkStart, chunkStop));
            
            for(unsigned int k=0; k<N; ++k)
            {
                start[N-1-k] = chunkStart[k];
                count[N-1-k] = buffer.shape(k);
            }
            if(offset == 1)
            {
                start[N] = 0;
                count[N] = numBandsOfType;
            }
            HDF5Handle filespace(H5Dget_space(datasetHandle),
                                 &H5Sclose, "HDF5File::write(): unable to create hyperslabs.");
            status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start.data(), NULL, count.data(), NULL);
            if(status < 0)
                break;
                
            HDF5Handle dataspace(H5Screate_simple(count.size(), count.data(), NULL),
                                 &H5Sclose, "HDF5File::write(): unable to create hyperslabs."); 
            status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, null.data(), NULL, count.data(), NULL);
            if(status < 0)
                break;
                
            status = H5Dwrite(datasetHandle, datatype, dataspace, filespace, H5P_DEFAULT, buffer.data());
            if(status < 0)
                break;
        }
    }
    vigra_postcondition(status >= 0,
        "HDF5File::write(): write to dataset '" + datasetName + "' via H5Dwrite() failed.");
}

/********************************************************************/

template<unsigned int N, class T, class Stride>
herr_t HDF5File::writeBlock_(HDF5HandleShared datasetHandle, 
                             typename MultiArrayShape<N>::type &blockOffset, 
                             const MultiArrayView<N, T, Stride> & array, 
                             const hid_t datatype, 
                             const int numBandsOfType)
{
    vigra_precondition(!isReadOnly(),
        "HDF5File::writeBlock(): file is read-only.");
        
    ArrayVector<hsize_t> boffset, bshape, bones(N+1, 1);
    hssize_t dimensions = getDatasetDimensions_(datasetHandle);
    if(numBandsOfType > 1)
    {
        vigra_precondition(N+1 == dimensions,
            "HDF5File::readBlock(): Array dimension disagrees with data dimension.");
        bshape.resize(N+1);
        boffset.resize(N+1);
        bshape[N] = numBandsOfType;
        boffset[N] = 0;
    }
    else
    {
        vigra_precondition(N == dimensions,
            "HDF5File::readBlock(): Array dimension disagrees with data dimension.");
        bshape.resize(N);
        boffset.resize(N);
    }

    for(int i = 0; i < N; ++i)
    {
        // vigra and hdf5 use different indexing
        bshape[N-1-i] = array.shape(i);
        boffset[N-1-i] = blockOffset[i];
    }

    // create a target dataspace in memory with the shape of the desired block
    HDF5Handle memspace_handle (H5Screate_simple(bshape.size(), bshape.data(), NULL),
                                &H5Sclose,
                                "Unable to get origin dataspace");

    // get file dataspace and select the desired block
    HDF5Handle dataspaceHandle (H5Dget_space(datasetHandle),&H5Sclose,"Unable to create target dataspace");
    H5Sselect_hyperslab(dataspaceHandle, H5S_SELECT_SET, 
                        boffset.data(), bones.data(), bones.data(), bshape.data());

    herr_t status = 0;
    if(array.isUnstrided())
    {
        // when the array is unstrided, we can read the data directly from the array buffer
        status = H5Dwrite( datasetHandle, datatype, memspace_handle, dataspaceHandle, H5P_DEFAULT, array.data());
    }
    else
    {
        // otherwise, we must copy the data into an unstrided extra buffer
        MultiArray<N, T> buffer(array);
        status = H5Dwrite( datasetHandle, datatype, memspace_handle, dataspaceHandle, H5P_DEFAULT, buffer.data());
    }
    return status;
}

/********************************************************************/

template<unsigned int N, class T, class Stride>
void HDF5File::write_attribute_(std::string name, 
                                const std::string & attribute_name,
                                const MultiArrayView<N, T, Stride> & array,
                                const hid_t datatype, 
                                const int numBandsOfType)
{
    vigra_precondition(!isReadOnly(),
        "HDF5File::writeAttribute(): file is read-only.");
        
    // shape of the array. Add one dimension, if array contains non-scalars.
    ArrayVector<hsize_t> shape(array.shape().begin(), array.shape().end());
    std::reverse(shape.begin(), shape.end());
    if(numBandsOfType > 1)
        shape.push_back(numBandsOfType);

    HDF5Handle dataspace(H5Screate_simple(shape.size(),
                                          shape.begin(), NULL),
                         &H5Sclose, "HDF5File::writeAttribute(): Can not"
                                    " create dataspace.");

    std::string errorMessage ("HDF5File::writeAttribute(): can not find "
                              "object '" + name + "'.");

    H5O_type_t h5_type = get_object_type_(name);
    bool is_group = h5_type == H5O_TYPE_GROUP;
    if (!is_group && h5_type != H5O_TYPE_DATASET)
        vigra_precondition(0, "HDF5File::writeAttribute(): object \""
                               + name + "\" is neither a group nor a "
                               "dataset.");
    // get parent object handle
    HDF5Handle object_handle(is_group
                                 ? openCreateGroup_(name)
                                 : getDatasetHandle_(name),
                             is_group
                                 ? &H5Gclose
                                 : &H5Dclose,
                             errorMessage.c_str());
    // create / open attribute
    bool exists = existsAttribute(name, attribute_name);
    HDF5Handle attributeHandle(exists
                               ? H5Aopen(object_handle,
                                         attribute_name.c_str(),
                                         H5P_DEFAULT)
                               : H5Acreate(object_handle,
                                           attribute_name.c_str(), datatype,
                                           dataspace, H5P_DEFAULT,
                                           H5P_DEFAULT),
                               &H5Aclose,
                               "HDF5File::writeAttribute(): Can not create"
                               " attribute.");
    herr_t status = 0;
    if(array.isUnstrided())
    {
        // write the data directly from the array data buffer
        status = H5Awrite(attributeHandle, datatype, array.data());
    }
    else
    {
        // write the data via an unstrided copy
        // (we assume that attributes are small arrays, so that the wasted memory is uncritical)
        MultiArray<N, T> buffer(array);
        status = H5Awrite(attributeHandle, datatype, buffer.data());
    }
    vigra_postcondition(status >= 0,
        "HDF5File::writeAttribute(): write to attribute '" + attribute_name + "' via H5Awrite() failed.");
}
    
/********************************************************************/

template<unsigned int N, class T, class Stride>
void HDF5File::read_(std::string datasetName, 
                     MultiArrayView<N, T, Stride> array, 
                     const hid_t datatype, const int numBandsOfType)
{
    //Prepare to read without using HDF5ImportInfo
    ArrayVector<hsize_t> dimshape = getDatasetShape(datasetName);

    std::string errorMessage ("HDF5File::read(): Unable to open dataset '" + datasetName + "'.");
    HDF5Handle datasetHandle(getDatasetHandle_(datasetName), &H5Dclose, errorMessage.c_str());

    // the object in the HDF5 file may have one additional dimension which we 
    // interprete as the pixel type's bands
    int offset = (numBandsOfType > 1)
                    ? 1
                    : 0;

    vigra_precondition((N + offset ) == MultiArrayIndex(dimshape.size()), 
        "HDF5File::read(): Array dimension disagrees with dataset dimension.");

    typename MultiArrayShape<N>::type shape;
    for(int k=offset; k < (int)dimshape.size(); ++k)
        shape[k-offset] = (MultiArrayIndex)dimshape[k];

    vigra_precondition(shape == array.shape(),
                       "HDF5File::read(): Array shape disagrees with dataset shape.");
    if (offset)
        vigra_precondition(dimshape[0] == static_cast<hsize_t>(numBandsOfType),
                           "HDF5File::read(): Band count doesn't match destination array compound type.");

    herr_t status = 0;
    if(array.isUnstrided())
    {
        // when the array is unstrided, we can read the data directly into the array buffer
        status = H5Dread(datasetHandle, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, array.data() );
    }
    else
    {
        // otherwise, we need an intermediate buffer

        ArrayVector<hsize_t> null(dimshape.size(), 0),
                             chunks(dimshape.size(), 1),
                             start(dimshape.size(), 0),
                             count(dimshape.size(), 1);

        HDF5Handle properties(H5Dget_create_plist(datasetHandle),
                             &H5Pclose, "HDF5File::read(): failed to get property list");
        if(H5D_CHUNKED == H5Pget_layout(properties))
        {
            // if the file is chunked, we use a buffer that matches the chunk size.
            H5Pget_chunk(properties, chunks.size(), chunks.data());
            std::reverse(chunks.begin(), chunks.end());
        }
        else
        {
            // otherwise, we compute a suitable chunk size.
            chunks[0] = numBandsOfType;
            MultiArrayIndex prod = 1;
            for(unsigned int k=0; k<N; ++k)
            {
                chunks[k+offset] = array.shape(k);
                prod *= array.shape(k);
                if(prod > 300000)
                    break;
            }
        }
        
        count[N-1-offset] = numBandsOfType;
        
        typedef typename MultiArrayShape<N>::type Shape;
        Shape chunkCount, chunkMaxShape;
        for(unsigned int k=offset; k<chunks.size(); ++k)
        {
            chunkMaxShape[k-offset] = chunks[k];
            chunkCount[k-offset] = (MultiArrayIndex)std::ceil(double(dimshape[k]) / chunks[k]);
        }
        
        typename CoupledIteratorType<N>::type chunkIter = createCoupledIterator(chunkCount),
                                              chunkEnd  = chunkIter.getEndIterator();
        for(; chunkIter != chunkEnd; ++chunkIter)
        {
            Shape chunkStart(chunkIter.point() * chunkMaxShape),
                  chunkStop(min(chunkStart + chunkMaxShape, array.shape()));
            MultiArray<N, T> buffer(chunkStop - chunkStart);
            
            for(unsigned int k=0; k<N; ++k)
            {
                start[N-1-k] = chunkStart[k];
                count[N-1-k] = buffer.shape(k);
            }
            if(offset == 1)
            {
                start[N] = 0;
                count[N] = numBandsOfType;
            }
            HDF5Handle filespace(H5Dget_space(datasetHandle),
                                 &H5Sclose, "HDF5File::read(): unable to create hyperslabs.");
            status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start.data(), NULL, count.data(), NULL);
            if(status < 0)
                break;
                
            HDF5Handle dataspace(H5Screate_simple(count.size(), count.data(), NULL),
                                 &H5Sclose, "HDF5File::read(): unable to create hyperslabs."); 
            status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, null.data(), NULL, count.data(), NULL);
            if(status < 0)
                break;
                
            status = H5Dread(datasetHandle, datatype, dataspace, filespace, H5P_DEFAULT, buffer.data());
            if(status < 0)
                break;
                
            array.subarray(chunkStart, chunkStop) = buffer;
        }
    }
    vigra_postcondition(status >= 0,
        "HDF5File::read(): read from dataset '" + datasetName + "' via H5Dread() failed.");
}

/********************************************************************/

template<unsigned int N, class T, class Stride>
herr_t HDF5File::readBlock_(HDF5HandleShared datasetHandle, 
                            typename MultiArrayShape<N>::type &blockOffset, 
                            typename MultiArrayShape<N>::type &blockShape, 
                            MultiArrayView<N, T, Stride> &array, 
                            const hid_t datatype, const int numBandsOfType)
{
    vigra_precondition(blockShape == array.shape(),
         "HDF5File::readBlock(): Array shape disagrees with block size.");

    ArrayVector<hsize_t> boffset, bshape, bones(N+1, 1);
    hssize_t dimensions = getDatasetDimensions_(datasetHandle);
    if(numBandsOfType > 1)
    {
        vigra_precondition(N+1 == dimensions,
            "HDF5File::readBlock(): Array dimension disagrees with data dimension.");
        bshape.resize(N+1);
        boffset.resize(N+1);
        bshape[N] = numBandsOfType;
        boffset[N] = 0;
    }
    else
    {
        vigra_precondition(N == dimensions,
            "HDF5File::readBlock(): Array dimension disagrees with data dimension.");
        bshape.resize(N);
        boffset.resize(N);
    }

    for(int i = 0; i < N; ++i)
    {
        // vigra and hdf5 use different indexing
        bshape[N-1-i] = blockShape[i];
        boffset[N-1-i] = blockOffset[i];
    }

    // create a target dataspace in memory with the shape of the desired block
    HDF5Handle memspace_handle(H5Screate_simple(bshape.size(), bshape.data(), NULL),
                               &H5Sclose,
                               "Unable to create target dataspace");

    // get file dataspace and select the desired block
    HDF5Handle dataspaceHandle(H5Dget_space(datasetHandle), &H5Sclose, 
                               "Unable to get dataspace");
    H5Sselect_hyperslab(dataspaceHandle, H5S_SELECT_SET, 
                        boffset.data(), bones.data(), bones.data(), bshape.data());

    herr_t status = 0;
    if(array.isUnstrided())
    {
        // when the array is unstrided, we can read the data directly into the array buffer
        status = H5Dread( datasetHandle, datatype, memspace_handle, dataspaceHandle, H5P_DEFAULT, array.data());
    }
    else
    {
        // otherwise, we need an unstrided extra buffer ...
        MultiArray<N, T> buffer(array.shape());
        status = H5Dread( datasetHandle, datatype, memspace_handle, dataspaceHandle, H5P_DEFAULT, buffer.data());
        // ... and must copy the values
        if(status >= 0)
            array = buffer;
    }
    return status;
}

/********************************************************************/

template<unsigned int N, class T, class Stride>
void HDF5File::read_attribute_(std::string datasetName, 
                               std::string attributeName, 
                               MultiArrayView<N, T, Stride> array, 
                               const hid_t datatype, const int numBandsOfType)
{
    std::string dataset_path = get_absolute_path(datasetName);
    // open Attribute handle
    std::string message = "HDF5File::readAttribute(): could not get handle for attribute '"+attributeName+"'' of object '"+dataset_path+"'.";
    HDF5Handle attr_handle (H5Aopen_by_name(fileHandle_,dataset_path.c_str(),attributeName.c_str(),H5P_DEFAULT,H5P_DEFAULT),&H5Aclose, message.c_str());

    // get Attribute dataspace
    message = "HDF5File::readAttribute(): could not get dataspace for attribute '"+attributeName+"'' of object '"+dataset_path+"'.";
    HDF5Handle attr_dataspace_handle (H5Aget_space(attr_handle),&H5Sclose,message.c_str());

    // obtain Attribute shape
    int raw_dims = H5Sget_simple_extent_ndims(attr_dataspace_handle);
    int dims = std::max(raw_dims, 1); // scalar attributes may be stored with raw_dims==0
    ArrayVector<hsize_t> dimshape(dims);
    if(raw_dims > 0)
        H5Sget_simple_extent_dims(attr_dataspace_handle, dimshape.data(), NULL);
    else
        dimshape[0] = 1;
    
    // invert the dimensions to guarantee VIGRA-compatible order
    std::reverse(dimshape.begin(), dimshape.end());

    int offset = (numBandsOfType > 1)
                    ? 1
                    : 0;
    message = "HDF5File::readAttribute(): Array dimension disagrees with dataset dimension.";
    // the object in the HDF5 file may have one additional dimension which we then interpret as the pixel type bands
    vigra_precondition((N + offset) == MultiArrayIndex(dims), message);

    for(int k=offset; k < (int)dimshape.size(); ++k)
        vigra_precondition(array.shape()[k-offset] == (MultiArrayIndex)dimshape[k],
                           "HDF5File::readAttribute(): Array shape disagrees with dataset shape");

    herr_t status = 0;
    if(array.isUnstrided())
    {
        // when the array is unstrided, we can read the data directly into the array buffer
        status = H5Aread( attr_handle, datatype, array.data());
    }
    else
    {
        // otherwise, we need an unstrided extra buffer ...
        // (we assume that attributes are small arrays, so that the wasted memory is uncritical)
        MultiArray<N, T> buffer(array.shape());
        status = H5Aread( attr_handle, datatype, buffer.data() );
        // ... and must copy the values
        if(status >= 0)
            array = buffer;
    }
    vigra_postcondition(status >= 0,
        "HDF5File::readAttribute(): read from attribute '" + attributeName + "' via H5Aread() failed.");
}

/********************************************************************/

/** \brief Read the data specified by the given \ref vigra::HDF5ImportInfo object
                and write the into the given 'array'.
                
    The array must have the correct number of dimensions and shape for the dataset 
    represented by 'info'. When the element type of 'array' differs from the stored element
    type, HDF5 will convert the type on the fly (except when the HDF5 version is 1.6 or below,
    in which case an error will result). Multi-channel element types (i.e. \ref vigra::RGBValue,
    \ref vigra::TinyVector, and \ref vigra::FFTWComplex) are recognized and handled correctly.
    
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

template<unsigned int N, class T, class StrideTag>
inline void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, T, StrideTag> array)
{
    readHDF5(info, array, 0, 0); // last two arguments are not used
}

template<unsigned int N, class T, class StrideTag>
void readHDF5(const HDF5ImportInfo &info, MultiArrayView<N, T, StrideTag> array, const hid_t datatype, const int numBandsOfType)
{
    HDF5File file(info.getFilePath(), HDF5File::OpenReadOnly);
    file.read(info.getPathInFile(), array);
    file.close();
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




/** \brief Store array data in an HDF5 file.
                
    The number of dimensions, shape and element type of the stored dataset is automatically 
    determined from the properties of the given \a array. Strided arrays are stored in an
    unstrided way, i.e. in contiguous scan-order. Multi-channel element types 
    (i.e. \ref vigra::RGBValue, \ref vigra::TinyVector and \ref vigra::FFTWComplex)
    are recognized and handled correctly
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

template<unsigned int N, class T, class StrideTag>
inline void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, T, StrideTag> & array)
{
    //last two arguments are not used
    writeHDF5(filePath, pathInFile, array, 0, 0);
}

template<unsigned int N, class T, class StrideTag>
void writeHDF5(const char* filePath, const char* pathInFile, const MultiArrayView<N, T, StrideTag> & array, const hid_t datatype, const int numBandsOfType)
{
    HDF5File file(filePath, HDF5File::Open);
    file.write(pathInFile, array);
    file.close();
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
