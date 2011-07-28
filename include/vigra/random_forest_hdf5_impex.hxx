/************************************************************************/
/*                                                                      */
/*       Copyright 2009 by Rahul Nair and  Ullrich Koethe               */
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


#ifndef VIGRA_RANDOM_FOREST_IMPEX_HDF5_HXX
#define VIGRA_RANDOM_FOREST_IMPEX_HDF5_HXX

#include "config.hxx"
#include "random_forest.hxx"
#include "hdf5impex.hxx"
#include <cstdio>
#include <string>

#ifdef HasHDF5

namespace vigra 
{

namespace detail
{


/** shallow search the hdf5 group for containing elements
 * returns negative value if unsuccessful
 * \param grp_id    hid_t containing path to group.
 * \param cont        reference to container that supports
 *                     insert(). valuetype of cont must be
 *                     std::string
 */
template<class Container>
bool find_groups_hdf5(hid_t grp_id, Container &cont)
{
    
    //get group info
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
    hsize_t size;
    H5Gget_num_objs(grp_id, &size);
#else
    hsize_t size;
    H5G_info_t ginfo;
    herr_t         status;    
    status = H5Gget_info (grp_id , &ginfo);
    if(status < 0)
        std::runtime_error("find_groups_hdf5():"
                           "problem while getting group info");
    size = ginfo.nlinks;
#endif
    for(hsize_t ii = 0; ii < size; ++ii)
    {
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
        ssize_t buffer_size = 
                H5Gget_objname_by_idx(grp_id, 
                                      ii, NULL, 0 ) + 1;
#else
        std::ptrdiff_t buffer_size =
                H5Lget_name_by_idx(grp_id, ".",
                                   H5_INDEX_NAME,
                                   H5_ITER_INC,
                                   ii, 0, 0, H5P_DEFAULT)+1;
#endif
        ArrayVector<char> buffer(buffer_size);
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
        buffer_size = 
                H5Gget_objname_by_idx(grp_id, 
                                      ii, buffer.data(), 
                                      (size_t)buffer_size );
#else
        buffer_size =
                H5Lget_name_by_idx(grp_id, ".",
                                   H5_INDEX_NAME,
                                   H5_ITER_INC,
                                   ii, buffer.data(),
                                   (size_t)buffer_size,
                                   H5P_DEFAULT);
#endif
        cont.insert(cont.end(), std::string(buffer.data()));
    }
    return true;
}


/** shallow search the hdf5 group for containing elements
 * returns negative value if unsuccessful
 * \param filename name of hdf5 file
 * \param groupname path in hdf5 file
 * \param cont        reference to container that supports
 *                     insert(). valuetype of cont must be
 *                     std::string
 */
template<class Container>
bool find_groups_hdf5(std::string filename, 
                              std::string groupname, 
                              Container &cont)
{
    //check if file exists
    FILE* pFile;
    pFile = std::fopen ( filename.c_str(), "r" );
    if ( pFile == NULL)
    {    
        return 0;
    }
    std::fclose(pFile);
    //open the file
    HDF5Handle file_id(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
                       &H5Fclose, "Unable to open HDF5 file");
    HDF5Handle grp_id;
    if(groupname == "")
    {
        grp_id = HDF5Handle(file_id, 0, "");
    }
    else
    {
        grp_id = HDF5Handle(H5Gopen(file_id, groupname.c_str(), H5P_DEFAULT),
                            &H5Gclose, "Unable to open group");

    }
    bool res =  find_groups_hdf5(grp_id, cont); 
    return res; 
}

VIGRA_EXPORT int get_number_of_digits(int in);

VIGRA_EXPORT std::string make_padded_number(int number, int max_number);

/** write a ArrayVector to a hdf5 dataset.
 */
template<class U, class T>
void write_array_2_hdf5(hid_t & id, 
                        ArrayVector<U> const & arr, 
                        std::string    const & name, 
                        T  type) 
{
    hsize_t size = arr.size(); 
    vigra_postcondition(H5LTmake_dataset (id, 
                                          name.c_str(), 
                                          1, 
                                          &size, 
                                          type, 
                                          arr.begin()) 
                        >= 0,
                        "write_array_2_hdf5():"
                        "unable to write dataset");
}


template<class U, class T>
void write_hdf5_2_array(hid_t & id, 
                        ArrayVector<U>       & arr, 
                        std::string    const & name, 
                        T  type) 
{    
    // The last three values of get_dataset_info can be NULL
    // my EFFING FOOT! that is valid for HDF5 1.8 but not for
    // 1.6 - but documented the other way around AAARRHGHGHH
    hsize_t size; 
    H5T_class_t a; 
    size_t b;
    vigra_postcondition(H5LTget_dataset_info(id, 
                                             name.c_str(), 
                                             &size, 
                                             &a, 
                                             &b) >= 0,
                        "write_hdf5_2_array(): "
                        "Unable to locate dataset");
    arr.resize((typename ArrayVector<U>::size_type)size);
    vigra_postcondition(H5LTread_dataset (id, 
                                          name.c_str(),
                                          type, 
                                          arr.data()) >= 0,
                        "write_array_2_hdf5():"
                        "unable to read dataset");
}

/*
inline void options_import_HDF5(hid_t & group_id,
                         RandomForestOptions & opt, 
                         std::string name)
{
    ArrayVector<double> serialized_options;
    write_hdf5_2_array(group_id, serialized_options,
                          name, H5T_NATIVE_DOUBLE); 
    opt.unserialize(serialized_options.begin(),
                      serialized_options.end());
}

inline void options_export_HDF5(hid_t & group_id,
                         RandomForestOptions const & opt, 
                         std::string name)
{
    ArrayVector<double> serialized_options(opt.serialized_size());
    opt.serialize(serialized_options.begin(),
                  serialized_options.end());
    write_array_2_hdf5(group_id, serialized_options,
                      name, H5T_NATIVE_DOUBLE); 
}
*/

struct MyT
{
    enum type { INT8 = 1,  INT16 = 2,  INT32 =3,  INT64=4, 
                  UINT8 = 5, UINT16 = 6, UINT32= 7, UINT64= 8,
                  FLOAT = 9, DOUBLE = 10, OTHER = 3294};
};



#define create_type_of(TYPE, ENUM) \
inline MyT::type type_of(TYPE)\
{\
    return MyT::ENUM; \
}
create_type_of(Int8, INT8)
create_type_of(Int16, INT16)
create_type_of(Int32, INT32)
create_type_of(Int64, INT64)
create_type_of(UInt8, UINT8)
create_type_of(UInt16, UINT16)
create_type_of(UInt32, UINT32)
create_type_of(UInt64, UINT64)
create_type_of(float, FLOAT)
create_type_of(double, DOUBLE)
#undef create_type_of

VIGRA_EXPORT MyT::type type_of_hid_t(hid_t group_id, std::string name);

VIGRA_EXPORT void options_import_HDF5(hid_t & group_id, 
                        RandomForestOptions  & opt, 
                        std::string name);

VIGRA_EXPORT void options_export_HDF5(hid_t & group_id, 
                         RandomForestOptions const & opt, 
                         std::string name);

template<class T>
void problemspec_import_HDF5(hid_t & group_id, 
                             ProblemSpec<T>  & param, 
                             std::string name)
{
    hid_t param_id = H5Gopen (group_id, 
                              name.c_str(), 
                              H5P_DEFAULT);

    vigra_postcondition(param_id >= 0, 
                        "problemspec_import_HDF5():"
                        " Unable to open external parameters");

    //get a map containing all the double fields
    std::set<std::string> ext_set;
    find_groups_hdf5(param_id, ext_set);
    std::map<std::string, ArrayVector <double> > ext_map;
    std::set<std::string>::iterator iter;
    if(ext_set.find(std::string("labels")) == ext_set.end())
        std::runtime_error("labels are missing");
    for(iter = ext_set.begin(); iter != ext_set.end(); ++ iter)
    {
        if(*iter != std::string("labels"))
        {
            ext_map[*iter] = ArrayVector<double>();
            write_hdf5_2_array(param_id, ext_map[*iter], 
                               *iter, H5T_NATIVE_DOUBLE);
        }
    }
    param.make_from_map(ext_map);
    //load_class_labels
    switch(type_of_hid_t(param_id,"labels" ))
    {
        #define SOME_CASE(type_, enum_) \
      case MyT::enum_ :\
        {\
            ArrayVector<type_> tmp;\
            write_hdf5_2_array(param_id, tmp, "labels", H5T_NATIVE_##enum_);\
            param.classes_(tmp.begin(), tmp.end());\
        }\
            break;
        SOME_CASE(UInt8,     UINT8);
        SOME_CASE(UInt16,     UINT16);
        SOME_CASE(UInt32,     UINT32);
        SOME_CASE(UInt64,     UINT64);
        SOME_CASE(Int8,      INT8);
        SOME_CASE(Int16,     INT16);
        SOME_CASE(Int32,     INT32);
        SOME_CASE(Int64,     INT64);
        SOME_CASE(double,     DOUBLE);
        SOME_CASE(float,     FLOAT);
        default:
            std::runtime_error("exportRF_HDF5(): unknown class type"); 
        #undef SOME_CASE
    }
    H5Gclose(param_id);
}

template<class T>
void problemspec_export_HDF5(hid_t & group_id, 
                             ProblemSpec<T> const & param, 
                             std::string name)
{
    hid_t param_id = H5Gcreate(group_id, name.c_str(), 
                                           H5P_DEFAULT, 
                                           H5P_DEFAULT, 
                                        H5P_DEFAULT);
    vigra_postcondition(param_id >= 0, 
                        "problemspec_export_HDF5():"
                        " Unable to create external parameters");

    //get a map containing all the double fields
    std::map<std::string, ArrayVector<double> > serialized_param;
    param.make_map(serialized_param);
    std::map<std::string, ArrayVector<double> >::iterator iter;
    for(iter = serialized_param.begin(); iter != serialized_param.end(); ++iter)
        write_array_2_hdf5(param_id, iter->second, iter->first, H5T_NATIVE_DOUBLE);
    
    //save class_labels
    switch(type_of(param.classes[0]))
    {
        #define SOME_CASE(type) \
        case MyT::type:\
            write_array_2_hdf5(param_id, param.classes, "labels", H5T_NATIVE_##type);\
            break;
        SOME_CASE(UINT8);
        SOME_CASE(UINT16);
        SOME_CASE(UINT32);
        SOME_CASE(UINT64);
        SOME_CASE(INT8);
        SOME_CASE(INT16);
        SOME_CASE(INT32);
        SOME_CASE(INT64);
        SOME_CASE(DOUBLE);
        SOME_CASE(FLOAT);
        default:
            std::runtime_error("exportRF_HDF5(): unknown class type"); 
        #undef SOME_CASE
    }
    H5Gclose(param_id);
}

VIGRA_EXPORT void dt_import_HDF5(hid_t & group_id,
                    detail::DecisionTree & tree,
                    std::string name);


VIGRA_EXPORT void dt_export_HDF5(hid_t & group_id,
                    detail::DecisionTree const & tree,
                    std::string name);
                    
} //namespace detail

template<class T>
bool rf_export_HDF5(RandomForest<T> const &rf, 
                    std::string filename, 
                    std::string pathname = "",
                    bool overwriteflag = false)
{ 
    using detail::make_padded_number;
    using detail::options_export_HDF5;
    using detail::problemspec_export_HDF5;
    using detail::dt_export_HDF5;
    
    hid_t file_id;
    //if file exists load it.
    FILE* pFile = std::fopen ( filename.c_str(), "r" );
    if ( pFile != NULL)
    {    
        std::fclose(pFile);
        file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, 
                                                H5P_DEFAULT);
    }
    else
    {
        //create a new file.
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, 
                                                    H5P_DEFAULT, 
                                                    H5P_DEFAULT);
    }
    vigra_postcondition(file_id >= 0, 
                        "rf_export_HDF5(): Unable to open file.");
    //std::cerr << pathname.c_str()

    //if the group already exists this will cause an error
    //we will have to use the overwriteflag to check for 
    //this, but i dont know how to delete groups...

    hid_t group_id = pathname== "" ?
                        file_id
                    :    H5Gcreate(file_id, pathname.c_str(), 
                                              H5P_DEFAULT, 
                                              H5P_DEFAULT, 
                                           H5P_DEFAULT);

    vigra_postcondition(group_id >= 0, 
                        "rf_export_HDF5(): Unable to create group");

    //save serialized options
        options_export_HDF5(group_id, rf.options(), "_options"); 
    //save external parameters
        problemspec_export_HDF5(group_id, rf.ext_param(), "_ext_param");
    //save trees
    
    int tree_count = rf.options_.tree_count_;
    for(int ii = 0; ii < tree_count; ++ii)
    {
        std::string treename =     "Tree_"  + 
                                make_padded_number(ii, tree_count -1);
        dt_export_HDF5(group_id, rf.tree(ii), treename); 
    }
    
    //clean up the mess
    if(pathname != "")
        H5Gclose(group_id);
    H5Fclose(file_id);

    return 1;
}


template<class T>
bool rf_import_HDF5(RandomForest<T> &rf, 
                    std::string filename, 
                    std::string pathname = "")
{ 
    using detail::find_groups_hdf5;
    using detail::options_import_HDF5;
    using detail::problemspec_import_HDF5;
    using detail::dt_export_HDF5;
    // check if file exists
    FILE* pFile = std::fopen ( filename.c_str(), "r" );
    if ( pFile == NULL)
        return 0;
    std::fclose(pFile);
    //open file
    hid_t file_id = H5Fopen (filename.c_str(), 
                             H5F_ACC_RDONLY, 
                             H5P_DEFAULT);
    
    vigra_postcondition(file_id >= 0, 
                        "rf_import_HDF5(): Unable to open file.");
    hid_t group_id = pathname== "" ?
                        file_id
                    :    H5Gopen (file_id, 
                                 pathname.c_str(), 
                                 H5P_DEFAULT);
    
    vigra_postcondition(group_id >= 0, 
                        "rf_export_HDF5(): Unable to create group");

    //get serialized options
        options_import_HDF5(group_id, rf.options_, "_options"); 
    //save external parameters
        problemspec_import_HDF5(group_id, rf.ext_param_, "_ext_param");
    // TREE SAVING TIME
    // get all groups in base path
    
    std::set<std::string> tree_set;
    std::set<std::string>::iterator iter; 
    find_groups_hdf5(filename, pathname, tree_set);
    
    for(iter = tree_set.begin(); iter != tree_set.end(); ++iter)
    {
        if((*iter)[0] != '_')
        {
            rf.trees_.push_back(detail::DecisionTree(rf.ext_param_));
            dt_import_HDF5(group_id, rf.trees_.back(), *iter); 
        }
    }
    
    //clean up the mess
    if(pathname != "")
        H5Gclose(group_id);
    H5Fclose(file_id);
    /*rf.tree_indices_.resize(rf.tree_count());
    for(int ii = 0; ii < rf.tree_count(); ++ii)
        rf.tree_indices_[ii] = ii; */
    return 1;
}
} // namespace vigra

#endif // HasHDF5

#endif // VIGRA_RANDOM_FOREST_HDF5_IMPEX_HXX

