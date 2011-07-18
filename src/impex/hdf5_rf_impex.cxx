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

#include "vigra/random_forest_hdf5_impex.hxx"
#include "vigra/multi_array.hxx"
#include <iostream>
#include <cstring> // rm...
#include <cstdio> // rm...

#include <string>
#include <sstream>
#include <iomanip>

namespace vigra {

namespace detail {

struct padded_number_string_data : public std::ostringstream
{
    unsigned w;
};
padded_number_string::padded_number_string(int n)
    : padded_number(new padded_number_string_data())
{
    (*padded_number) << (n - 1);
    padded_number->w = padded_number->str().size();
}
std::string padded_number_string::operator()(int k) const
{
    padded_number->str("");
    (*padded_number) << std::setw(padded_number->w) << std::setfill('0') << k;
    return padded_number->str();
}
padded_number_string::~padded_number_string()
{
    delete padded_number;
}

MyT::type type_of_hid_t(hid_t group_id, std::string name)
{
    hid_t m_dataset_handle = 
    H5Dopen(group_id, name.c_str(), H5P_DEFAULT);
    hid_t datatype = H5Dget_type(m_dataset_handle);
    H5T_class_t dataclass = H5Tget_class(datatype);
    size_t datasize  = H5Tget_size(datatype);
    H5T_sign_t datasign  = H5Tget_sign(datatype);
    MyT::type result = MyT::OTHER; 
    if(dataclass == H5T_FLOAT)
    {
        if(datasize == 4)
            result = MyT::FLOAT;
        else if(datasize == 8)
            result = MyT::DOUBLE;
    }
    else if(dataclass == H5T_INTEGER)    
    {
        if(datasign == H5T_SGN_NONE)
        {
            if(datasize ==  1)
                result = MyT::UINT8;
            else if(datasize == 2)
                result = MyT::UINT16;
            else if(datasize == 4)
                result = MyT::UINT32;
            else if(datasize == 8)
                result = MyT::UINT64;
        }
        else
        {
            if(datasize ==  1)
                result = MyT::INT8;
            else if(datasize == 2)
                result = MyT::INT16;
            else if(datasize == 4)
                result = MyT::INT32;
            else if(datasize == 8)
                result = MyT::INT64;
        }
    }
    H5Tclose(datatype);
    H5Dclose(m_dataset_handle);
    return result;
}

    
void options_import_HDF5(hid_t & group_id, 
                        RandomForestOptions  & opt, 
                        std::string name)
{
    hid_t opt_id = H5Gopen (group_id, 
                              name.c_str(), 
                              H5P_DEFAULT);

    vigra_postcondition(opt_id >= 0, 
                        "problemspec_import_HDF5():"
                        " Unable to open external opteters");

    //get a map containing all the double fields
    std::set<std::string> ext_set;
    find_groups_hdf5(opt_id, ext_set);
    std::map<std::string, ArrayVector <double> > ext_map;
    std::set<std::string>::iterator iter;
    if(ext_set.find(std::string("labels")) == ext_set.end())
        std::runtime_error("labels are missing");
    for(iter = ext_set.begin(); iter != ext_set.end(); ++ iter)
    {
        if(*iter != std::string("labels"))
        {
            ext_map[*iter] = ArrayVector<double>();
            write_hdf5_2_array(opt_id, ext_map[*iter], 
                               *iter, H5T_NATIVE_DOUBLE);
        }
    }
    opt.make_from_map(ext_map);
    H5Gclose(opt_id);
}

void options_export_HDF5(HDF5File & h5context, RandomForestOptions const & opt,
                         std::string name)
{
    h5context.cd_mk(name);
    rf_export_map_to_HDF5(h5context, opt);
    h5context.cd_up();
}

void dt_import_HDF5(    hid_t & group_id,
                             detail::DecisionTree & tree,
                            std::string name)
{
    //check if ext_param was written and write it if not
       if(tree.ext_param_.actual_msample_ == 0)
    {
        problemspec_import_HDF5(group_id, tree.ext_param_, "_ext_param");
        tree.classCount_ = tree.ext_param_.class_count_;
    }
    
    hid_t tree_id =H5Gopen (group_id, name.c_str(), H5P_DEFAULT);
    //write down topology
    write_hdf5_2_array(tree_id, 
                       tree.topology_, 
                       "topology", 
                       H5T_NATIVE_INT);
    //write down parameters
    write_hdf5_2_array(tree_id, 
                       tree.parameters_, 
                       "parameters", 
                       H5T_NATIVE_DOUBLE);
    H5Gclose(tree_id);
}

void dt_export_HDF5(HDF5File & h5context,
                    detail::DecisionTree const & tree,
                    std::string name)
{
    // make the folder for the tree.
    h5context.cd_mk(name);
    // write down topology
    h5context.write("topology", tree.topology_);
    // write down parameters
    h5context.write("parameters", tree.parameters_);
    h5context.cd_up();
}

} // namespace detail
} // namespace vigra

#endif // HasHDF5
