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
#include <cstring>
#include <cstdio>

namespace vigra {

namespace detail
{

int get_number_of_digits(int in)
{
    int num = 0;
    int i = 1; 
    while(double(in) / double(i) >= 1)
    {
        i *= 10;
        num += 1; 
    }
    if(num == 0)
        num = 1;
    return num; 
}

std::string make_padded_number(int number, int max_number)
{
    int max_digit_ct = get_number_of_digits(max_number);
    char buffer [50];
    std::sprintf(buffer, "%d", number);
    std::string padding = "";
    std::string numeral = buffer;
    int digit_ct = get_number_of_digits(number); 
    for(int gg = 0; gg < max_digit_ct - digit_ct; ++ gg)
        padding = padding + "0";
    return padding + numeral;
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

void options_export_HDF5(hid_t & group_id, 
                             RandomForestOptions const & opt, 
                             std::string name)
{
    hid_t opt_id = H5Gcreate(group_id, name.c_str(), 
                                           H5P_DEFAULT, 
                                           H5P_DEFAULT, 
                                        H5P_DEFAULT);
    vigra_postcondition(opt_id >= 0, 
                        "problemspec_export_HDF5():"
                        " Unable to create external opteters");

    //get a map containing all the double fields
    std::map<std::string, ArrayVector<double> > serialized_opt;
    opt.make_map(serialized_opt);
    std::map<std::string, ArrayVector<double> >::iterator iter;
    for(iter = serialized_opt.begin(); iter != serialized_opt.end(); ++iter)
        write_array_2_hdf5(opt_id, iter->second, iter->first, H5T_NATIVE_DOUBLE);
    H5Gclose(opt_id);
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


void dt_export_HDF5(    hid_t & group_id,
                            detail::DecisionTree const & tree,
                            std::string name)
{
    //check if ext_param was written and write it if not
    hid_t e_id = H5Gopen (group_id, 
                          "_ext_param", 
                          H5P_DEFAULT);
       if(e_id < 0)
    {
        problemspec_export_HDF5(group_id,
                                tree.ext_param_, 
                                "_ext_param"); 
    }
    else H5Gclose(e_id);
    
    //make the folder for the tree.
    hid_t tree_id =    H5Gcreate(group_id, name.c_str(), 
                                           H5P_DEFAULT, 
                                           H5P_DEFAULT, 
                                        H5P_DEFAULT);
    //write down topology
    write_array_2_hdf5(tree_id, 
                       tree.topology_, 
                       "topology", 
                       H5T_NATIVE_INT);
    //write down parameters
    write_array_2_hdf5(tree_id, 
                       tree.parameters_, 
                       "parameters", 
                       H5T_NATIVE_DOUBLE);
    H5Gclose(tree_id);
}

}} // namespace vigra::detail

#endif // HasHDF5
