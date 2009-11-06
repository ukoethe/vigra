/************************************************************************/
/*                                                                      */
/*       Copyright 2009 by Rahul Nair and  Ullrich Koethe               */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#include "random_forest.hxx"
#include <cstdio>
#include <string>
#include <hdf5.h>

#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
	#define H5Gcreate(a,b,c,d,e) H5Gcreate(a,b, 1);
	#define H5Gopen(a,b,c) H5Gopen(a,b)
	# include <H5LT.h>
#else
	# include <hdf5_hl.h>
#endif

namespace vigra 
{

	



namespace detail
{

/** shallow search the hdf5 group for containing elements
 * returns negative value if unsuccessful
 * \param filename name of hdf5 file
 * \param groupname path in hdf5 file
 * \param cont		reference to container that supports
 * 					insert(). valuetype of cont must be
 * 					std::string
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
	herr_t 		status;	
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
		ssize_t buffer_size =
				H5Lget_name_by_idx(grp_id, ".",
								   H5_INDEX_NAME,
								   H5_ITER_INC,
								   ii, 0, 0, H5P_DEFAULT)+1;
#endif
		char* buffer = new char[buffer_size];
#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 6)
		buffer_size = 
				H5Gget_objname_by_idx(grp_id, 
									  ii, buffer, 
									  (size_t)buffer_size );
#else
		buffer_size =
				H5Lget_name_by_idx(grp_id, ".",
								   H5_INDEX_NAME,
								   H5_ITER_INC,
								   ii, buffer,
								   (size_t)buffer_size,
								   H5P_DEFAULT);
#endif
		cont.insert(cont.end(), std::string(buffer));
		delete [] buffer;
	}
	return 1;
}


template<class Container>
bool find_groups_hdf5(std::string filename, 
							  std::string groupname, 
							  Container &cont)
{
	//check if file exists
	FILE* pFile;
	pFile = fopen ( filename.c_str(), "r" );
	if ( pFile == NULL)
	{	
		return 0;
	}
	//open the file
	hid_t file_id =  H5Fopen(filename.c_str(),
				  			 H5F_ACC_RDONLY,
				  			 H5P_DEFAULT);
	hid_t grp_id = groupname == "" ?
		file_id
		:H5Gopen(file_id, groupname.c_str(), H5P_DEFAULT);

	bool res =  find_groups_hdf5(grp_id, cont); 
	if(groupname != "")
		H5Gclose(grp_id);
	H5Fclose(file_id);
	return res; 
}

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
	arr.resize(size);
	vigra_postcondition(H5LTread_dataset (id, 
										  name.c_str(),
										  type, 
										  arr.data()) >= 0,
						"write_array_2_hdf5():"
						"unable to read dataset");
}


void options_import_HDF5(hid_t & group_id,
						 RandomForestOptions & opt, 
						 std::string name)
{
	ArrayVector<double> serialized_options;
	write_hdf5_2_array(group_id, serialized_options,
					  	name, H5T_NATIVE_DOUBLE); 
	opt.unserialize(serialized_options.begin(),
				  	serialized_options.end());
}

void options_export_HDF5(hid_t & group_id,
						 RandomForestOptions const & opt, 
						 std::string name)
{
	hsize_t size = opt.serialized_size();
	ArrayVector<double> serialized_options(opt.serialized_size());
	opt.serialize(serialized_options.begin(),
				  serialized_options.end());
	write_array_2_hdf5(group_id, serialized_options,
					  name, H5T_NATIVE_DOUBLE); 
}

void problemspec_import_HDF5(hid_t & group_id, 
							 ProblemSpec  & param, 
							 std::string name)
{
	hid_t param_id = H5Gopen (group_id, 
							  name.c_str(), 
							  H5P_DEFAULT);

	vigra_postcondition(param_id >= 0, 
						"problemspec_export_HDF5():"
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
	switch(param.class_type_)
	{
		#define SOME_CASE(type_, hdf5_type_) \
		case ProblemSpec::type_##_t:\
		{\
			ArrayVector<type_> tmp;\
			write_hdf5_2_array(param_id, tmp, "labels", hdf5_type_);\
			param.classes_(tmp.begin(), tmp.end());\
		}\
			break;
		SOME_CASE(UInt8, 	H5T_NATIVE_UINT8);
		SOME_CASE(UInt16, 	H5T_NATIVE_UINT16);
		SOME_CASE(UInt32, 	H5T_NATIVE_UINT32);
		SOME_CASE(UInt64, 	H5T_NATIVE_UINT64);
		SOME_CASE(Int8,  	H5T_NATIVE_INT8);
		SOME_CASE(Int16, 	H5T_NATIVE_INT16);
		SOME_CASE(Int32, 	H5T_NATIVE_INT32);
		SOME_CASE(Int64, 	H5T_NATIVE_INT64);
		SOME_CASE(double, 	H5T_NATIVE_DOUBLE);
		SOME_CASE(float, 	H5T_NATIVE_FLOAT);
		default:
			std::runtime_error("exportRF_HDF5(): unknown class type"); 
		#undef SOME_CASE
	}
	H5Gclose(param_id);
}

void problemspec_export_HDF5(hid_t & group_id, 
							 ProblemSpec const & param, 
							 std::string name)
{
	hsize_t		size = 1;
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
	switch(param.class_type_)
	{
		#define SOME_CASE(type_, hdf5_type_) \
		case ProblemSpec::type_##_t:\
			write_array_2_hdf5(param_id, param.type_##_classes_, "labels", hdf5_type_);\
			break;
		SOME_CASE(UInt8, 	H5T_NATIVE_UINT8);
		SOME_CASE(UInt16, 	H5T_NATIVE_UINT16);
		SOME_CASE(UInt32, 	H5T_NATIVE_UINT32);
		SOME_CASE(UInt64, 	H5T_NATIVE_UINT64);
		SOME_CASE(Int8,  	H5T_NATIVE_INT8);
		SOME_CASE(Int16, 	H5T_NATIVE_INT16);
		SOME_CASE(Int32, 	H5T_NATIVE_INT32);
		SOME_CASE(Int64, 	H5T_NATIVE_INT64);
		SOME_CASE(double, 	H5T_NATIVE_DOUBLE);
		SOME_CASE(float, 	H5T_NATIVE_FLOAT);
		case ProblemSpec::UNKNOWN:
			break;
		default:
			std::runtime_error("exportRF_HDF5(): unknown class type"); 
		#undef SOME_CASE
	}
	H5Gclose(param_id);
}

void dt_import_HDF5(hid_t & group_id,
				RF_Traits::DecisionTree_t & tree,
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


void dt_export_HDF5(hid_t & group_id,
				RF_Traits::DecisionTree_t const & tree,
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
	hid_t tree_id =	H5Gcreate(group_id, name.c_str(), 
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
} //namespace detail

template<class T>
bool rf_export_HDF5(RandomForest<T> &rf, 
					std::string filename, 
					std::string pathname = "",
					bool overwriteflag = false)
{ 
	using detail::make_padded_number;
	using detail::options_export_HDF5;
	using detail::problemspec_export_HDF5;
	using detail::dt_export_HDF5;
	//if file exists delete it.
	FILE* pFile = fopen ( filename.c_str(), "r" );
	if ( pFile != NULL && !overwriteflag)
		return 0;
	else if(pFile != 0 &&std::remove(filename.c_str()) != 0)
		return 0;
	
	//create a new file and group.
	hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, 
												H5P_DEFAULT, 
												H5P_DEFAULT);
	vigra_postcondition(file_id >= 0, 
						"rf_export_HDF5(): Unable to open file.");
    hid_t group_id = pathname== "" ?
						file_id
					:	H5Gcreate(file_id, pathname.c_str(), 
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
		std::string treename = 	"Tree_"  + 
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
	//if file exists delete it.
	FILE* pFile = fopen ( filename.c_str(), "r" );
	if ( pFile == NULL)
		return 0;
	
	//open file
    hid_t file_id = H5Fopen (filename.c_str(), 
							 H5F_ACC_RDONLY, 
							 H5P_DEFAULT);
	
	vigra_postcondition(file_id >= 0, 
						"rf_import_HDF5(): Unable to open file.");
    hid_t group_id = pathname== "" ?
						file_id
					:	H5Gopen (file_id, 
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
			rf.trees_.push_back(RF_Traits::DecisionTree_t(rf.ext_param_));
			dt_import_HDF5(group_id, rf.trees_.back(), *iter); 
		}
	}
	
	//clean up the mess
	if(pathname != "")
		H5Gclose(group_id);
	H5Fclose(file_id);

	return 1;
}
} // namespace vigra


#endif // VIGRA_RANDOM_FOREST_HDF5_IMPEX_HXX

