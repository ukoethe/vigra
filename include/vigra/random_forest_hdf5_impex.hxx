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


#ifndef VIGRA_RANDOM_FOREST_IMPEX_HXX
#define VIGRA_RANDOM_FOREST_IMPEX_HXX

#include "random_forest.hxx"
#include "hdf5impex.hxx"
#include <cstdio>

namespace vigra 
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
bool find_groups_hdf5(std::string filename, 
							  std::string groupname, 
							  Container &cont)
{
	//check if file exists
	FILE* pFile;
	pFile = fopen ( filename.c_str(), "r" );
	if ( pFile == NULL)
	{	
		return -1;
	}
	
	//open the file
	hid_t file_id =  H5Fopen(filename.c_str(),
				  			 H5F_ACC_RDONLY,
				  			 H5P_DEFAULT);
	//get group info
	herr_t 		status;	
	H5G_info_t 	finfo;
	status = H5Gget_info (file_id , &finfo);
	if(status < 0)
		std::runtime_error("find_groups_hdf5():"
						   "problem while getting file info");

	hid_t grp_id = groupname == "" ?
		file_id
		:H5Gopen(file_id, groupname.c_str(), H5P_DEFAULT);
	H5G_info_t ginfo;
	status = H5Gget_info (grp_id , &ginfo);

	if(status < 0)
		std::runtime_error("find_groups_hdf5():"
						   "problem while getting group info");

	for(hsize_t ii = 0; ii < ginfo.nlinks; ++ii)
	{

		ssize_t buffer_size =
				H5Lget_name_by_idx(grp_id, ".",
								   H5_INDEX_NAME,
								   H5_ITER_INC,
								   ii, 0, 0, H5P_DEFAULT)+1;
		char* buffer = new char[buffer_size];
		buffer_size =
				H5Lget_name_by_idx(grp_id, ".",
								   H5_INDEX_NAME,
								   H5_ITER_INC,
								   ii, buffer,
								   (size_t)buffer_size,
								   H5P_DEFAULT);
		cont.insert(cont.end(), std::string(buffer));
		delete [] buffer;
	}
	if(groupname != "")
		H5Gclose(grp_id);
	H5Fclose(file_id);
}
	



	template<class T>
	void rf_import_HDF5(RandomForest<T> &rf, 
					   std::string filename, 
					   std::string pathname = "")
	{
		typedef MultiArrayShape<2>::type Shp; 	
		// get serialized options_ if they are there
		try
		{
			std::string opt_name = pathname + std::string("_options");
		
			HDF5ImportInfo opt_info(filename.c_str(), opt_name.c_str());
    		MultiArray<2 , double> opt(Shp(opt_info.shapeOfDimension(0), 
				 						   opt_info.shapeOfDimension(1)));
			loadFromHDF5File(opt_info, opt);
			rf.options_.unserialize(opt.begin(), opt.end());
		}
		catch(char const* i)
		{}

		//get external parameters
		std::string ext_path = pathname + "_ext_param/";
		std::set<std::string> ext_set;
		find_groups_hdf5(filename, ext_path, ext_set);
		std::map<std::string, MultiArray<2, double> > ext_map;
		std::set<std::string>::iterator iter;
		if(ext_set.find(std::string("labels")) == ext_set.end())
			std::runtime_error("labels are missing");
		for(iter = ext_set.begin(); iter != ext_set.end(); ++ iter)
		{
			if(		*iter != std::string("weights") 
			   && 	*iter != std::string("labels"))
			{
				ext_map[*iter] = MultiArray<2, double>(Shp(1,1));
				std::string name = pathname + "_ext_param/" + *iter;
				HDF5ImportInfo info(filename.c_str(), name.c_str());
				loadFromHDF5File(info, ext_map[*iter]);
			}
		}
		rf.ext_param_.make_from_map(ext_map);
		{
			std::string name = pathname + "_ext_param/weights";
			HDF5ImportInfo info(filename.c_str(), name.c_str());
			MultiArray<2, double> c_w(Shp(info.shapeOfDimension(0),
										  info.shapeOfDimension(1)));
			loadFromHDF5File(info, c_w);
			rf.ext_param_.class_weights(c_w.begin(), c_w.end());
		}
	
		{
			std::string name = pathname + "_ext_param/labels";
			HDF5ImportInfo info(filename.c_str(), name.c_str());
			switch(info.pixelType())
			{
			#define TYPECASE(ENUM_, REAL_)\
				case HDF5ImportInfo::ENUM_ :\
				{\
					MultiArray<2, REAL_> c_(Shp(info.shapeOfDimension(0),\
										  info.shapeOfDimension(1)));\
					loadFromHDF5File(info, c_);\
					rf.ext_param_.classes_(c_.begin(), c_.end());\
					break;\
				}
				TYPECASE(UINT8, UInt8);
				TYPECASE(UINT16, UInt16);
				TYPECASE(UINT32, UInt32);
				TYPECASE(UINT64, UInt64);
				TYPECASE(INT8, Int8);
				TYPECASE(INT16, Int16);
				TYPECASE(INT32, Int32);
				TYPECASE(INT64, Int64);
				TYPECASE(FLOAT, float);
				TYPECASE(DOUBLE, double);
				default:
					std::runtime_error("rf_import_hdf5() : Unknown class type");
			#undef TYPECASE
			}
		}

		// TREE SAVING TIME
		// get all groups in base path
		std::set<std::string> tree_set;
		find_groups_hdf5(filename, pathname, tree_set);
		for(iter = tree_set.begin(); iter != tree_set.end(); ++iter)
		{
			if((*iter)[0] != '_')
			{
				typedef typename RandomForest<T>::DecisionTree_t DecisionTree;
				rf.trees_.push_back(DecisionTree(rf.ext_param_));
				//get topology
				std::string name = pathname + *iter + "/topology";
				HDF5ImportInfo info(filename.c_str(), name.c_str());
				rf.trees_.back().topology_.resize(info.shapeOfDimension(0) *
												 info.shapeOfDimension(1));
				MultiArrayView<2, Int32> topo(Shp(info.shapeOfDimension(0),
										    	  info.shapeOfDimension(1)),
											  rf.trees_.back().topology_.data());
				loadFromHDF5File(info, topo);

				//get parameters;
				std::string pname = pathname + *iter + "/parameters";
				HDF5ImportInfo pinfo(filename.c_str(), pname.c_str());
				rf.trees_.back().parameters_.resize(pinfo.shapeOfDimension(0) *
												 pinfo.shapeOfDimension(1));
				MultiArrayView<2, double> param(Shp(pinfo.shapeOfDimension(0),
										    	  pinfo.shapeOfDimension(1)),
											  rf.trees_.back().parameters_.data());
				loadFromHDF5File(pinfo, param);
			}
		}
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
	template<class T>
	void rf_export_HDF5(RandomForest<T> &rf, 
					   std::string filename, 
					   std::string pathname= "")
	{ 
		typedef MultiArrayShape<2>::type Shp;
		int a = get_number_of_digits(rf.options_.tree_count_-1);
		//save trees
		for(int ii = 0; ii < rf.options_.tree_count_; ++ii)
		{
			std::stringstream sout; 
			int bb = get_number_of_digits(ii); 
			sout << pathname << "Tree_"; 
			for(int gg = 0; gg < a - bb; ++ gg)
				sout << 0;
			sout << ii << "/";
			std::string topo_name_  = sout.str() + std::string("topology");
			std::string param_name_ = sout.str() + std::string("parameters");
			MultiArrayView<2, Int32> topology(Shp(rf.tree(ii).topology_.size(), 1),
											  rf.tree(ii).topology_.data());
			MultiArrayView<2, double> param(Shp(rf.tree(ii).parameters_.size(), 1),
											rf.tree(ii).parameters_.data());
			writeToHDF5File(filename.c_str(), topo_name_.c_str(),  topology);
			writeToHDF5File(filename.c_str(), param_name_.c_str(), param);
		}

		//save serialized options
		MultiArray<2, double> serialized_options(Shp(rf.options_.serialized_size(), 1));
		rf.options_.serialize(serialized_options.begin(),
							  serialized_options.end()); 
		std::string opt_name = pathname + std::string("_options");
		writeToHDF5File(filename.c_str(), opt_name.c_str(), serialized_options);

		//external parameters are save in _ext_param - unlike options they 
		//are not serialized
		
		//save class weights
		if(rf.ext_param_.is_weighted)
		{
			MultiArrayView<2, double> class_weights(Shp(rf.ext_param_.class_weights_.size(), 1), 
													rf.ext_param_.class_weights_.data());
			std::string weights_name = pathname + std::string("_ext_param/weights");
			writeToHDF5File(filename.c_str(), weights_name.c_str(), class_weights);
		}

		//save class Labels
		switch(rf.ext_param_.class_type_)
		{
			#define SOME_CASE(type_) \
 			case ProblemSpec::type_##_t:\
			{\
				MultiArrayView<2, type_> labels(Shp(rf.ext_param_.type_##_classes_.size(), 1), \
												rf.ext_param_.type_##_classes_.data());\
				std::string label_name = pathname + std::string("_ext_param/labels"); \
				writeToHDF5File(filename.c_str(), label_name.c_str(), labels);\
				break;\
			}
			SOME_CASE(UInt8);
			SOME_CASE(UInt16);
			SOME_CASE(UInt32);
			SOME_CASE(UInt64);
			SOME_CASE(Int8);
			SOME_CASE(Int16);
			SOME_CASE(Int32);
			SOME_CASE(Int64);
			SOME_CASE(double);
			SOME_CASE(float);
			case ProblemSpec::UNKNOWN:
				break;
			default:
				std::runtime_error("exportRF_HDF5(): unknown class type"); 
		}

		//save other stuff
		
		std::map<std::string, MultiArray<2,double> > serialized_param;
		rf.ext_param_.make_map(serialized_param);
		std::map<std::string, MultiArray<2, double> >::iterator iter;
		for(iter = serialized_param.begin(); iter != serialized_param.end(); ++iter)
		{
			std::string name = pathname + std::string("_ext_param/") + (*iter).first;
			writeToHDF5File(filename.c_str(), name.c_str(), (*iter).second);
		}
	}


} // namespace vigra


#endif // VIGRA_RANDOM_FOREST_IMPEX_HXX

