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

	template<class T>
	void importRF_HDF5(RandomForest<T> &rf, 
					   std::string filename, 
					   std::string pathname = "")
	{
		typedef MultiArrayShape<2, T> Shp; 	
		// get serialized options_ if they are there
		
		// get ext_param
		std::string name = pathname + std::string("_options");
		HDF5ImportInfo info(filename, name); 

		// get all folders , if they start with _ ignore, else make trees. 
	}

	template<class T>
	void exportRF_HDF5(RandomForest<T> &rf, 
					   std::string filename, 
					   std::string pathname= "")
	{ 
		typedef MultiArrayShape<2, T> Shp;
		//save trees
		for(int ii = 0; ii < rf.options_.tree_count_; ++ii)
		{
			std::stringstream sout; 
			sout << pathname << "Tree_" << ii;
			std::string topo_name_  = sout.str() + std::string("topology");
			std::string param_name_ = sout.str() + std::string("parameters");
			MultiArrayView<2, Int32> topology(Shp(rf.tree(ii).topology_.size(), 1),
											  rf.tree(ii).topology_.data());
			MultiArrayView<2, double> param(Shp(rf.tree(ii).parameters_size(), 1),
											rf.tree(ii).parameters_.data());
			writeToHDF5File(filename.c_str(), topo_name_.c_str(),  topology);
			writeToHDF5File(filename.c_str(), param_name_.c_str(), param);
		}

		//save serialized options
		MultiArray<2, double> serialized_options_(Shp(rf.options_.serialized_size(), 1));
		rf.options_.serialize(serialized_options_.begin(),
							  serialized_options_.end()); 
		std::string opt_name = pathname + std::string("_options");
		writeToHDF5File(filename.c_str(), opt_name.c_str(), serialized_options);

		//external parameters are save in _ext_param - unlike options they 
		//are not serialized
		
		//save class weights
		if(rf.ext_param_.is_weighted)
			MultiArrayView<2, double> class_weights(Shp(rf.ext_param_.class_weights_.size(), 1), 
													rf.ext_param_.class_weights_.data());
		std::string weights_name = pathname + std::string("/_ext_param/weights");
		writeToHDF5File(filename.c_str(), weights_name, class_weights);

		//save class Labels
		switch(class_type_)
		{
			#define SOME_CASE(type_) \
			class type_##_t:\
			{\
				MultiArrayView<2, type_> labels(Shp(type_##_classes_.size(), 1), \
												type_##_classes_);\
				std::string label_name = pathname + std::string("/_ext_param/labels"); \
				writeToHDF5File(filename.c_str(), label_name, labels);\
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
			class UNKNOWN:
				break;
			default:
				VIGRA_ERROR("exportRF_HDF5(): unknown class type"); 
		}

		//save other stuff
		
		std::map<string, double> serialized_param;
		rf.ext_param_.make_map(serialized_param);
		std::map<string, double>::iterator iter;
		for(iter = serialized_param.begin(); iter != serialized_param.end(); ++iter)
		{
			std::string name = pathname + std::string("/_ext_param") + iter->first;
			writeToHDF5File(filename.c_str(), name, iter->second);
		}
	}


} // namespace vigra


#endif // VIGRA_RANDOM_FOREST_IMPEX_HXX

