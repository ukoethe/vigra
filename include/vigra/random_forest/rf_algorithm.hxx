/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by Rahul Nair                             */
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
#define VIGRA_RF_ALGORTIHM_HXX

#include <vector>

namespace vigra
{
namespace rf
{
/** This namespace contains all algorithms developed for feature 
 * selection
 *
 */
namespace algorithms
{

namespace detail
{
    /** create a MultiArray containing only columns supplied between iterators
        b and e
    */
	template<class OrigMultiArray,
			 class Iter,
			 class DestMultiArray>
	void choose(OrigMultiArray 	const & in,
				Iter		   	const & b,
				Iter   			const & e,
				DestMultiArray  	  & out)
	{
		int columnCount = std::distance(b, e);
		int rowCount 	= in.shape(0);
		out.resize(MultiArrayShape<2>::type(rowCount, columnCount));
		for(Iter iter = b, int ii = 0; iter != e; ++iter, ++ii)
		{
			columnVector(out, ii) = columnVector(in, *iter);
		}
	}
}
/** Structure to hold Variable Selection results
 */
class VariableSelectionResult
{
	bool initialized;

  public:
	VariableSelectionResult()
	: initialized(false)
	{}

	typedef std::vector<int> FeatureList_t

	/** list of features. 
	 */
	FeatureList_t selected;
	
	/** vector of size (number of features)
	 *
	 * the i-th entry encodes the error rate obtained
	 * while using features [0 - i](including i) 
	 *
	 * if the i-th entry is -1 then no error rate was obtained
	 * this may happen if more than one feature is added to the
	 * selected list in one step of the algorithm.
	 *
	 * during initialisation error[m+n-1] is always filled
	 */
	std::vector<double> error;
	

	/** errorrate using no features
	 */
	double no_features;

	/**initialization routine. Will be called only once in the lifetime
	 * of a VariableSelectionResult. Subsequent calls will not reinitialize
	 * member variables.
	 *
	 * This is intended, to allow continuing variable selection at a point 
	 * stopped in an earlier iteration. 
	 *
	 * returns true if initialization was successful and false if 
	 * the object was already initialized before.
	 */
	template<class FeatureT, 
			 class ResponseT, 
			 class ErrorRateCallBack = RFErrorCallback>
	bool init(FeatureT const & all_features,
			  ResponseT const & response,
			  ErrorRateCallBack errorcallback = RFErrorCallback())
	{
		if(initialized)
		{
			return false;
		}
		// calculate error with all features
		error.resize(all_features.shape(1), -1);
		error.back() = errorcallback(all_features, response);

		// calculate error rate if no features are chosen 
		// corresponds to max(prior probability) of the classes
		std::map<typename ResponseT::value_type, int> 	res_map;
		std::vector<int> 								cts;
		int 											counter = 0;
		for(int ii = 0; ii < response.shape(0); ++ii)
		{
			if(res_map.find(response(ii, 0)) == res_map.end())
			{
				res_map[response(ii, 0)] = counter;
				++counter;
			}
			cnts[res_map[response(ii,0)]] +=1;
		}
		no_features = double(*(std::max_element(cts.begin(),
												 cts.end())))
					/ double(response.shape(0));

		//init not_selected vector;
		not_selected.resize(all_features.shape(1), 0);
		for(int ii = 0; ii < not_selected.size(); ++ii)
		{
			not_selected[ii] = ii;
		}
		initialized = true;
		return true;
	}
};

/** Standard random forest Errorrate callback functor
 *
 * returns the random forest error estimate when invoked. 
 */
class RFErrorCallback
{
	RandomForestOptions options;
	
	public:
	/** Default constructor
	 *
	 * optionally supply options to the random forest classifier
	 * \sa RandomForestOptions
	 */
	RFErrorCallback(RandomForestOptions opt = RandomForestOptions())
	: options(opt)
	{}

	/** returns the RF OOB error estimate given features and 
	 * labels
	 */
	template<class Feature_t, class Response_t>
	double operator() (Feature_t const & features,
					   Response_t const & response)
	{
		RandomForest<> 			rf(options);
		visitors::OOB_Error		oob;
		rf.learn(features, 
				 response, 
				 visitors::create_visitor(oob));
		return oob.oob_breiman;
	}
};

	
/** Perform forward selection
 *
 * \param features	IN: 	n x p matrix containing n instances with p attributes/features
 * 							used in the variable selection algorithm
 * \param response  IN: 	n x 1 matrix containing the corresponding response
 * \param result	IN/OUT: VariableSelectionResult struct which will contain the results
 * 							of the algorithm. 
 * 							Features between result.selected.begin() and result.pivot will
 * 							be left untouched.
 * 							\sa VariableSelectionResult
 * \param errorcallback
 * 					IN, OPTIONAL: 
 * 							Functor that returns the error rate given a set of 
 * 							features and labels. Default is the RandomForest OOB Error.
 *
 * Forward selection subsequently chooses the next feature that decreases the Error rate most.
 *
 * usage:
 * \code
 * 		MultiArray<2, double> 	features = createSomeFeatures();
 * 		MultiArray<2, int>		labels   = createCorrespondingLabels();
 * 		VariableSelectionResult  result;
 * 		forward_selection(features, labels, result);
 * \endcode
 * To use forward selection but ensure that a specific feature e.g. feature 5 is always 
 * included one would do the following
 *
 * \code
 * 		VariableSelectionResult result;
 * 		result.init(features, labels);
 * 		std::swap(result.selected[0], result.selected[5]);
 * 		result.setPivot(1);
 * 		forward_selection(features, labels, result);
 * \endcode
 *
 * \sa VariableSelectionResult
 *
 */					
template<class FeatureT, class ResponseT, class ErrorRateCallBack = RFErrorCallback>
void forward_selection(FeatureT 		 const & features,
					   ResponseT 		 const & response,
					   VariableSelectionResult & result,
					   ErrorRateCallBack 		 errorcallback
					   								= RFErrorCallback())
{
	VariableSelectionResult::FeatureList_t & selected 		= result.selected;
	VariableSelectionResult::Pivot_t	   & pivot			= result.pivot;	
	int featureCount = features.shape(1);
	// initialize result struct if in use for the first time
	if(!result.init(features, response, errorcallback))
	{
		//result is being reused just ensure that the number of features is
		//the same.
		vigra_precondition(selected.size() == featureCount,
						   "forward_selection(): Number of features in Feature "
						   "matrix and number of features in previously used "
						   "result struct mismatch!");
	}
	

	int not_selected_size = std::distance(pivot, selected.end());
	while(not_selected_size > 0)
	{
		std::vector<int> current_errors;
		VariableSeletionResult::Pivot_t next = pivot;
		for(int ii = 0; ii < not_selected_size; ++ii, ++next)
		{
			std::swap(*pivot, *next);
			MultiArray<2, double> cur_feats;
			detail::choose( features, 
							selected.begin(), 
							pivot, 
							cur_feats);
			double error = errorcallback(cur_feats, labels);
			current_errors.push_back(error);
			std::swap(*pivot, *next);
		}
		int pos = std::distance(current_errors.begin(),
								std::min_element(current_errors.begin(),
								  				 current_errors.end()));
		next = pivot;
		std::advance(next, pos);
		std::swap(*pivot, *next);
		errors[std::distance(selected.begin, pivot)] = current_errors[pos];
		not_selected_size = std::distance(pivot, selected.end());
		++pivot;
	}
}


/** Perform backward elimination
 *
 * \param features	IN: 	n x p matrix containing n instances with p attributes/features
 * 							used in the variable selection algorithm
 * \param response  IN: 	n x 1 matrix containing the corresponding response
 * \param result	IN/OUT: VariableSelectionResult struct which will contain the results
 * 							of the algorithm. 
 * 							Features between result.pivot and result.selected.end() will
 * 							be left untouched.
 * 							\sa VariableSelectionResult
 * \param errorcallback
 * 					IN, OPTIONAL: 
 * 							Functor that returns the error rate given a set of 
 * 							features and labels. Default is the RandomForest OOB Error.
 *
 * Backward elimination subsequently eliminates features that have the least influence
 * on the error rate
 *
 * usage:
 * \code
 * 		MultiArray<2, double> 	features = createSomeFeatures();
 * 		MultiArray<2, int>		labels   = createCorrespondingLabels();
 * 		VariableSelectionResult  result;
 * 		backward_elimination(features, labels, result);
 * \endcode
 * To use backward elimination but ensure that a specific feature e.g. feature 5 is always 
 * excluded one would do the following:
 *
 * \code
 * 		VariableSelectionResult result;
 * 		result.init(features, labels);
 * 		std::swap(result.selected[result.selected.size()-1], result.selected[5]);
 * 		result.setPivot(result.selected.size()-1);
 * 		backward_elimination(features, labels, result);
 * \endcode
 *
 * \sa VariableSelectionResult
 *
 */					
template<class FeatureT, class ResponseT, class ErrorRateCallBack = RFErrorCallback>
void backward_elimination(FeatureT 		 	const & features,
					   	  ResponseT 		const & response,
					      VariableSelectionResult & result,
					      ErrorRateCallBack 		errorcallback
					   								= RFErrorCallback())
{
	int featureCount = features.shape(1);
	VariableSelectionResult::FeatureList_t & selected 		= result.selected;
	VariableSelectionResult::Pivot_t	   & pivot			= result.pivot;	
	
	
	// initialize result struct if in use for the first time
	if(!result.init(features, response, errorcallback))
	{
		//result is being reused just ensure that the number of features is
		//the same.
		vigra_precondition(selected.size()+not_selected.size() == featureCount,
						   "backward_elimination(): Number of features in Feature "
						   "matrix and number of features in previously used "
						   "result struct mismatch!");
	}

	int selected_size = std::distance(selected.begin(), pivot);
	while(selected_size > 0)
	{
		VariableSeletionResult::Pivot_t next = selected.begin();
		std::vector<int> current_errors;
		for(int ii = 0; ii < selected_size; ++ii, ++next)
		{
			std::swap(*pivot, *next);
			MultiArray<2, double> cur_feats;
			detail::choose( features, 
							selected.begin(), 
							pivot, 
							cur_feats);
			double error = errorcallback(cur_feats, labels);
			current_errors.push_back(error);
			std::swap(*pivot, *next);
		}
		int pos = std::distance(current_errors.begin(),
								std::max_element(current_errors.begin(),
								  				 current_errors.end()));
		next = selected.begin();
		std::advance(next, pos);
		std::swap(*pivot, *next);
		errors[std::distance(selected.begin(), pivot)] = current_errors[pos];
		selected_size = std::distance(selected.begin(), pivot);
		--pivot;
	}
}


}//namespace algorithms
}//namespace rf
}//namespace vigra
