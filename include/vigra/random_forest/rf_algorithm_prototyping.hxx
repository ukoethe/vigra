#ifndef VIGRA_RF_ALGORITHM_HXX
#define VIGRA_RF_ALGORTIHM_HXX
/* First idea for algorithms class.
 * delete file once all methods and classes have been ported
 */
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

/** Container for output
 */
typedef std::vector<int> SelFeats_t;


/** Variable Selection Error summary
 */
struct VarSelectOutput
{
	double before;
	double after;
	

	/** Error rate before variable
	 * selection
	 */
	double before()
	{
		return before;
	}


	/** Error rate after variable selection
	 */
	double after()
	{
		return after;
	}
	

	/** before - after
	 */
	double improvement()
	{
		return before - after;
	}
};


/** Classifier typedef.
 *  has been kept as typedef as somebody may want to apply the
 *  same algorithms with different classifiers - made to simplify
 *  porting when necessary
 */
typedef RandomForest<>   ClassifierT;

/** Perform forward selection using Random FOrests
 * 
 * \param features 	Matrix containing the features used.
 * \param label	   	Matrix containing the corresponding labels
 * \param selected_feats 
 * 					- output. Linear container which will contain
 * 					the selected features.
 * \output Selection summary.
 * \param tolerance between best solution and selected solution 
 * with less features (fraction of error with all features)
 *
 * 	FeatureT and LabelT should be vigra::MultiArray compatible 
 * 	SelFeatsT should be a back insertion containerr i.e. std::vector
 * 	or vigra::ArrayVector
 */
template<class FeatureT, class LabelT>
VarSelectOutput forward_select(		FeatureT & features, 
									LabelT & labels, 
									SelFeatsT & selected_feats,
									double tolerance = 0.0)
{
	int featureCount = features.shape(1);
	std::vector<int> selected;
	std::vector<int> not_selected;
	for(int ii = 0; ii < featureCount; ++ii)
	{
		not_selected.push_back(ii);
	}
	while(not_selected.size() != 0)
	{
		std::vector<int> current_errors(not_selected.size(), 1);
		for(int ii = 0; ii < not_selected.size(); ++ii)
		{
			selected.push_back(not_selected[ii]);
			MultiArray<2, double> cur_feats = choose( features, 
													  selected.begin(), 
											    	  selected.end());
			selected.pop_back();
			visitors::OOB_Error oob;
			visitors::RandomForestProgressVisitor progress;
			ClassifierT classifier;
			classifier.learn(cur_feats, 
							 labels, 
							 create_visitor(oob, progress));
			current_errors.push_back(oob.oob_breiman);
		}
		int pos = std::min_element(current_errors.begin(),
								   current_errors.end()) 
				  - 	current_errors.begin();
		selected.push_back(not_selected[pos]);
		errors.push_back(current_errors[pos]);
		not_selected.erase(pos);
	}
}

/** Perform backward elimination using Random Forests
 * 
 * \param features 	Matrix containing the features used.
 * \param label	   	Matrix containing the corresponding labels
 * \param selected_feats 
 * 					- output. Linear container which will contain
 * 					the selected features.
 * \output Selection summary.
 * \param tolerance between best solution and selected solution 
 * with less features (fraction of error with all features)
 *
 * 	FeatureT and LabelT should be vigra::MultiArray compatible 
 * 	SelFeatsT should be a back insertion containerr i.e. std::vector
 * 	or vigra::ArrayVector
 */
template<class FeatureT, class LabelT>
VarSelectOutput backward_eliminate(	FeatureT & features, 
									LabelT & labels, 
									SelFeatsT & selected_feats,
									double tolerance = 0.0);
{
	int featureCount = features.shape(1);
	std::vector<int> selected;
	std::vector<int> not_selected;
	for(int ii = 0; ii < featureCount; ++ii)
	{
		selected.push_back(ii);
	}
	while(selected.size() != 0)
	{
		std::vector<int> current_errors(not_selected.size(), 1);
		for(int ii = 0; ii < not_selected.size(); ++ii)
		{
			selected.push_back(not_selected[ii]);
			MultiArray<2, double> cur_feats = choose( features, 
													  selected.begin(), 
											    	  selected.end());
			selected.pop_back();
			visitors::OOB_Error oob;
			visitors::RandomForestProgressVisitor progress;
			ClassifierT classifier;
			classifier.learn(cur_feats, 
							 labels, 
							 create_visitor(oob, progress));
			current_errors.push_back(oob.oob_breiman);
		}
		int pos = std::min_element(current_errors.begin(),
								   current_errors.end()) 
				  - 	current_errors.begin();
		selected.push_back(not_selected[pos]);
		errors.push_back(current_errors[pos]);
		not_selected.erase(pos);
	}

}


/** Perform rank selection using Random Forests and a fixed predefined
 *  ranking.
 * 
 * \param features 	Matrix containing the features used.
 * \param label	   	Matrix containing the corresponding labels
 * \param ranking	ranking of features by relevance (i.e. rf variable
 * 					importance measures)
 * \param selected_feats 
 * 					- output. Linear container which will contain
 * 					the selected features.
 * \output Selection summary.
 * \param tolerance between best solution and selected solution 
 * with less features (fraction of error with all features)
 *
 * 	FeatureT and LabelT should be vigra::MultiArray compatible 
 * 	SelFeatsT should be a back insertion containerr i.e. std::vector
 * 	or vigra::ArrayVector
 */
template<class FeatureT, class LabelT, class RankingT>
VarSelectOutput rank_select(			FeatureT & features, 
										LabelT & labels, 
										RankingT & ranking,
										SelFeatsT & selected_feats,
										double tolerance = 0.0);
{
	typename RankingT::iterator iter = ranking.begin();
	for(; iter != ranking.end(); ++iter)
	{
		MultiArray<2, double> cur_feats = choose( features, 
												  ranking.begin(), 
										    	  iter);
		ClassifierT classifier;
		classifier.learn(cur_feats, 
						 labels, 
						 create_visitor(oob, progress));
		errors.push_back(oob.oob_breiman);

	}

}

}//namespace algorithms
}//namespace rf
}//namespace vigra

#undef //VIGRA_RF_ALGORITHM_HXX
