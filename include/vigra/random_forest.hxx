/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by  Ullrich Koethe and Rahul Nair         */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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


#ifndef VIGRA_RANDOM_FOREST_HXX
#define VIGRA_RANDOM_FOREST_HXX

#include <algorithm>
#include <map>
#include <set>
#include <numeric>
#include "mathutil.hxx"
#include "array_vector.hxx"
#include "sized_int.hxx"
#include "matrix.hxx"
#include "random.hxx"
#include "functorexpression.hxx"
#include "random_forest/rf_common.hxx"
#include "random_forest/rf_nodeproxy.hxx"
#include "random_forest/rf_split.hxx"
#include "random_forest/rf_decisionTree.hxx"
#include "random_forest/rf_visitors.hxx"
#include "random_forest/rf_region.hxx"
#include "random_forest/rf_sampling.hxx"
#include "random_forest/rf_preprocessing.hxx"
namespace vigra
{

namespace detail
{

/** todo - remove and make the labels parameter in the sampling options
 * const*/
class staticMultiArrayViewHelper
{
    public:
        static vigra::MultiArrayView<2, Int32> array;
    public:
        friend SamplingOptions
        createSamplingOptions(vigra::RandomForestOptions& RF_opt,
							  vigra::MultiArrayView<2, int> & labels);
};


/**\brief sampling option factory function
 */
SamplingOptions make_sampler_opt ( RF_Traits::Options_t		& RF_opt,
                                   MultiArrayView<2, Int32> & labels
										= staticMultiArrayViewHelper::array)
{
    SamplingOptions return_opt;
    return_opt.sample_with_replacement = RF_opt.sample_with_replacement_;
    if(labels.data() != 0)
    {
        if(RF_opt.stratification_method_ == RF_EQUAL)
            return_opt
				.sampleClassesIndividually(
					ArrayVectorView<int>(labels.size(),
										 labels.data()));
        else if(RF_opt.stratification_method_ == RF_PROPORTIONAL)
            return_opt
				.sampleStratified(
					ArrayVectorView<int>(labels.size(),
										 labels.data()));
    }
    return return_opt;
}
}//namespace detail

/** Random Forest class
 *
 * \tparam <PrprocessorTag = ClassificationTag> Class used to preprocess
 * 			the input while learning and predicting. Currently Available:
 * 			ClassificationTag and RegressionTag. It is recommended to use
 * 			Splitfunctor::Preprocessor_t while using custom splitfunctors
 * 			as they may need the data to be in a different format. 
 * 			\sa Preprocessor, How to make a Split Functor
 *
 * 
 *
*/
template <class PreprocessorTag = ClassificationTag>
class RandomForest
{

  public:
	//public typedefs
	typedef RF_Traits::Options_t 			Options_t;
	typedef	RF_Traits::DecisionTree_t 		DecisionTree_t;
	typedef RF_Traits::ProblemSpec_t 		ProblemSpec_t;
	typedef RF_Traits::Default_Split_t 		Default_Split_t;
	typedef	RF_Traits::Default_Stop_t 		Default_Stop_t;
	typedef	RF_Traits::Default_Visitor_t 	Default_Visitor_t;
  protected:
	/**\brief helper function used to choose the right
	 * value of split, early stopping and visitor
	 */



  public:

	//problem independent data.
	Options_t			 			options_;
	//problem dependent data members - is only set if
	//a copy constructor, some sort of import
	//function or the learn function is called
    ArrayVector<DecisionTree_t>
	   								trees_;
	ProblemSpec_t 					ext_param_;


    void reset()
    {
        ext_param_.clear();
        trees_.clear();
    }

  public:

	/** \name Contructors
	 * Note: No copy Constructor specified as no pointers are manipulated
	 * in this class
	 */
	/*\{*/
	/**\brief default constructor
	 *
	 * \param options 	general options to the Random Forest. Must be of Type
	 * 				  	Options_t
	 * \param ext_param problem specific values that can be supplied 
	 * 					additionally. (class weights , labels etc)
	 * \sa 	ProblemSpec_t
	 *
	 *
	 * simple usage for classification (regression is not yet supported):
	 * \code
	 * 		typedef xxx	feature_t \\ replace xxx wit whichever type
	 * 		typedef yyy label_t	  \\ meme chose. 
	 * 		MultiArrayView<2, feature_t> f = get_some_features();
	 * 		MultiArrayView<2, label_t>   l = get_some_labels)(
	 * 		RandomForest<> rf()
	 * 		double oob_error = rf.learn(f, l);
	 *		
	 *		MultiArrayView<2, feature_t> pf = get_some_unknown_features();
	 *		MultiArrayView<2, label_t> prediction 
	 *										= allocate_space_for_response();
	 *		MultiArrayView<2, double> prob  = allocate_space_for_probability();
	 *		
	 *		rf.predict_labels(pf, prediction);
	 *		rf.predict_probabilities(pf, prob);
	 *
	 * \endcode
	 *
	 * Classification Random Forest specific:
	 * The Random Forest remembers the LabelType (without specifying it as a 
	 * template parameter during construction):
	 * - Only numeral and native C types (IntX, UIntX, float, double)
	 * 	 are suppported as label type.
	 * - Types are automatically converted to the right type give in the 
	 *   predict function - This may cause problems if the original type
	 *   was double - and you were stu... unvigilant enough to name the 
	 *   labels	3.4 , 3.2 etc.
	 *   This conversion is only done once during construction and causes
	 *   no type casting overhead during prediction.
	 */
    RandomForest(Options_t const & options = Options_t(), 
				 ProblemSpec_t const & ext_param = ProblemSpec_t())
    :
		options_(options),
		ext_param_(ext_param)
    {}

	/**\brief Create RF from external source
	 *
	 * \param ext_param Extrinsic parameters that specify the problem e.g.
	 * 					ClassCount, featureCount etc.
	 * \param tree_top	Iterator to a Container where the topology_ data
	 * 					of the trees are stored.
	 * \param tree_par  iterator to a Container where the parameters_ data
	 * 					of the trees are stored.
	 * \param options  	(optional) specify options used to train the original
	 * 					Random forest. This parameter is not used anywhere
	 * 					during prediction and thus is optional.
	 *
	 * TODO:
	 * Note: This constructor may be replaced by a Constructor using
	 * NodeProxy iterators to encapsulate the underlying data type.
	 */
    template<class TreeIterator, class WeightIterator>
    RandomForest(int 				treeCount,
                  TreeIterator 			trees,
				  WeightIterator 		weights,
				  ProblemSpec_t const & problem_spec,
				  Options_t const & 	options = Options_t())
    :
		trees_(treeCount, DecisionTree_t(problem_spec)),
      	ext_param_(problem_spec),
		options_(options)
    {
        for(unsigned int k=0; k<treeCount; ++k, ++trees, ++weights)
        {
            trees_[k].topology_ = *trees;
            trees_[k].parameters_ = *weights;
        }
    }

	/*\}*/


	/** \name Data Access
	 * data access interface - usage of member objects is deprecated
	 * (I like the word deprecated)
	 */

	/*\{*/


	/**\brief return external parameters for viewing
	 * \return ProblemSpec_t
	 */
	ProblemSpec_t const & ext_param() const
    {
        vigra_precondition(ext_param_.used() == true,
           "RandomForest::ext_param(): "
		   "Random forest has not been trained yet.");
        return ext_param_;
    }

	/**\brief set external parameters
	 *
	 *  \param in external parameters to be set
	 *
	 * set external parameters explicitly if Random Forest has not been
	 * trained the preprocessor will either ignore filling values set this
	 * way or will throw an exception if values specified manually do not
	 * match the value calculated during the preparation step.
	 * \sa Option_t::presupplied_ext_param member for further details.
	 */
	void set_ext_param(ProblemSpec_t const & in)
	{
		vigra_precondition(ext_param_.used() == false,
			"RandomForest::set_ext_param():"
			"Random forest has been trained! Call reset()"
			"before specifying new extrinsic parameters.");
	}

	/**\brief access random forest options
	 *
	 * \return random forest options
	 */
	Options_t & set_options()
	{
		return options;
	}


	/**\brief access const random forest options
	 *
	 * \return const Option_t
	 */
	Options_t const & options() const
	{
		return options_;
	}

	/**\brief access const trees
	 */
	DecisionTree_t const & tree(int index) const
	{
		return trees_[index];
	}

	/**\brief access trees
	 */
	DecisionTree_t & tree(int index)
	{
		return trees_[index];
	}

	/*\}*/


	/**\name Learning
	 * Following functions differ in the degree of customization
	 * allowed
	 */
	/*\{*/
	/**\brief learn on data with custom config and random number generator
     *
	 * \param features 	a N x M matrix containing N samples with M
	 * 					features
	 * \param response 	a N x D matrix containing the corresponding
	 * 					response. Current split functors assume D to
	 * 					be 1 and ignore any additional columns.
	 * 					This is not enforced to allow future support
	 * 					for uncertain labels, label independent strata etc.
	 * 					The Preprocessor specified during construction
	 * 					should be able to handle features and labels
	 * 					features and the labels.
	 * 	\sa    	SplitFunctor, Preprocessing
	 *
	 * \param visitor 	visitor which is to be applied after each split,
	 * 					tree and at the end. Use RF_Default for using
	 * 					default value.
	 * \sa 		visitor
	 * \param split		split functor to be used to calculate each split
	 * 					use RF_Default() for using default value.
	 * \param stop
	 * 					predicate to be used to calculate each split
	 * 					use RF_Default() for using default value.
	 * \param random	RandomNumberGenerator to be used. Use
	 * 					RF_Default() to use default value.
	 * \return			oob_error.
	 *
	 *\sa OOB_Visitor, VariableImportanceVisitor 
	 *
	 */
	template <class U, class C1,
			 class U2,class C2,
			 class Split_t,
			 class Stop_t,
			 class Visitor_t,
			 class Random_t>
	double learn(		MultiArrayView<2, U, C1> const  & 	features,
						MultiArrayView<2, U2,C2> const  & 	response,
						Visitor_t 							visitor,
						Split_t 							split,
						Stop_t 								stop,
						Random_t 				 const 	&	random);

	template <class U, class C1,
			 class U2,class C2,
			 class Split_t,
			 class Stop_t,
			 class Visitor_t>
	double learn(		MultiArrayView<2, U, C1> const  & 	features,
						MultiArrayView<2, U2,C2> const  & 	response,
						Visitor_t 							visitor,
						Split_t 							split,
						Stop_t 								stop)

	{
        RandomNumberGenerator<> rnd = RandomNumberGenerator<>(RandomSeed);
		return learn(features, response,
					 visitor, split, stop,
					 rnd);
	}



	/**\brief learn on data with default configuration
	 *
	 * \param features 	a N x M matrix containing N samples with M
	 * 					features
	 * \param labels 	a N x D matrix containing the corresponding
	 * 					N labels. Current split functors assume D to
	 * 					be 1 and ignore any additional columns.
	 * 					this is not enforced to allow future support
	 * 					for uncertain labels.
	 * \return 			out of bag error estimate.
	 *
	 * learning is done with:
	 *
	 * \sa GiniSplit, EarlyStoppingStd, OOB_Visitor
	 *
	 * - Randomly seeded random number generator
	 * - default gini split functor as described by Breiman
	 * - default The standard early stopping criterion
	 * - the oob visitor, whose value is returned.
	 */
    template <class U, class C1, class U2,class C2, class Visitor_t>
    double learn(   MultiArrayView<2, U, C1> const  & features,
                    MultiArrayView<2, U2,C2> const  & labels,
					Visitor_t 						  visitor)
    {
        return learn(features, 
					 labels, 
					 visitor, 
					 rf_default(), 
					 rf_default());
    }

	template <class U, class C1, class U2,class C2, 
			  class Visitor_t, class Split_t>
    double learn(   MultiArrayView<2, U, C1> const  & features,
                    MultiArrayView<2, U2,C2> const  & labels,
					Visitor_t						  visitor,
					Split_t							  split)
    {
        return learn(features, 
					 labels, 
					 visitor, 
					 split, 
					 rf_default());
    }

	template <class U, class C1, class U2,class C2>
    double learn(   MultiArrayView<2, U, C1> const  & features,
                    MultiArrayView<2, U2,C2> const  & labels)
    {
        return learn(features, 
					 labels, 
					 rf_default(), 
					 rf_default(), 
					 rf_default());
    }
	/*\}*/



	/**\name prediction
	 */
	/*\{*/
	/** \brief predict a label given a feature.
	 *
	 * \param features: a 1 by featureCount matrix containing
	 *        data point to be predicted (this only works in
	 *        classification setting)
	 * \return double value representing class. You can use the
	 *         predictLabels() function together with the
	 *         rf.external_parameter().class_type_ attribute
	 *         to get back the same type used during learning. 
	 */
    template <class U, class C>
    double predictLabel(MultiArrayView<2, U, C>const & features);

	/** \brief predict a label with features and class priors
	 *
	 * \param features: same as above.
	 * \param priors:   iterator to prior weighting of classes
	 * \return sam as above.
	 */
    template <class U, class C, class Iterator>
    double predictLabel(MultiArrayView<2, U, C> const & features,
                                Iterator priors) ;

	/** \brief predict multiple labels with given features
	 *
	 * \param features: a n by featureCount matrix containing
	 *        data point to be predicted (this only works in
	 *        classification setting)
	 * \param labels: a n by 1 matrix passed by reference to store
     *        output.
	 */
    template <class U, class C1, class T, class C2>
    void predictLabels(MultiArrayView<2, U, C1>const & features,
                       MultiArrayView<2, T, C2> & labels)
    {
        vigra_precondition(features.shape(0) == labels.shape(0),
            "RandomForest::predictLabels(): Label array has wrong size.");
        for(int k=0; k<features.shape(0); ++k)
            labels(k,0) = predictLabel(rowVector(features, k));
    }


    /** \brief predict the class probabilities for multiple labels
     *
     *  \param features same as above
     *  \param prob a n x class_count_ matrix. passed by reference to
     *  save class probabilities
     */
    template <class U, class C1, class T, class C2>
    void predictProbabilities(MultiArrayView<2, U, C1>const & 	features,
                              MultiArrayView<2, T, C2> & 		prob) ;



	/*\}*/

};


template <class PreprocessorTag>
template <class U, class C1,
		 class U2,class C2,
		 class Split_t,
		 class Stop_t,
		 class Visitor_t,
		 class Random_t>
double RandomForest<PreprocessorTag>::
					 learn( MultiArrayView<2, U, C1> const  & 	features,
							MultiArrayView<2, U2,C2> const  & 	response,
							Visitor_t 							visitor_,
							Split_t 					    	split_,
							Stop_t 								stop_,
							Random_t 				 const 	&	random)
{
	using namespace rf;
	this->reset();
	//typedefs
	typedef typename Split_t::StackEntry_t 			StackEntry_t;
    typedef 		 UniformIntRandomFunctor<Random_t>
													RandFunctor_t;

    // See rf_preprocessing.hxx for more info on this
    typedef Processor<PreprocessorTag, U, C1, U2, C2> Preprocessor_t;
	
	// default values and initialization
	// Value Chooser chooses second argument as value if first argument
	// is of type RF_DEFAULT. (thanks to template magic - don't care about
	// it - just smile and wave.
	
	#define RF_CHOOSER(type_) detail::Value_Chooser<type_, Default_##type_>	
	Default_Stop_t default_stop(options_);
	typename RF_CHOOSER(Stop_t)::type stop
			= RF_CHOOSER(Stop_t)::choose(stop_, default_stop); 
	Default_Split_t default_split;
	typename RF_CHOOSER(Split_t)::type split 
			= RF_CHOOSER(Split_t)::choose(split_, default_split); 
	StopVisiting stopvisiting;
	OOB_Visitor	 oob;
	VisitorNode<OOB_Visitor, typename RF_CHOOSER(Visitor_t)::type>
		visitor(oob, RF_CHOOSER(Visitor_t)::choose(visitor_, stopvisiting));
	#undef RF_CHOOSER

    // Make stl compatible random functor.
    RandFunctor_t			randint		( random);


	// Preprocess the data to get something the split functor can work
	// with. Also fill the ext_param structure by preprocessing
	// option parameters that could only be completely evaluated
	// when the training data is known.
	Preprocessor_t preprocessor( 	features, response,
                                    options_, ext_param_);

    // Give the Split functor information about the data.
    split.set_external_parameters(ext_param_);
    stop.set_external_parameters(ext_param_);


	//initialize trees.
	trees_.resize(options_.tree_count_	, DecisionTree_t(ext_param_));

	/**\todo 	replace this crappy class out. It uses function pointers.
	 * 			and is making code slower according to me
	 */
    Sampler<RandFunctor_t >	sampler(ext_param().row_count_,
									ext_param().actual_msample_,
								   	detail::make_sampler_opt(options_,
                                                     preprocessor.strata()),
								    randint);
	visitor.visit_at_beginning(*this);
	// THE MAIN EFFING RF LOOP - YEAY DUDE!
    for(int ii = 0; ii < (int)trees_.size(); ++ii)
    {
		//initialize First region/node/stack entry
        sampler
			.sample();

		StackEntry_t
			first_stack_entry( 	sampler.used_indices().begin(),
					   			sampler.used_indices().end(),
					   			ext_param_.class_count_);
        first_stack_entry
			.set_oob_range(		sampler.unused_indices().begin(),
								sampler.unused_indices().end());
		trees_[ii]
			.learn(				preprocessor.features(),
					   			preprocessor.response(),
                        		first_stack_entry,
                				split,
								stop,
								visitor,
								randint);
		visitor
			.visit_after_tree(	*this,
                                preprocessor,
								sampler,
								first_stack_entry,
								ii);
	}

	visitor.visit_at_end(*this, preprocessor);

	return 	visitor.return_val();
}







template <class Tag>
template <class U, class C>
double RandomForest<Tag>
    ::predictLabel(MultiArrayView<2, U, C> const & features)
{
    vigra_precondition(columnCount(features) >= ext_param_.column_count_,
        "RandomForestn::predictLabel():"
            " Too few columns in feature matrix.");
    vigra_precondition(rowCount(features) == 1,
        "RandomForestn::predictLabel():"
            " Feature matrix must have a singlerow.");
    Matrix<double>  prob(1, ext_param_.class_count_);
    double          d;
    predictProbabilities(features, prob);
    ext_param_.to_classlabel(argMax(prob), d);
    return d;
}


//Same thing as above with priors for each label !!!
template <class PreprocessorTag>
template <class U, class C, class Iterator>
double RandomForest<PreprocessorTag>
    ::predictLabel( MultiArrayView<2, U, C> const & features,
                    Iterator                        priors)
{
    using namespace functor;
    vigra_precondition(columnCount(features) >= ext_param_.column_count_,
        "RandomForestn::predictLabel(): Too few columns in feature matrix.");
    vigra_precondition(rowCount(features) == 1,
        "RandomForestn::predictLabel():"
        " Feature matrix must have a single row.");
    Matrix<double>  prob(1,ext_param_.class_count_);
    predictProbabilities(features, prob);
    std::transform( prob.begin(), prob.end(),
                    priors, prob.begin(),
                    Arg1()*Arg2());
    double          d;
    ext_param_.to_classlabel(argMax(prob), d);
    return d;
}

template <class PreprocessorTag>
template <class U, class C1, class T, class C2>
void RandomForest<PreprocessorTag>
    ::predictProbabilities(MultiArrayView<2, U, C1>const &  features,
                           MultiArrayView<2, T, C2> &       prob)
{

	//Features are n xp
	//prob is n x NumOfLabel probability for each feature in each class

    vigra_precondition(rowCount(features) == rowCount(prob),
      "RandomForestn::predictProbabilities():"
        " Feature matrix and probability matrix size mismatch.");

	// num of features must be bigger than num of features in Random forest training
	// but why bigger?
    vigra_precondition( columnCount(features) >= ext_param_.column_count_,
      "RandomForestn::predictProbabilities():"
        " Too few columns in feature matrix.");
    vigra_precondition( columnCount(prob)
                        == (MultiArrayIndex)ext_param_.class_count_,
      "RandomForestn::predictProbabilities():"
      " Probability matrix must have as many columns as there are classes.");

	//Classify for each row.
    for(int row=0; row < rowCount(features); ++row)
    {
	//contains the weights returned by a single tree???
	//thought that one tree has only one vote???
	//Pruning???
        ArrayVector<double>::const_iterator weights;

        //totalWeight == totalVoteCount!
	double totalWeight = 0.0;

	//Set each VoteCount = 0 - prob(row,l) contains vote counts until
	//further normalisation
        for(unsigned int l=0; l<ext_param_.class_count_; ++l)
            prob(row, l) = 0.0;

	//Let each tree classify...
        for(unsigned int k=0; k<options_.tree_count_; ++k)
        {
		//get weights predicted by single tree
            weights = trees_[k].predict(rowVector(features, row));

		//update votecount.
            for(unsigned int l=0; l<ext_param_.class_count_; ++l)
            {
                prob(row, l) += weights[l];
                //every weight in totalWeight.
		totalWeight += weights[l];
            }
        }

	//Normalise votes in each row by total VoteCount (totalWeight
        for(unsigned int l=0; l<ext_param_.class_count_; ++l)
                prob(row, l) /= totalWeight;
    }

}


} // namespace vigra

#endif // VIGRA_RANDOM_FOREST_HXX
