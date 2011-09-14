/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by  Ullrich Koethe and Rahul Nair         */
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


#ifndef VIGRA_RANDOM_FOREST_HXX
#define VIGRA_RANDOM_FOREST_HXX

#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <list>
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
#include "sampling.hxx"
#include "random_forest/rf_preprocessing.hxx"
#include "random_forest/rf_online_prediction_set.hxx"
#include "random_forest/rf_earlystopping.hxx"
#include "random_forest/rf_ridge_split.hxx"
namespace vigra
{

/** \addtogroup MachineLearning Machine Learning

    This module provides classification algorithms that map 
    features to labels or label probabilities.
    Look at the RandomForest class first for a overview of most of the 
    functionality provided as well as use cases. 
**/
//@{

namespace detail
{



/* \brief sampling option factory function
 */
inline SamplerOptions make_sampler_opt ( RandomForestOptions     & RF_opt)
{
    SamplerOptions return_opt;
    return_opt.withReplacement(RF_opt.sample_with_replacement_);
    return_opt.stratified(RF_opt.stratification_method_ == RF_EQUAL);
    return return_opt;
}
}//namespace detail

/** Random Forest class
 *
 * \tparam <PrprocessorTag = ClassificationTag> Class used to preprocess
 *          the input while learning and predicting. Currently Available:
 *          ClassificationTag and RegressionTag. It is recommended to use
 *          Splitfunctor::Preprocessor_t while using custom splitfunctors
 *          as they may need the data to be in a different format. 
 *          \sa Preprocessor
 * 	
 * 	simple usage for classification (regression is not yet supported):
 * 	look at RandomForest::learn() as well as RandomForestOptions() for additional
 * 	options. 
 *
 * 	\code
 *  using namespace vigra;
 *  using namespace rf;
 *  typedef xxx feature_t; \\ replace xxx with whichever type
 *  typedef yyy label_t;   \\ likewise 
 *  
 *  // allocate the training data
 *  MultiArrayView<2, feature_t> f = get_training_features();
 *  MultiArrayView<2, label_t>   l = get_training_labels();
 *  
 *  RandomForest<> rf;
 *
 *  // construct visitor to calculate out-of-bag error
 *  visitors::OOB_Error oob_v;
 *
 *  // perform training
 *  rf.learn(f, l, visitors::create_visitor(oob_v));
 *
 *  std::cout << "the out-of-bag error is: " << oob_v.oob_breiman << "\n";
 *      
 *  // get features for new data to be used for prediction
 *  MultiArrayView<2, feature_t> pf = get_features();
 *
 *  // allocate space for the response (pf.shape(0) is the number of samples)
 *  MultiArrayView<2, label_t> prediction(pf.shape(0), 1);
 *  MultiArrayView<2, double> prob(pf.shape(0), rf.class_count());
 *      
 *  // perform prediction on new data
 *  rf.predict_labels(pf, prediction);
 *  rf.predict_probabilities(pf, prob);
 *
 * 	\endcode
 *
 * 	Additional information such as Variable Importance measures are accessed
 * 	via Visitors defined in rf::visitors. 
 *  Have a look at rf::split for other splitting methods.
 *
*/
template <class LabelType = double , class PreprocessorTag = ClassificationTag >
class RandomForest
{

  public:
    //public typedefs
    typedef RandomForestOptions             Options_t;
    typedef detail::DecisionTree            DecisionTree_t;
    typedef ProblemSpec<LabelType>          ProblemSpec_t;
    typedef GiniSplit                       Default_Split_t;
    typedef EarlyStoppStd                   Default_Stop_t;
    typedef rf::visitors::StopVisiting      Default_Visitor_t;
    typedef  DT_StackEntry<ArrayVectorView<Int32>::iterator>
                    StackEntry_t;
    typedef LabelType                       LabelT; 
  protected:

    /** optimisation for predictLabels
     * */
    mutable MultiArray<2, double> garbage_prediction_;

  public:

    //problem independent data.
    Options_t                                   options_;
    //problem dependent data members - is only set if
    //a copy constructor, some sort of import
    //function or the learn function is called
    ArrayVector<DecisionTree_t>                 trees_;
    ProblemSpec_t                               ext_param_;
    /*mutable ArrayVector<int>                    tree_indices_;*/
    rf::visitors::OnlineLearnVisitor            online_visitor_;


    void reset()
    {
        ext_param_.clear();
        trees_.clear();
    }

  public:

    /** \name Constructors
     * Note: No copy Constructor specified as no pointers are manipulated
     * in this class
     */
    /*\{*/
    /**\brief default constructor
     *
     * \param options   general options to the Random Forest. Must be of Type
     *                  Options_t
     * \param ext_param problem specific values that can be supplied 
     *                  additionally. (class weights , labels etc)
     * \sa  RandomForestOptions, ProblemSpec
     *
     */
    RandomForest(Options_t const & options = Options_t(), 
                 ProblemSpec_t const & ext_param = ProblemSpec_t())
    :
        options_(options),
        ext_param_(ext_param)/*,
        tree_indices_(options.tree_count_,0)*/
    {
        /*for(int ii = 0 ; ii < int(tree_indices_.size()); ++ii)
            tree_indices_[ii] = ii;*/
    }

    /**\brief Create RF from external source
     * \param treeCount Number of trees to add.
     * \param topology_begin     
     * 					Iterator to a Container where the topology_ data
     *                  of the trees are stored.
     *                  Iterator should support at least treeCount forward 
     *                  iterations. (i.e. topology_end - topology_begin >= treeCount
     * \param parameter_begin  
     * 					iterator to a Container where the parameters_ data
     *                  of the trees are stored. Iterator should support at 
     *                  least treeCount forward iterations.
     * \param problem_spec 
     *                  Extrinsic parameters that specify the problem e.g.
     *                  ClassCount, featureCount etc.
     * \param options   (optional) specify options used to train the original
     *                  Random forest. This parameter is not used anywhere
     *                  during prediction and thus is optional.
     *
     */
     /* TODO: This constructor may be replaced by a Constructor using
     * NodeProxy iterators to encapsulate the underlying data type.
     */
    template<class TopologyIterator, class ParameterIterator>
    RandomForest(int                	   treeCount,
                  TopologyIterator         topology_begin,
                  ParameterIterator        parameter_begin,
                  ProblemSpec_t const & problem_spec,
                  Options_t const &     options = Options_t())
    :
        trees_(treeCount, DecisionTree_t(problem_spec)),
        ext_param_(problem_spec),
        options_(options)
    {
        for(unsigned int k=0; k<treeCount; ++k, ++topology_begin, ++parameter_begin)
        {
            trees_[k].topology_ = *topology_begin;
            trees_[k].parameters_ = *parameter_begin;
        }
    }

    /*\}*/


    /** \name Data Access
     * data access interface - usage of member variables is deprecated
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
     * set external parameters explicitly. 
     * If Random Forest has not been trained the preprocessor will 
     * either ignore filling values set this way or will throw an exception 
     * if values specified manually do not match the value calculated 
     & during the preparation step.
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

    /**\brief return number of features used while 
     * training.
     */
    int feature_count() const
    {
      return ext_param_.column_count_;
    }
    
    
    /**\brief return number of features used while 
     * training.
     *
     * deprecated. Use feature_count() instead.
     */
    int column_count() const
    {
      return ext_param_.column_count_;
    }

    /**\brief return number of classes used while 
     * training.
     */
    int class_count() const
    {
      return ext_param_.class_count_;
    }

    /**\brief return number of trees
     */
    int tree_count() const
    {
      return options_.tree_count_;
    }


    
    template<class U,class C1,
        class U2, class C2,
        class Split_t,
        class Stop_t,
        class Visitor_t,
        class Random_t>
    void onlineLearn(   MultiArrayView<2,U,C1> const & features,
                        MultiArrayView<2,U2,C2> const & response,
                        int new_start_index,
                        Visitor_t visitor_,
                        Split_t split_,
                        Stop_t stop_,
                        Random_t & random,
                        bool adjust_thresholds=false);

    template <class U, class C1, class U2,class C2>
    void onlineLearn(   MultiArrayView<2, U, C1> const  & features,
                        MultiArrayView<2, U2,C2> const  & labels,int new_start_index,bool adjust_thresholds=false)
    {
        RandomNumberGenerator<> rnd = RandomNumberGenerator<>(RandomSeed);
        onlineLearn(features, 
                    labels, 
                    new_start_index,
                    rf_default(), 
                    rf_default(), 
                    rf_default(),
                    rnd,
                    adjust_thresholds);
    }

    template<class U,class C1,
        class U2, class C2,
        class Split_t,
        class Stop_t,
        class Visitor_t,
        class Random_t>
    void reLearnTree(MultiArrayView<2,U,C1> const & features,
                     MultiArrayView<2,U2,C2> const & response,
                     int treeId,
                     Visitor_t visitor_,
                     Split_t split_,
                     Stop_t stop_,
                     Random_t & random);

    template<class U, class C1, class U2, class C2>
    void reLearnTree(MultiArrayView<2, U, C1> const & features,
                     MultiArrayView<2, U2, C2> const & labels,
                     int treeId)
    {
        RandomNumberGenerator<> rnd = RandomNumberGenerator<>(RandomSeed);
        reLearnTree(features,
                    labels,
                    treeId,
                    rf_default(),
                    rf_default(),
                    rf_default(),
                    rnd);
    }


    /**\name Learning
     * Following functions differ in the degree of customization
     * allowed
     */
    /*\{*/
    /**\brief learn on data with custom config and random number generator
     *
     * \param features  a N x M matrix containing N samples with M
     *                  features
     * \param response  a N x D matrix containing the corresponding
     *                  response. Current split functors assume D to
     *                  be 1 and ignore any additional columns.
     *                  This is not enforced to allow future support
     *                  for uncertain labels, label independent strata etc.
     *                  The Preprocessor specified during construction
     *                  should be able to handle features and labels
     *                  features and the labels.
     *  				see also: SplitFunctor, Preprocessing
     *
     * \param visitor   visitor which is to be applied after each split,
     *                  tree and at the end. Use rf_default for using
     *                  default value. (No Visitors)
     * 					see also: rf::visitors
     * \param split     split functor to be used to calculate each split
     *                  use rf_default() for using default value. (GiniSplit)
     * 					see also:  rf::split 
     * \param stop
     *                  predicate to be used to calculate each split
     *                  use rf_default() for using default value. (EarlyStoppStd)
     * \param random    RandomNumberGenerator to be used. Use
     *                  rf_default() to use default value.(RandomMT19337)
     *
     *
     */
    template <class U, class C1,
             class U2,class C2,
             class Split_t,
             class Stop_t,
             class Visitor_t,
             class Random_t>
    void learn( MultiArrayView<2, U, C1> const  &   features,
                MultiArrayView<2, U2,C2> const  &   response,
                Visitor_t                           visitor,
                Split_t                             split,
                Stop_t                              stop,
                Random_t                 const  &   random);

    template <class U, class C1,
             class U2,class C2,
             class Split_t,
             class Stop_t,
             class Visitor_t>
    void learn( MultiArrayView<2, U, C1> const  &   features,
                MultiArrayView<2, U2,C2> const  &   response,
                Visitor_t                           visitor,
                Split_t                             split,
                Stop_t                              stop)

    {
        RandomNumberGenerator<> rnd = RandomNumberGenerator<>(RandomSeed);
        learn(  features, 
                response,
                visitor, 
                split, 
                stop,
                rnd);
    }

    template <class U, class C1, class U2,class C2, class Visitor_t>
    void learn( MultiArrayView<2, U, C1> const  & features,
                MultiArrayView<2, U2,C2> const  & labels,
                Visitor_t                         visitor)
    {
        learn(  features, 
                labels, 
                visitor, 
                rf_default(), 
                rf_default());
    }

    template <class U, class C1, class U2,class C2, 
              class Visitor_t, class Split_t>
    void learn(   MultiArrayView<2, U, C1> const  & features,
                  MultiArrayView<2, U2,C2> const  & labels,
                  Visitor_t                         visitor,
                  Split_t                           split)
    {
        learn(  features, 
                labels, 
                visitor, 
                split, 
                rf_default());
    }

    /**\brief learn on data with default configuration
     *
     * \param features  a N x M matrix containing N samples with M
     *                  features
     * \param labels    a N x D matrix containing the corresponding
     *                  N labels. Current split functors assume D to
     *                  be 1 and ignore any additional columns.
     *                  this is not enforced to allow future support
     *                  for uncertain labels.
     *
     * learning is done with:
     *
     * \sa rf::split, EarlyStoppStd
     *
     * - Randomly seeded random number generator
     * - default gini split functor as described by Breiman
     * - default The standard early stopping criterion
     */
    template <class U, class C1, class U2,class C2>
    void learn(   MultiArrayView<2, U, C1> const  & features,
                    MultiArrayView<2, U2,C2> const  & labels)
    {
        learn(  features, 
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
     * \param stop: early stopping criterion
     * \return double value representing class. You can use the
     *         predictLabels() function together with the
     *         rf.external_parameter().class_type_ attribute
     *         to get back the same type used during learning. 
     */
    template <class U, class C, class Stop>
    LabelType predictLabel(MultiArrayView<2, U, C>const & features, Stop & stop) const;

    template <class U, class C>
    LabelType predictLabel(MultiArrayView<2, U, C>const & features)
    {
        return predictLabel(features, rf_default()); 
    } 
    /** \brief predict a label with features and class priors
     *
     * \param features: same as above.
     * \param prior:   iterator to prior weighting of classes
     * \return sam as above.
     */
    template <class U, class C>
    LabelType predictLabel(MultiArrayView<2, U, C> const & features,
                                ArrayVectorView<double> prior) const;

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
                       MultiArrayView<2, T, C2> & labels) const
    {
        vigra_precondition(features.shape(0) == labels.shape(0),
            "RandomForest::predictLabels(): Label array has wrong size.");
        for(int k=0; k<features.shape(0); ++k)
            labels(k,0) = detail::RequiresExplicitCast<T>::cast(predictLabel(rowVector(features, k), rf_default()));
    }

    template <class U, class C1, class T, class C2, class Stop>
    void predictLabels(MultiArrayView<2, U, C1>const & features,
                       MultiArrayView<2, T, C2> & labels,
                       Stop                     & stop) const
    {
        vigra_precondition(features.shape(0) == labels.shape(0),
            "RandomForest::predictLabels(): Label array has wrong size.");
        for(int k=0; k<features.shape(0); ++k)
            labels(k,0) = detail::RequiresExplicitCast<T>::cast(predictLabel(rowVector(features, k), stop));
    }
    /** \brief predict the class probabilities for multiple labels
     *
     *  \param features same as above
     *  \param prob a n x class_count_ matrix. passed by reference to
     *  save class probabilities
     *  \param stop earlystopping criterion
     *  \sa EarlyStopping
     */
    template <class U, class C1, class T, class C2, class Stop>
    void predictProbabilities(MultiArrayView<2, U, C1>const &   features,
                              MultiArrayView<2, T, C2> &        prob,
                              Stop                     &        stop) const;
    template <class T1,class T2, class C>
    void predictProbabilities(OnlinePredictionSet<T1> &  predictionSet,
                               MultiArrayView<2, T2, C> &       prob);

    /** \brief predict the class probabilities for multiple labels
     *
     *  \param features same as above
     *  \param prob a n x class_count_ matrix. passed by reference to
     *  save class probabilities
     */
    template <class U, class C1, class T, class C2>
    void predictProbabilities(MultiArrayView<2, U, C1>const &   features,
                              MultiArrayView<2, T, C2> &        prob)  const
    {
        predictProbabilities(features, prob, rf_default()); 
    }   

    template <class U, class C1, class T, class C2>
    void predictRaw(MultiArrayView<2, U, C1>const &   features,
                    MultiArrayView<2, T, C2> &        prob)  const;


    /*\}*/

};


template <class LabelType, class PreprocessorTag>
template<class U,class C1,
    class U2, class C2,
    class Split_t,
    class Stop_t,
    class Visitor_t,
    class Random_t>
void RandomForest<LabelType, PreprocessorTag>::onlineLearn(MultiArrayView<2,U,C1> const & features,
                                                             MultiArrayView<2,U2,C2> const & response,
                                                             int new_start_index,
                                                             Visitor_t visitor_,
                                                             Split_t split_,
                                                             Stop_t stop_,
                                                             Random_t & random,
                                                             bool adjust_thresholds)
{
    online_visitor_.activate();
    online_visitor_.adjust_thresholds=adjust_thresholds;

    using namespace rf;
    //typedefs
    typedef Processor<PreprocessorTag,LabelType,U,C1,U2,C2> Preprocessor_t;
    typedef          UniformIntRandomFunctor<Random_t>
                                                    RandFunctor_t;
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
    rf::visitors::StopVisiting stopvisiting;
    typedef  rf::visitors::detail::VisitorNode
                <rf::visitors::OnlineLearnVisitor, 
                 typename RF_CHOOSER(Visitor_t)::type> 
                                                        IntermedVis; 
    IntermedVis
        visitor(online_visitor_, RF_CHOOSER(Visitor_t)::choose(visitor_, stopvisiting));
    #undef RF_CHOOSER

    // Preprocess the data to get something the split functor can work
    // with. Also fill the ext_param structure by preprocessing
    // option parameters that could only be completely evaluated
    // when the training data is known.
    ext_param_.class_count_=0;
    Preprocessor_t preprocessor(    features, response,
                                    options_, ext_param_);

    // Make stl compatible random functor.
    RandFunctor_t           randint     ( random);

    // Give the Split functor information about the data.
    split.set_external_parameters(ext_param_);
    stop.set_external_parameters(ext_param_);


    //Create poisson samples
    PoissonSampler<RandomTT800> poisson_sampler(1.0,vigra::Int32(new_start_index),vigra::Int32(ext_param().row_count_));

    //TODO: visitors for online learning
    //visitor.visit_at_beginning(*this, preprocessor);

    // THE MAIN EFFING RF LOOP - YEAY DUDE!
    for(int ii = 0; ii < (int)trees_.size(); ++ii)
    {
        online_visitor_.tree_id=ii;
        poisson_sampler.sample();
        std::map<int,int> leaf_parents;
        leaf_parents.clear();
        //Get all the leaf nodes for that sample
        for(int s=0;s<poisson_sampler.numOfSamples();++s)
        {
            int sample=poisson_sampler[s];
            online_visitor_.current_label=preprocessor.response()(sample,0);
            online_visitor_.last_node_id=StackEntry_t::DecisionTreeNoParent;
            int leaf=trees_[ii].getToLeaf(rowVector(features,sample),online_visitor_);


            //Add to the list for that leaf
            online_visitor_.add_to_index_list(ii,leaf,sample);
            //TODO: Class count?
            //Store parent
            if(Node<e_ConstProbNode>(trees_[ii].topology_,trees_[ii].parameters_,leaf).prob_begin()[preprocessor.response()(sample,0)]!=1.0)
            {
                leaf_parents[leaf]=online_visitor_.last_node_id;
            }
        }


        std::map<int,int>::iterator leaf_iterator;
        for(leaf_iterator=leaf_parents.begin();leaf_iterator!=leaf_parents.end();++leaf_iterator)
        {
            int leaf=leaf_iterator->first;
            int parent=leaf_iterator->second;
            int lin_index=online_visitor_.trees_online_information[ii].exterior_to_index[leaf];
            ArrayVector<Int32> indeces;
            indeces.clear();
            indeces.swap(online_visitor_.trees_online_information[ii].index_lists[lin_index]);
            StackEntry_t stack_entry(indeces.begin(),
                                     indeces.end(),
                                     ext_param_.class_count_);


            if(parent!=-1)
            {
                if(NodeBase(trees_[ii].topology_,trees_[ii].parameters_,parent).child(0)==leaf)
                {
                    stack_entry.leftParent=parent;
                }
                else
                {
                    vigra_assert(NodeBase(trees_[ii].topology_,trees_[ii].parameters_,parent).child(1)==leaf,"last_node_id seems to be wrong");
                    stack_entry.rightParent=parent;
                }
            }
            //trees_[ii].continueLearn(preprocessor.features(),preprocessor.response(),stack_entry,split,stop,visitor,randint,leaf);
            trees_[ii].continueLearn(preprocessor.features(),preprocessor.response(),stack_entry,split,stop,visitor,randint,-1);
            //Now, the last one moved onto leaf
            online_visitor_.move_exterior_node(ii,trees_[ii].topology_.size(),ii,leaf);
            //Now it should be classified correctly!
        }

        /*visitor
            .visit_after_tree(  *this,
                                preprocessor,
                                poisson_sampler,
                                stack_entry,
                                ii);*/
    }

    //visitor.visit_at_end(*this, preprocessor);
    online_visitor_.deactivate();
}

template<class LabelType, class PreprocessorTag>
template<class U,class C1,
    class U2, class C2,
    class Split_t,
    class Stop_t,
    class Visitor_t,
    class Random_t>
void RandomForest<LabelType, PreprocessorTag>::reLearnTree(MultiArrayView<2,U,C1> const & features,
                 MultiArrayView<2,U2,C2> const & response,
                 int treeId,
                 Visitor_t visitor_,
                 Split_t split_,
                 Stop_t stop_,
                 Random_t & random)
{
    using namespace rf;
    
    
    typedef          UniformIntRandomFunctor<Random_t>
                                                    RandFunctor_t;

    // See rf_preprocessing.hxx for more info on this
    ext_param_.class_count_=0;
    typedef Processor<PreprocessorTag,LabelType, U, C1, U2, C2> Preprocessor_t;
    
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
    rf::visitors::StopVisiting stopvisiting;
    typedef  rf::visitors::detail::VisitorNode
                <rf::visitors::OnlineLearnVisitor, 
                typename RF_CHOOSER(Visitor_t)::type> IntermedVis; 
    IntermedVis
        visitor(online_visitor_, RF_CHOOSER(Visitor_t)::choose(visitor_, stopvisiting));
    #undef RF_CHOOSER
    vigra_precondition(options_.prepare_online_learning_,"reLearnTree: Re learning trees only makes sense, if online learning is enabled");
    online_visitor_.activate();

    // Make stl compatible random functor.
    RandFunctor_t           randint     ( random);

    // Preprocess the data to get something the split functor can work
    // with. Also fill the ext_param structure by preprocessing
    // option parameters that could only be completely evaluated
    // when the training data is known.
    Preprocessor_t preprocessor(    features, response,
                                    options_, ext_param_);

    // Give the Split functor information about the data.
    split.set_external_parameters(ext_param_);
    stop.set_external_parameters(ext_param_);

    /**\todo    replace this crappy class out. It uses function pointers.
     *          and is making code slower according to me.
     *          Comment from Nathan: This is copied from Rahul, so me=Rahul
     */
    Sampler<Random_t > sampler(preprocessor.strata().begin(),
                               preprocessor.strata().end(),
                               detail::make_sampler_opt(options_)
                                        .sampleSize(ext_param().actual_msample_),
                                    random);
    //initialize First region/node/stack entry
    sampler
        .sample();

    StackEntry_t
        first_stack_entry(  sampler.sampledIndices().begin(),
                            sampler.sampledIndices().end(),
                            ext_param_.class_count_);
    first_stack_entry
        .set_oob_range(     sampler.oobIndices().begin(),
                            sampler.oobIndices().end());
    online_visitor_.reset_tree(treeId);
    online_visitor_.tree_id=treeId;
    trees_[treeId].reset();
    trees_[treeId]
        .learn( preprocessor.features(),
                preprocessor.response(),
                first_stack_entry,
                split,
                stop,
                visitor,
                randint);
    visitor
        .visit_after_tree(  *this,
                            preprocessor,
                            sampler,
                            first_stack_entry,
                            treeId);

    online_visitor_.deactivate();
}

template <class LabelType, class PreprocessorTag>
template <class U, class C1,
         class U2,class C2,
         class Split_t,
         class Stop_t,
         class Visitor_t,
         class Random_t>
void RandomForest<LabelType, PreprocessorTag>::
                     learn( MultiArrayView<2, U, C1> const  &   features,
                            MultiArrayView<2, U2,C2> const  &   response,
                            Visitor_t                           visitor_,
                            Split_t                             split_,
                            Stop_t                              stop_,
                            Random_t                 const  &   random)
{
    using namespace rf;
    //this->reset();
    //typedefs
    typedef          UniformIntRandomFunctor<Random_t>
                                                    RandFunctor_t;

    // See rf_preprocessing.hxx for more info on this
    typedef Processor<PreprocessorTag,LabelType, U, C1, U2, C2> Preprocessor_t;

    vigra_precondition(features.shape(0) == response.shape(0),
        "RandomForest::learn(): shape mismatch between features and response.");
    
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
    rf::visitors::StopVisiting stopvisiting;
    typedef  rf::visitors::detail::VisitorNode<
                rf::visitors::OnlineLearnVisitor, 
                typename RF_CHOOSER(Visitor_t)::type> IntermedVis; 
    IntermedVis
        visitor(online_visitor_, RF_CHOOSER(Visitor_t)::choose(visitor_, stopvisiting));
    #undef RF_CHOOSER
    if(options_.prepare_online_learning_)
        online_visitor_.activate();
    else
        online_visitor_.deactivate();


    // Make stl compatible random functor.
    RandFunctor_t           randint     ( random);


    // Preprocess the data to get something the split functor can work
    // with. Also fill the ext_param structure by preprocessing
    // option parameters that could only be completely evaluated
    // when the training data is known.
    Preprocessor_t preprocessor(    features, response,
                                    options_, ext_param_);

    // Give the Split functor information about the data.
    split.set_external_parameters(ext_param_);
    stop.set_external_parameters(ext_param_);


    //initialize trees.
    trees_.resize(options_.tree_count_  , DecisionTree_t(ext_param_));

    Sampler<Random_t > sampler(preprocessor.strata().begin(),
                               preprocessor.strata().end(),
                               detail::make_sampler_opt(options_)
                                        .sampleSize(ext_param().actual_msample_),
                                    random);

    visitor.visit_at_beginning(*this, preprocessor);
    // THE MAIN EFFING RF LOOP - YEAY DUDE!
    
    for(int ii = 0; ii < (int)trees_.size(); ++ii)
    {
        //initialize First region/node/stack entry
        sampler
            .sample();  
        StackEntry_t
            first_stack_entry(  sampler.sampledIndices().begin(),
                                sampler.sampledIndices().end(),
                                ext_param_.class_count_);
        first_stack_entry
            .set_oob_range(     sampler.oobIndices().begin(),
                                sampler.oobIndices().end());
        trees_[ii]
            .learn(             preprocessor.features(),
                                preprocessor.response(),
                                first_stack_entry,
                                split,
                                stop,
                                visitor,
                                randint);
        visitor
            .visit_after_tree(  *this,
                                preprocessor,
                                sampler,
                                first_stack_entry,
                                ii);
    }

    visitor.visit_at_end(*this, preprocessor);
    // Only for online learning?
    online_visitor_.deactivate();
}




template <class LabelType, class Tag>
template <class U, class C, class Stop>
LabelType RandomForest<LabelType, Tag>
    ::predictLabel(MultiArrayView<2, U, C> const & features, Stop & stop) const
{
    vigra_precondition(columnCount(features) >= ext_param_.column_count_,
        "RandomForestn::predictLabel():"
            " Too few columns in feature matrix.");
    vigra_precondition(rowCount(features) == 1,
        "RandomForestn::predictLabel():"
            " Feature matrix must have a singlerow.");
    typedef MultiArrayShape<2>::type Shp;
    garbage_prediction_.reshape(Shp(1, ext_param_.class_count_), 0.0);
    LabelType          d;
    predictProbabilities(features, garbage_prediction_, stop);
    ext_param_.to_classlabel(argMax(garbage_prediction_), d);
    return d;
}


//Same thing as above with priors for each label !!!
template <class LabelType, class PreprocessorTag>
template <class U, class C>
LabelType RandomForest<LabelType, PreprocessorTag>
    ::predictLabel( MultiArrayView<2, U, C> const & features,
                    ArrayVectorView<double> priors) const
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
                    priors.begin(), prob.begin(),
                    Arg1()*Arg2());
    LabelType          d;
    ext_param_.to_classlabel(argMax(prob), d);
    return d;
}

template<class LabelType,class PreprocessorTag>
template <class T1,class T2, class C>
void RandomForest<LabelType,PreprocessorTag>
    ::predictProbabilities(OnlinePredictionSet<T1> &  predictionSet,
                          MultiArrayView<2, T2, C> &       prob)
{
    //Features are n xp
    //prob is n x NumOfLabel probability for each feature in each class
    
    vigra_precondition(rowCount(predictionSet.features) == rowCount(prob),
                       "RandomFroest::predictProbabilities():"
                       " Feature matrix and probability matrix size mismatch.");
    // num of features must be bigger than num of features in Random forest training
    // but why bigger?
    vigra_precondition( columnCount(predictionSet.features) >= ext_param_.column_count_,
      "RandomForestn::predictProbabilities():"
        " Too few columns in feature matrix.");
    vigra_precondition( columnCount(prob)
                        == (MultiArrayIndex)ext_param_.class_count_,
      "RandomForestn::predictProbabilities():"
      " Probability matrix must have as many columns as there are classes.");
    prob.init(0.0);
    //store total weights
    std::vector<T1> totalWeights(predictionSet.indices[0].size(),0.0);
    //Go through all trees
    int set_id=-1;
    for(int k=0; k<options_.tree_count_; ++k)
    {
        set_id=(set_id+1) % predictionSet.indices[0].size();
        typedef std::set<SampleRange<T1> > my_set;
        typedef typename my_set::iterator set_it;
        //typedef std::set<std::pair<int,SampleRange<T1> > >::iterator set_it;
        //Build a stack with all the ranges we have
        std::vector<std::pair<int,set_it> > stack;
        stack.clear();
        for(set_it i=predictionSet.ranges[set_id].begin();
             i!=predictionSet.ranges[set_id].end();++i)
            stack.push_back(std::pair<int,set_it>(2,i));
        //get weights predicted by single tree
        int num_decisions=0;
        while(!stack.empty())
        {
            set_it range=stack.back().second;
            int index=stack.back().first;
            stack.pop_back();
            ++num_decisions;

            if(trees_[k].isLeafNode(trees_[k].topology_[index]))
            {
                ArrayVector<double>::iterator weights=Node<e_ConstProbNode>(trees_[k].topology_,
                                                                            trees_[k].parameters_,
                                                                            index).prob_begin();
                for(int i=range->start;i!=range->end;++i)
                {
                    //update votecount.
                    for(int l=0; l<ext_param_.class_count_; ++l)
                    {
                        prob(predictionSet.indices[set_id][i], l) += (T2)weights[l];
                        //every weight in totalWeight.
                        totalWeights[predictionSet.indices[set_id][i]] += (T1)weights[l];
                    }
                }
            }

            else
            {
                if(trees_[k].topology_[index]!=i_ThresholdNode)
                {
                    throw std::runtime_error("predicting with online prediction sets is only supported for RFs with threshold nodes");
                }
                Node<i_ThresholdNode> node(trees_[k].topology_,trees_[k].parameters_,index);
                if(range->min_boundaries[node.column()]>=node.threshold())
                {
                    //Everything goes to right child
                    stack.push_back(std::pair<int,set_it>(node.child(1),range));
                    continue;
                }
                if(range->max_boundaries[node.column()]<node.threshold())
                {
                    //Everything goes to the left child
                    stack.push_back(std::pair<int,set_it>(node.child(0),range));
                    continue;
                }
                //We have to split at this node
                SampleRange<T1> new_range=*range;
                new_range.min_boundaries[node.column()]=FLT_MAX;
                range->max_boundaries[node.column()]=-FLT_MAX;
                new_range.start=new_range.end=range->end;
                int i=range->start;
                while(i!=range->end)
                {
                    //Decide for range->indices[i]
                    if(predictionSet.features(predictionSet.indices[set_id][i],node.column())>=node.threshold())
                    {
                        new_range.min_boundaries[node.column()]=std::min(new_range.min_boundaries[node.column()],
                                                                    predictionSet.features(predictionSet.indices[set_id][i],node.column()));
                        --range->end;
                        --new_range.start;
                        std::swap(predictionSet.indices[set_id][i],predictionSet.indices[set_id][range->end]);

                    }
                    else
                    {
                        range->max_boundaries[node.column()]=std::max(range->max_boundaries[node.column()],
                                                                 predictionSet.features(predictionSet.indices[set_id][i],node.column()));
                        ++i;
                    }
                }
                //The old one ...
                if(range->start==range->end)
                {
                    predictionSet.ranges[set_id].erase(range);
                }
                else
                {
                    stack.push_back(std::pair<int,set_it>(node.child(0),range));
                }
                //And the new one ...
                if(new_range.start!=new_range.end)
                {
                    std::pair<set_it,bool> new_it=predictionSet.ranges[set_id].insert(new_range);
                    stack.push_back(std::pair<int,set_it>(node.child(1),new_it.first));
                }
            }
        }
        predictionSet.cumulativePredTime[k]=num_decisions;
    }
    for(unsigned int i=0;i<totalWeights.size();++i)
    {
        double test=0.0;
        //Normalise votes in each row by total VoteCount (totalWeight
        for(int l=0; l<ext_param_.class_count_; ++l)
        {
            test+=prob(i,l);
            prob(i, l) /= totalWeights[i];
        }
        assert(test==totalWeights[i]);
        assert(totalWeights[i]>0.0);
    }
}

template <class LabelType, class PreprocessorTag>
template <class U, class C1, class T, class C2, class Stop_t>
void RandomForest<LabelType, PreprocessorTag>
    ::predictProbabilities(MultiArrayView<2, U, C1>const &  features,
                           MultiArrayView<2, T, C2> &       prob,
                           Stop_t                   &       stop_) const
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

    #define RF_CHOOSER(type_) detail::Value_Chooser<type_, Default_##type_> 
    Default_Stop_t default_stop(options_);
    typename RF_CHOOSER(Stop_t)::type & stop
            = RF_CHOOSER(Stop_t)::choose(stop_, default_stop); 
    #undef RF_CHOOSER 
    stop.set_external_parameters(ext_param_, tree_count());
    prob.init(NumericTraits<T>::zero());
    /* This code was originally there for testing early stopping
     * - we wanted the order of the trees to be randomized
    if(tree_indices_.size() != 0)
    {
       std::random_shuffle(tree_indices_.begin(),
                           tree_indices_.end()); 
    }
    */
    //Classify for each row.
    for(int row=0; row < rowCount(features); ++row)
    {
        ArrayVector<double>::const_iterator weights;

        //totalWeight == totalVoteCount!
        double totalWeight = 0.0;

        //Let each tree classify...
        for(int k=0; k<options_.tree_count_; ++k)
        {
            //get weights predicted by single tree
            weights = trees_[k /*tree_indices_[k]*/].predict(rowVector(features, row));

            //update votecount.
            int weighted = options_.predict_weighted_;
            for(int l=0; l<ext_param_.class_count_; ++l)
            {
                double cur_w = weights[l] * (weighted * (*(weights-1))
                                           + (1-weighted));
                prob(row, l) += (T)cur_w;
                //every weight in totalWeight.
                totalWeight += cur_w;
            }
            if(stop.after_prediction(weights, 
                                     k,
                                     rowVector(prob, row),
                                     totalWeight))
            {
                break;
            }
        }

        //Normalise votes in each row by total VoteCount (totalWeight
        for(int l=0; l< ext_param_.class_count_; ++l)
        {
            prob(row, l) /= detail::RequiresExplicitCast<T>::cast(totalWeight);
        }
    }

}

template <class LabelType, class PreprocessorTag>
template <class U, class C1, class T, class C2>
void RandomForest<LabelType, PreprocessorTag>
    ::predictRaw(MultiArrayView<2, U, C1>const &  features,
                           MultiArrayView<2, T, C2> &       prob) const
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

    #define RF_CHOOSER(type_) detail::Value_Chooser<type_, Default_##type_> 
    prob.init(NumericTraits<T>::zero());
    /* This code was originally there for testing early stopping
     * - we wanted the order of the trees to be randomized
    if(tree_indices_.size() != 0)
    {
       std::random_shuffle(tree_indices_.begin(),
                           tree_indices_.end()); 
    }
    */
    //Classify for each row.
    for(int row=0; row < rowCount(features); ++row)
    {
        ArrayVector<double>::const_iterator weights;

        //totalWeight == totalVoteCount!
        double totalWeight = 0.0;

        //Let each tree classify...
        for(int k=0; k<options_.tree_count_; ++k)
        {
            //get weights predicted by single tree
            weights = trees_[k /*tree_indices_[k]*/].predict(rowVector(features, row));

            //update votecount.
            int weighted = options_.predict_weighted_;
            for(int l=0; l<ext_param_.class_count_; ++l)
            {
                double cur_w = weights[l] * (weighted * (*(weights-1))
                                           + (1-weighted));
                prob(row, l) += (T)cur_w;
                //every weight in totalWeight.
                totalWeight += cur_w;
            }
        }
    }
    prob/= options_.tree_count_;

}

//@}

} // namespace vigra

#include "random_forest/rf_algorithm.hxx"
#endif // VIGRA_RANDOM_FOREST_HXX
