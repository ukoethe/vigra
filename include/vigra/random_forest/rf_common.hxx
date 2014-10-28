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


#ifndef VIGRA_RF_COMMON_HXX
#define VIGRA_RF_COMMON_HXX

namespace vigra
{


struct ClassificationTag
{};

struct RegressionTag
{};

namespace detail
{
    class RF_DEFAULT;
}
inline detail::RF_DEFAULT& rf_default();
namespace detail
{

/* \brief singleton default tag class -
 *
 *  use the rf_default() factory function to use the tag.
 *  \sa RandomForest<>::learn();
 */
class RF_DEFAULT
{
    private:
        RF_DEFAULT()
        {}
    public:
        friend RF_DEFAULT& ::vigra::rf_default();

        /** ok workaround for automatic choice of the decisiontree
         * stackentry.
         */
};

/* \brief chooses between default type and type supplied
 * 
 * This is an internal class and you shouldn't really care about it.
 * Just pass on used in RandomForest.learn()
 * Usage:
 *\code
 *      // example: use container type supplied by user or ArrayVector if 
 *      //          rf_default() was specified as argument;
 *      template<class Container_t>
 *      void do_some_foo(Container_t in)
 *      {
 *          typedef ArrayVector<int>    Default_Container_t;
 *          Default_Container_t         default_value;
 *          Value_Chooser<Container_t,  Default_Container_t> 
 *                      choose(in, default_value);
 *
 *          // if the user didn't care and the in was of type 
 *          // RF_DEFAULT then default_value is used.
 *          do_some_more_foo(choose.value());
 *      }
 *      Value_Chooser choose_val<Type, Default_Type>
 *\endcode
 */
template<class T, class C>
class Value_Chooser
{
public:
    typedef T type;
    static T & choose(T & t, C &)
    {
        return t; 
    }
};

template<class C>
class Value_Chooser<detail::RF_DEFAULT, C>
{
public:
    typedef C type;
    
    static C & choose(detail::RF_DEFAULT &, C & c)
    {
        return c; 
    }
};




} //namespace detail


/**\brief factory function to return a RF_DEFAULT tag
 * \sa RandomForest<>::learn()
 */
detail::RF_DEFAULT& rf_default()
{
    static detail::RF_DEFAULT result;
    return result;
}

/** tags used with the RandomForestOptions class
 * \sa RF_Traits::Option_t
 */
enum RF_OptionTag   { RF_EQUAL,
                      RF_PROPORTIONAL,
                      RF_EXTERNAL,
                      RF_NONE,
                      RF_FUNCTION,
                      RF_LOG,
                      RF_SQRT,
                      RF_CONST,
                      RF_ALL};


/** \addtogroup MachineLearning 
**/
//@{

/**\brief Options object for the random forest
 *
 * usage:
 * RandomForestOptions a =  RandomForestOptions()
 *                              .param1(value1)
 *                              .param2(value2)
 *                              ...
 *
 * This class only contains options/parameters that are not problem
 * dependent. The ProblemSpec class contains methods to set class weights
 * if necessary.
 *
 * Note that the return value of all methods is *this which makes
 * concatenating of options as above possible.
 */
class RandomForestOptions
{
  public:
    /**\name sampling options*/
    /*\{*/
    // look at the member access functions for documentation
    double  training_set_proportion_;
    int     training_set_size_;
    int (*training_set_func_)(int);
    RF_OptionTag
        training_set_calc_switch_;

    bool    sample_with_replacement_;
    RF_OptionTag
            stratification_method_;


    /**\name general random forest options
     *
     * these usually will be used by most split functors and
     * stopping predicates
     */
    /*\{*/
    RF_OptionTag    mtry_switch_;
    int     mtry_;
    int (*mtry_func_)(int) ;

    bool predict_weighted_; 
    int tree_count_;
    int min_split_node_size_;
    int max_tree_depth_;
    bool prepare_online_learning_;
    /*\}*/

    typedef ArrayVector<double> double_array;
    typedef std::map<std::string, double_array> map_type;

    int serialized_size() const
    {
        return 12;
    }
    

    bool operator==(RandomForestOptions & rhs) const
    {
        bool result = true;
        #define COMPARE(field) result = result && (this->field == rhs.field); 
        COMPARE(training_set_proportion_);
        COMPARE(training_set_size_);
        COMPARE(training_set_calc_switch_);
        COMPARE(sample_with_replacement_);
        COMPARE(stratification_method_);
        COMPARE(mtry_switch_);
        COMPARE(mtry_);
        COMPARE(tree_count_);
        COMPARE(min_split_node_size_);
        COMPARE(predict_weighted_);
        #undef COMPARE

        return result;
    }
    bool operator!=(RandomForestOptions & rhs_) const
    {
        return !(*this == rhs_);
    }
    template<class Iter>
    void unserialize(Iter const & begin, Iter const & end)
    {
        Iter iter = begin;
        vigra_precondition(static_cast<int>(end - begin) == serialized_size(), 
                           "RandomForestOptions::unserialize():"
                           "wrong number of parameters");
        #define PULL(item_, type_) item_ = type_(*iter); ++iter;
        PULL(training_set_proportion_, double);
        PULL(training_set_size_, int);
        ++iter; //PULL(training_set_func_, double);
        PULL(training_set_calc_switch_, (RF_OptionTag)int);
        PULL(sample_with_replacement_, 0 != );
        PULL(stratification_method_, (RF_OptionTag)int);
        PULL(mtry_switch_, (RF_OptionTag)int);
        PULL(mtry_, int);
        ++iter; //PULL(mtry_func_, double);
        PULL(tree_count_, int);
        PULL(min_split_node_size_, int);
        PULL(predict_weighted_, 0 !=);
        #undef PULL
    }
    template<class Iter>
    void serialize(Iter const &  begin, Iter const & end) const
    {
        Iter iter = begin;
        vigra_precondition(static_cast<int>(end - begin) == serialized_size(), 
                           "RandomForestOptions::serialize():"
                           "wrong number of parameters");
        #define PUSH(item_) *iter = double(item_); ++iter;
        PUSH(training_set_proportion_);
        PUSH(training_set_size_);
        if(training_set_func_ != 0)
        {
            PUSH(1);
        }
        else
        {
            PUSH(0);
        }
        PUSH(training_set_calc_switch_);
        PUSH(sample_with_replacement_);
        PUSH(stratification_method_);
        PUSH(mtry_switch_);
        PUSH(mtry_);
        if(mtry_func_ != 0)
        {
            PUSH(1);
        }
        else
        {
            PUSH(0);
        }
        PUSH(tree_count_);
        PUSH(min_split_node_size_);
        PUSH(predict_weighted_);
        #undef PUSH
    }
    
    void make_from_map(map_type & in) // -> const: .operator[] -> .find
    {
        #define PULL(item_, type_) item_ = type_(in[#item_][0]); 
        #define PULLBOOL(item_, type_) item_ = type_(in[#item_][0] > 0); 
        PULL(training_set_proportion_,double);
        PULL(training_set_size_, int);
        PULL(mtry_, int);
        PULL(tree_count_, int);
        PULL(min_split_node_size_, int);
        PULLBOOL(sample_with_replacement_, bool);
        PULLBOOL(prepare_online_learning_, bool);
        PULLBOOL(predict_weighted_, bool);
        
        PULL(training_set_calc_switch_, (RF_OptionTag)(int));

        PULL(stratification_method_, (RF_OptionTag)(int));
        PULL(mtry_switch_, (RF_OptionTag)(int));
        
        /*don't pull*/
        //PULL(mtry_func_!=0, int);
        //PULL(training_set_func,int);
        #undef PULL
        #undef PULLBOOL
    }
    void make_map(map_type & in) const
    {
        #define PUSH(item_, type_) in[#item_] = double_array(1, double(item_));
        #define PUSHFUNC(item_, type_) in[#item_] = double_array(1, double(item_!=0));
        PUSH(training_set_proportion_,double);
        PUSH(training_set_size_, int);
        PUSH(mtry_, int);
        PUSH(tree_count_, int);
        PUSH(min_split_node_size_, int);
        PUSH(sample_with_replacement_, bool);
        PUSH(prepare_online_learning_, bool);
        PUSH(predict_weighted_, bool);
        
        PUSH(training_set_calc_switch_, RF_OptionTag);
        PUSH(stratification_method_, RF_OptionTag);
        PUSH(mtry_switch_, RF_OptionTag);
        
        PUSHFUNC(mtry_func_, int);
        PUSHFUNC(training_set_func_,int);
        #undef PUSH
        #undef PUSHFUNC
    }


    /**\brief create a RandomForestOptions object with default initialisation.
     *
     * look at the other member functions for more information on default
     * values
     */
    RandomForestOptions()
    :
        training_set_proportion_(1.0),
        training_set_size_(0),
        training_set_func_(0),
        training_set_calc_switch_(RF_PROPORTIONAL),
        sample_with_replacement_(true),
        stratification_method_(RF_NONE),
        mtry_switch_(RF_SQRT),
        mtry_(0),
        mtry_func_(0),
        predict_weighted_(false),
        tree_count_(256),
        min_split_node_size_(1),
        max_tree_depth_(0),
        prepare_online_learning_(false)
    {}

    /**\brief specify stratification strategy
     *
     * default: RF_NONE
     * possible values: RF_EQUAL, RF_PROPORTIONAL,
     *                  RF_EXTERNAL, RF_NONE
     * RF_EQUAL:        get equal amount of samples per class.
     * RF_PROPORTIONAL: sample proportional to fraction of class samples
     *                  in population
     * RF_EXTERNAL:     strata_weights_ field of the ProblemSpec_t object
     *                  has been set externally. (defunct)
     */
    RandomForestOptions & use_stratification(RF_OptionTag in)
    {
        vigra_precondition(in == RF_EQUAL ||
                           in == RF_PROPORTIONAL ||
                           in == RF_EXTERNAL ||
                           in == RF_NONE,
                           "RandomForestOptions::use_stratification()"
                           "input must be RF_EQUAL, RF_PROPORTIONAL,"
                           "RF_EXTERNAL or RF_NONE");
        stratification_method_ = in;
        return *this;
    }

    RandomForestOptions & prepare_online_learning(bool in)
    {
        prepare_online_learning_=in;
        return *this;
    }

    /**\brief sample from training population with or without replacement?
     *
     * <br> Default: true
     */
    RandomForestOptions & sample_with_replacement(bool in)
    {
        sample_with_replacement_ = in;
        return *this;
    }

    /**\brief  specify the fraction of the total number of samples 
     * used per tree for learning. 
     *
     * This value should be in [0.0 1.0] if sampling without
     * replacement has been specified.
     *
     * <br> default : 1.0
     */
    RandomForestOptions & samples_per_tree(double in)
    {
        training_set_proportion_ = in;
        training_set_calc_switch_ = RF_PROPORTIONAL;
        return *this;
    }

    /**\brief directly specify the number of samples per tree
     */
    RandomForestOptions & samples_per_tree(int in)
    {
        training_set_size_ = in;
        training_set_calc_switch_ = RF_CONST;
        return *this;
    }

    /**\brief use external function to calculate the number of samples each
     *        tree should be learnt with.
     *
     * \param in function pointer that takes the number of rows in the
     *           learning data and outputs the number samples per tree.
     */
    RandomForestOptions & samples_per_tree(int (*in)(int))
    {
        training_set_func_ = in;
        training_set_calc_switch_ = RF_FUNCTION;
        return *this;
    }
    
    /**\brief weight each tree with number of samples in that node
     */
    RandomForestOptions & predict_weighted()
    {
        predict_weighted_ = true;
        return *this;
    }

    /**\brief use built in mapping to calculate mtry
     *
     * Use one of the built in mappings to calculate mtry from the number
     * of columns in the input feature data.
     * \param in possible values: RF_LOG, RF_SQRT or RF_ALL
     *           <br> default: RF_SQRT.
     */
    RandomForestOptions & features_per_node(RF_OptionTag in)
    {
        vigra_precondition(in == RF_LOG ||
                           in == RF_SQRT||
                           in == RF_ALL,
                           "RandomForestOptions()::features_per_node():"
                           "input must be of type RF_LOG or RF_SQRT");
        mtry_switch_ = in;
        return *this;
    }

    /**\brief Set mtry to a constant value
     *
     * mtry is the number of columns/variates/variables randomly chosen
     * to select the best split from.
     *
     */
    RandomForestOptions & features_per_node(int in)
    {
        mtry_ = in;
        mtry_switch_ = RF_CONST;
        return *this;
    }

    /**\brief use a external function to calculate mtry
     *
     * \param in function pointer that takes int (number of columns
     *           of the and outputs int (mtry)
     */
    RandomForestOptions & features_per_node(int(*in)(int))
    {
        mtry_func_ = in;
        mtry_switch_ = RF_FUNCTION;
        return *this;
    }

    /** How many trees to create?
     *
     * <br> Default: 255.
     */
    RandomForestOptions & tree_count(int in)
    {
        tree_count_ = in;
        return *this;
    }

    /**\brief Number of examples required for a node to be split.
     *
     *  When the number of examples in a node is below this number,
     *  the node is not split even if class separation is not yet perfect.
     *  Instead, the node returns the proportion of each class
     *  (among the remaining examples) during the prediction phase.
     *  <br> Default: 1 (complete growing)
     */
    RandomForestOptions & min_split_node_size(int in)
    {
        min_split_node_size_ = in;
        return *this;
    }

    /**\brief maximum depth of a decision tree.
     *
     *  When the depth of a node is equal to this number
     *  the node is not split even if class separation is not yet perfect.
     *  Instead, the node returns the proportion of each class
     *  (among the remaining examples) during the prediction phase.
     *  <br> Default: 0 (complete growing)
     */
    RandomForestOptions & max_tree_depth(int in)
    {
        max_tree_depth_ = in;
        return *this;
    }
};


/** \brief problem types 
 */
enum Problem_t{REGRESSION, CLASSIFICATION, CHECKLATER};


/** \brief problem specification class for the random forest.
 *
 * This class contains all the problem specific parameters the random
 * forest needs for learning. Specification of an instance of this class
 * is optional as all necessary fields will be computed prior to learning
 * if not specified.
 *
 * if needed usage is similar to that of RandomForestOptions
 */

template<class LabelType = double>
class ProblemSpec
{


public:

    /** \brief  problem class
     */

    typedef LabelType       Label_t;
    ArrayVector<Label_t>    classes;
    typedef ArrayVector<double>                 double_array;
    typedef std::map<std::string, double_array> map_type;

    int                     column_count_;    // number of features
    int                     class_count_;     // number of classes
    int                     row_count_;       // number of samples

    int                     actual_mtry_;     // mtry used in training
    int                     actual_msample_;  // number if in-bag samples per tree

    Problem_t               problem_type_;    // classification or regression
    
    int used_;                                // this ProblemSpec is valid
    ArrayVector<double>     class_weights_;   // if classes have different importance
    int                     is_weighted_;     // class_weights_ are used
    double                  precision_;       // termination criterion for regression loss
    int                     response_size_; 
        
    template<class T> 
    void to_classlabel(int index, T & out) const
    {
        out = T(classes[index]);
    }
    template<class T> 
    int to_classIndex(T index) const
    {
        return std::find(classes.begin(), classes.end(), index) - classes.begin();
    }

    #define EQUALS(field) field(rhs.field)
    ProblemSpec(ProblemSpec const & rhs)
    : 
        EQUALS(column_count_),
        EQUALS(class_count_),
        EQUALS(row_count_),
        EQUALS(actual_mtry_),
        EQUALS(actual_msample_),
        EQUALS(problem_type_),
        EQUALS(used_),
        EQUALS(class_weights_),
        EQUALS(is_weighted_),
        EQUALS(precision_),
        EQUALS(response_size_)
    {
        std::back_insert_iterator<ArrayVector<Label_t> >
                        iter(classes);
        std::copy(rhs.classes.begin(), rhs.classes.end(), iter); 
    }
    #undef EQUALS
    #define EQUALS(field) field(rhs.field)
    template<class T>
    ProblemSpec(ProblemSpec<T> const & rhs)
    : 
        EQUALS(column_count_),
        EQUALS(class_count_),
        EQUALS(row_count_),
        EQUALS(actual_mtry_),
        EQUALS(actual_msample_),
        EQUALS(problem_type_),
        EQUALS(used_),
        EQUALS(class_weights_),
        EQUALS(is_weighted_),
        EQUALS(precision_),
        EQUALS(response_size_)
    {
        std::back_insert_iterator<ArrayVector<Label_t> >
                        iter(classes);
        std::copy(rhs.classes.begin(), rhs.classes.end(), iter); 
    }
    #undef EQUALS

    #define EQUALS(field) (this->field = rhs.field);
    ProblemSpec & operator=(ProblemSpec const & rhs)
    {
        EQUALS(column_count_);
        EQUALS(class_count_);
        EQUALS(row_count_);
        EQUALS(actual_mtry_);
        EQUALS(actual_msample_);
        EQUALS(problem_type_);
        EQUALS(used_);
        EQUALS(is_weighted_);
        EQUALS(precision_);
        EQUALS(response_size_)
        class_weights_.clear();
        std::back_insert_iterator<ArrayVector<double> >
                        iter2(class_weights_);
        std::copy(rhs.class_weights_.begin(), rhs.class_weights_.end(), iter2); 
        classes.clear();
        std::back_insert_iterator<ArrayVector<Label_t> >
                        iter(classes);
        std::copy(rhs.classes.begin(), rhs.classes.end(), iter); 
        return *this;
    }

    template<class T>
    ProblemSpec<Label_t> & operator=(ProblemSpec<T> const & rhs)
    {
        EQUALS(column_count_);
        EQUALS(class_count_);
        EQUALS(row_count_);
        EQUALS(actual_mtry_);
        EQUALS(actual_msample_);
        EQUALS(problem_type_);
        EQUALS(used_);
        EQUALS(is_weighted_);
        EQUALS(precision_);
        EQUALS(response_size_)
        class_weights_.clear();
        std::back_insert_iterator<ArrayVector<double> >
                        iter2(class_weights_);
        std::copy(rhs.class_weights_.begin(), rhs.class_weights_.end(), iter2); 
        classes.clear();
        std::back_insert_iterator<ArrayVector<Label_t> >
                        iter(classes);
        std::copy(rhs.classes.begin(), rhs.classes.end(), iter); 
        return *this;
    }
    #undef EQUALS

    template<class T>
    bool operator==(ProblemSpec<T> const & rhs)
    {
        bool result = true;
        #define COMPARE(field) result = result && (this->field == rhs.field);
        COMPARE(column_count_);
        COMPARE(class_count_);
        COMPARE(row_count_);
        COMPARE(actual_mtry_);
        COMPARE(actual_msample_);
        COMPARE(problem_type_);
        COMPARE(is_weighted_);
        COMPARE(precision_);
        COMPARE(used_);
        COMPARE(class_weights_);
        COMPARE(classes);
        COMPARE(response_size_)
        #undef COMPARE
        return result;
    }

    bool operator!=(ProblemSpec & rhs)
    {
        return !(*this == rhs);
    }


    size_t serialized_size() const
    {
        return 10 + class_count_ *int(is_weighted_+1);
    }


    template<class Iter>
    void unserialize(Iter const & begin, Iter const & end)
    {
        Iter iter = begin;
        vigra_precondition(end - begin >= 10, 
                           "ProblemSpec::unserialize():"
                           "wrong number of parameters");
        #define PULL(item_, type_) item_ = type_(*iter); ++iter;
        PULL(column_count_,int);
        PULL(class_count_, int);

        vigra_precondition(end - begin >= 10 + class_count_, 
                           "ProblemSpec::unserialize(): 1");
        PULL(row_count_, int);
        PULL(actual_mtry_,int);
        PULL(actual_msample_, int);
        PULL(problem_type_, Problem_t);
        PULL(is_weighted_, int);
        PULL(used_, int);
        PULL(precision_, double);
        PULL(response_size_, int);
        if(is_weighted_)
        {
            vigra_precondition(end - begin == 10 + 2*class_count_, 
                               "ProblemSpec::unserialize(): 2");
            class_weights_.insert(class_weights_.end(),
                                  iter, 
                                  iter + class_count_);
            iter += class_count_; 
        }
        classes.insert(classes.end(), iter, end);
        #undef PULL
    }


    template<class Iter>
    void serialize(Iter const & begin, Iter const & end) const
    {
        Iter iter = begin;
        vigra_precondition(end - begin == serialized_size(), 
                           "RandomForestOptions::serialize():"
                           "wrong number of parameters");
        #define PUSH(item_) *iter = double(item_); ++iter;
        PUSH(column_count_);
        PUSH(class_count_)
        PUSH(row_count_);
        PUSH(actual_mtry_);
        PUSH(actual_msample_);
        PUSH(problem_type_);
        PUSH(is_weighted_);
        PUSH(used_);
        PUSH(precision_);
        PUSH(response_size_);
        if(is_weighted_)
        {
            std::copy(class_weights_.begin(),
                      class_weights_.end(),
                      iter);
            iter += class_count_; 
        }
        std::copy(classes.begin(),
                  classes.end(),
                  iter);
        #undef PUSH
    }

    void make_from_map(map_type & in) // -> const: .operator[] -> .find
    {
        #define PULL(item_, type_) item_ = type_(in[#item_][0]); 
        PULL(column_count_,int);
        PULL(class_count_, int);
        PULL(row_count_, int);
        PULL(actual_mtry_,int);
        PULL(actual_msample_, int);
        PULL(problem_type_, (Problem_t)int);
        PULL(is_weighted_, int);
        PULL(used_, int);
        PULL(precision_, double);
        PULL(response_size_, int);
        class_weights_ = in["class_weights_"];
        #undef PULL
    }
    void make_map(map_type & in) const
    {
        #define PUSH(item_) in[#item_] = double_array(1, double(item_));
        PUSH(column_count_);
        PUSH(class_count_)
        PUSH(row_count_);
        PUSH(actual_mtry_);
        PUSH(actual_msample_);
        PUSH(problem_type_);
        PUSH(is_weighted_);
        PUSH(used_);
        PUSH(precision_);
        PUSH(response_size_);
        in["class_weights_"] = class_weights_;
        #undef PUSH
    }
    
    /**\brief set default values (-> values not set)
     */
    ProblemSpec()
    :   column_count_(0),
        class_count_(0),
        row_count_(0),
        actual_mtry_(0),
        actual_msample_(0),
        problem_type_(CHECKLATER),
        used_(false),
        is_weighted_(false),
        precision_(0.0),
        response_size_(1)
    {}


    ProblemSpec & column_count(int in)
    {
        column_count_ = in;
        return *this;
    }

    /**\brief supply with class labels -
     * 
     * the preprocessor will not calculate the labels needed in this case.
     */
    template<class C_Iter>
    ProblemSpec & classes_(C_Iter begin, C_Iter end)
    {
        classes.clear();
        int size = end-begin;
        for(int k=0; k<size; ++k, ++begin)
            classes.push_back(detail::RequiresExplicitCast<LabelType>::cast(*begin));
        class_count_ = size;
        return *this;
    }

    /** \brief supply with class weights  -
     *
     * this is the only case where you would really have to 
     * create a ProblemSpec object.
     */
    template<class W_Iter>
    ProblemSpec & class_weights(W_Iter begin, W_Iter end)
    {
        class_weights_.clear();
        class_weights_.insert(class_weights_.end(), begin, end);
        is_weighted_ = true;
        return *this;
    }



    void clear()
    {
        used_ = false; 
        classes.clear();
        class_weights_.clear();
        column_count_ = 0 ;
        class_count_ = 0;
        actual_mtry_ = 0;
        actual_msample_ = 0;
        problem_type_ = CHECKLATER;
        is_weighted_ = false;
        precision_   = 0.0;
        response_size_ = 0;

    }

    bool used() const
    {
        return used_ != 0;
    }
};


//@}



/**\brief Standard early stopping criterion
 *
 * Stop if region.size() < min_split_node_size_;
 */
class EarlyStoppStd
{
    public:
    int min_split_node_size_;
    int max_tree_depth_;

    template<class Opt>
    EarlyStoppStd(Opt opt)
    :   min_split_node_size_(opt.min_split_node_size_),
        max_tree_depth_(opt.max_tree_depth_ == 0 ?
                            std::numeric_limits<int>::max()
                            : opt.max_tree_depth_)
    {}

    template<class T>
    void set_external_parameters(ProblemSpec<T>const  &, int /* tree_count */ = 0, bool /* is_weighted_ */ = false)
    {}

    template<class Region>
    bool operator()(Region& region)
    {
        return (region.size() < min_split_node_size_)||(region.depth() >= max_tree_depth_);
    }

    template<class WeightIter, class T, class C>
    bool after_prediction(WeightIter,  int /* k */, MultiArrayView<2, T, C> /* prob */, double /* totalCt */)
    {
        return false; 
    }
};


} // namespace vigra

#endif //VIGRA_RF_COMMON_HXX
