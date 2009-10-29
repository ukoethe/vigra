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


#ifndef VIGRA_RF_COMMON_HXX
#define VIGRA_RF_COMMON_HXX

namespace vigra
{

// FORWARD DECLARATIONS
// TODO : DECIDE WHETHER THIS IS A GOOD IDEA
struct                  ClassificationTag{};
struct                  RegressionTag{};

class                   GiniSplit;
class                   StopVisiting;
template<class T>
    class               OOB_Visitor;
class                   RandomForestOptions;
class                   ProblemSpec;
template<class Tag>
    class               RandomForest;
class 					EarlyStoppStd;
namespace detail
{
	class              	RF_DEFAULT;
	class               DecisionTree;
}
detail::RF_DEFAULT&     rf_default();
template <class T>
    class               DT_Region;

/**\brief Traits Class for the Random Forest
 *
 * refer to the typedefs in this class when using the default
 * objects as the names may change in future.
 */
class RF_Traits
{
public:
	typedef RandomForestOptions 	Options_t;
	typedef	detail::DecisionTree   	DecisionTree_t;
	typedef ProblemSpec				ProblemSpec_t;
	typedef ClassificationTag		Preprocessor_t;
	typedef GiniSplit	            Default_Split_t;
	typedef EarlyStoppStd           Default_Stop_t;
	typedef OOB_Visitor<StopVisiting>
									Default_Visitor_t;
	typedef StopVisiting			StopVisiting_t;

};

/**\brief Standard early stopping criterion
 *
 * Stop if region.size() < min_split_node_size_;
 */
class EarlyStoppStd
{
    public:
	int min_split_node_size_;

	template<class Opt>
	EarlyStoppStd(Opt opt)
	:	min_split_node_size_(opt.min_split_node_size_)
	{}

	void set_external_parameters(RF_Traits::ProblemSpec_t &prob)
	{}

    template<class Region>
    bool operator()(Region& region)
    {
        return region.size() < min_split_node_size_;
    }
};




namespace detail
{

/**\brief singleton default tag class -
 *
 *  use the rf_default() factory function to use the tag.
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
    	typedef DT_Region<ArrayVectorView<Int32>::iterator>
					StackEntry_t;
};

/**\brief chooses between default type and type supplied
 * 
 * used in RandomForest.learn()
 * Usage:
 *\code
 * 		// example: use container type supplied by user or ArrayVector if 
 * 		// 			rf_default() was specified as argument;
 *		template<class Container_t>
 *		void do_some_foo(Container_t in)
 *		{
 *			typedef ArrayVector<int> 	Default_Container_t;
 *			Default_Container_t			default_value;
 *			Value_Chooser<Container_t, 	Default_Container_t> 
 *						choose(in, default_value);
 *
 *			// if the user didn't care and the in was of type 
 *			// RF_DEFAULT then default_value is used.
 *			do_some_more_foo(choose.value());
 *		}
 *		Value_Chooser choose_val<Type, Default_Type>
 *\endcode
 */
template<class T, class C>
class Value_Chooser
{
    public:
    T & value_;
    Value_Chooser(T &val, C  default_val)
    :   value_(val)
    {}
    T& value()
    {
        return value_;
    }
};

template<class C>
class Value_Chooser<detail::RF_DEFAULT, C>
{
    public:
  C  value_;
  Value_Chooser(detail::RF_DEFAULT val, C default_val)
  :     value_(default_val)
  {}

  C& value()
  {
      return value_;
  }
};




} //namespace detail


/**\brief factory function to return a RF_DEFAULT tag
 *
 */
detail::RF_DEFAULT& rf_default()
{
	static detail::RF_DEFAULT result;
	return result;
};

/** tags used with the RandomForestOptions class
 *
 */
enum RF_OptionTag   { RF_EQUAL,
					  RF_PROPORTIONAL,
					  RF_EXTERNAL,
					  RF_NONE,
					  RF_FUNCTION,
					  RF_LOG,
					  RF_SQRT,
					  RF_CONST};


/**\brief Options object for the random forest
 *
 * usage:
 * RandomForestOptions a =  RandomForestOptions()
 *                          	.param1(value1)
 *								.param2(value2)
 *								...
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
	double 	training_set_proportion_;
	int 	training_set_size_;
   	int (*training_set_func_)(int);
	RF_OptionTag
		training_set_calc_switch_;

    bool 	sample_with_replacement_;
    RF_OptionTag
	  		stratification_method_;


	/**\name general random forest options
	 *
	 * these usually will be used by most split functors and
	 * stopping predicates
	 */
	/*\{*/
	RF_OptionTag 	mtry_switch_;
	int 	mtry_;
   	int (*mtry_func_)(int) ;


	int tree_count_;
	int min_split_node_size_;
	/*\}*/

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
		tree_count_(256),
		min_split_node_size_(1)
	{}

	/**\brief specify stratification strategy
	 *
	 * default: RF_NONE
	 * possible values: RF_EQUAL, RF_PROPORTIONAL,
	 * 					RF_EXTERNAL, RF_NONE
	 * RF_EQUAL: 		get equal amount of samples per class.
	 * RF_PROPORTIONAL:	sample proportional to fraction of class samples
	 * 					in population
	 * RF_EXTERNAL:	    strata_weights_ field of the ProblemSpec_t object
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
	 * 		  tree should be learnt with.
	 *
	 * \param in function pointer that takes the number of rows in the
	 * 			 learning data and outputs the number samples per tree.
	 */
	RandomForestOptions & samples_per_tree(int (*in)(int))
	{
		training_set_func_ = in;
		training_set_calc_switch_ = RF_FUNCTION;
		return *this;
	}

	/**\brief use built in mapping to calculate mtry
	 *
	 * Use one of the built in mappings to calculate mtry from the number
	 * of columns in the input feature data.
	 * \param in possible values: RF_LOG or RF_SQRT
	 * 			 <br> default: RF_SQRT.
	 */
	RandomForestOptions & features_per_node(RF_OptionTag in)
	{
		vigra_precondition(in == RF_LOG ||
						   in == RF_SQRT,
						   "RandomForestOptions()::features_per_node():"
						   "input must be of type RF_LOG or RF_SQRT");
		mtry_switch_ = in;
		return *this;
	}

	/**\brief Set mtry to a constant value
	 *
	 * mtry is the number of columns/variates/variables randomly choosen
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
	 * 			 of the and outputs int (mtry)
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
	 * 	When the number of examples in a node is below this number,
	 * 	the node is not	split even if class separation is not yet perfect.
	 * 	Instead, the node returns the proportion of each class
	 * 	(among the remaining examples) during the prediction phase.
	 * 	<br> Default: 1 (complete growing)
	 */
	RandomForestOptions & min_split_node_size(int in)
	{
		min_split_node_size_ = in;
		return *this;
	}
};


/** \brief problem specification class for the random forest.
 *
 * This class contains all the problem specific parameters the random
 * forest needs for learning. Specification of an instance of this class
 * is optional as all necessary fields will be computed prior to learning
 * if not specified.
 *
 * if needed usage is similar to that of RandomForestOptions
 */

class ProblemSpec
{


public:

	/** \brief 	problem class
	 */
	enum Problem_t{REGRESSION, CLASSIFICATION, CHECKLATER};

	/** \brief allowed types
	 *
	 * depending on the class label type the class labels
	 * will be internally stored as double, Int64 or std::string
	 */

	enum Types_t
				{UInt8_t, UInt16_t, UInt32_t, UInt64_t,
			   	 Int8_t, Int16_t, Int32_t, Int64_t,
			     double_t, float_t, UNKNOWN};

private:
   // dirty little helper functions -
   // I really see no reason to keep the class label type as a template
   // parameter just because of the predict function.
   // As numClasses << 100 this solution will not cause much overhead I
   // guess.
   //
   // todo it would be an alternative to just save all class types as
   // strings or doubles - whatever.
#define make_dlh(typee_)\
private:\
	ArrayVector<typee_>			typee_##_classes_;\
	Types_t type_of(typee_ in)\
	{\
		return typee_##_t;\
	}\
public:\
	void to_classlabel(int index, typee_  & out) const\
	{\
		out = typee_(typee_##_classes_[index]);\
	}\
private:

	make_dlh(UInt8);
	make_dlh(UInt16);
	make_dlh(UInt32);
	make_dlh(UInt64);
	make_dlh(Int8);
	make_dlh(Int16);
	make_dlh(Int32);
	make_dlh(Int64);
	make_dlh(double);
	make_dlh(float);

#undef make_dlh

public:

	int 		column_count_;
	int 		class_count_;
    int      row_count_;

	int 		actual_mtry_;
	int 		actual_msample_;

	Problem_t 	problem_type_;
	Types_t 	class_type_;



	ArrayVector<double> 		class_weights_;
	bool is_weighted;

	/**\brief set default values (-> values not set)
	 */
	ProblemSpec()
	:	column_count_(0),
		class_count_(0),
		row_count_(0),
		actual_mtry_(0),
		actual_msample_(0),
		problem_type_(CHECKLATER),
		class_type_(UNKNOWN),
		is_weighted(false)
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
		for(C_Iter iter = begin; iter!= end; ++iter)
		{
			UInt8_classes_.push_back(UInt8(*iter));
			UInt16_classes_.push_back(UInt16(*iter));
			UInt32_classes_.push_back(UInt32(*iter));
			UInt64_classes_.push_back(UInt64(*iter));
			Int8_classes_.push_back(Int8(*iter));
			Int16_classes_.push_back(Int16(*iter));
			Int32_classes_.push_back(Int32(*iter));
			Int64_classes_.push_back(Int64(*iter));
			double_classes_.push_back(double(*iter));
			float_classes_.push_back(float(*iter));
		}
		class_type_ = type_of(*begin);
		class_count_ = end-begin;
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
		for(W_Iter iter = begin; iter!= end; ++iter)
		{
			class_weights_.push_back(*iter);
		}
		is_weighted = true;
		return *this;
	}



	void clear()
	{
		used_ = false; 
		UInt8_classes_.clear();
		UInt16_classes_.clear();
		UInt32_classes_.clear();
		UInt64_classes_.clear();
		Int8_classes_.clear();
		Int16_classes_.clear();
		Int32_classes_.clear();
		Int64_classes_.clear();
		double_classes_.clear();
		float_classes_.clear();
		class_weights_.clear();
		column_count_ = 0 ;
		class_count_ = 0;
		actual_mtry_ = 0;
		actual_msample_ = 0;
		problem_type_ = CHECKLATER;
		class_type_ = UNKNOWN;
		is_weighted = false;

	}

	bool used_;
	bool used() const
	{
    	return used_;
	}
};



} // namespace vigra

#endif //VIGRA_RF_COMMON_HXX
