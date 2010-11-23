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

#ifndef VIGRA_RF_PREPROCESSING_HXX
#define VIGRA_RF_PREPROCESSING_HXX

#include <limits>
#include "rf_common.hxx"

namespace vigra
{

/** Class used while preprocessing  (currently used only during learn)
 *
 * This class is internally used by the Random Forest learn function. 
 * Different split functors may need to process the data in different manners
 * (i.e., regression labels that should not be touched and classification 
 * labels that must be converted into a integral format)
 *
 * This Class only exists in specialized versions, where the Tag class is 
 * fixed. 
 *
 * The Tag class is determined by Splitfunctor::Preprocessor_t . Currently
 * it can either be ClassificationTag or RegressionTag.  look At the 
 * RegressionTag specialisation for the basic interface if you ever happen
 * to care.... - or need some sort of vague new preprocessor.  
 * new preprocessor ( Soft labels or whatever)
 */
template<class Tag, class LabelType, class T1, class C1, class T2, class C2>
class Processor;

namespace detail
{

	/** Common helper function used in all Processors. 
	 * This function analyses the options struct and calculates the real 
	 * values needed for the current problem (data)
	 */
	template<class T>
	void fill_external_parameters(RandomForestOptions & options,
								  ProblemSpec<T> & ext_param)
	{
	    // set correct value for mtry
	    switch(options.mtry_switch_)
	    {
	        case RF_SQRT:
                ext_param.actual_mtry_ =
                    int(std::floor(
                            std::sqrt(double(ext_param.column_count_))
                            + 0.5));
                break;
            case RF_LOG:
                // this is in Breimans original paper
                ext_param.actual_mtry_ =
                    int(1+(std::log(double(ext_param.column_count_))
						   /std::log(2.0)));
                break;
            case RF_FUNCTION:
                ext_param.actual_mtry_ =
                    options.mtry_func_(ext_param.column_count_);
                break;
			case RF_ALL:
				ext_param.actual_mtry_ = ext_param.column_count_;
				break;
            default:
                ext_param.actual_mtry_ =
                    options.mtry_;
	    }
	    // set correct value for msample
	    switch(options.training_set_calc_switch_)
	    {
	        case RF_CONST:
                ext_param.actual_msample_ =
                    options.training_set_size_;
                break;
            case RF_PROPORTIONAL:
                ext_param.actual_msample_ =
                    (int)std::ceil(  options.training_set_proportion_ *
                                     ext_param.row_count_);
                    break;
            case RF_FUNCTION:
                ext_param.actual_msample_ =
                    options.training_set_func_(ext_param.row_count_);
                break;
            default:
                vigra_precondition(1!= 1, "unexpected error");

	    }

	}
	
	/** Returns true if MultiArray contains NaNs
	 */
	template<unsigned int N, class T, class C>
	bool contains_nan(MultiArrayView<N, T, C> const & in)
	{
		for(int ii = 0; ii < in.size(); ++ii)
			if(in[ii] != in[ii])
				return true;
		return false; 
	}
	
	/** Returns true if MultiArray contains Infs
	 */
	template<unsigned int N, class T, class C>
	bool contains_inf(MultiArrayView<N, T, C> const & in)
	{
		 if(!std::numeric_limits<T>::has_infinity)
			 return false;
		 for(int ii = 0; ii < in.size(); ++ii)
			if(in[ii] == std::numeric_limits<T>::infinity())
				return true;
		 return false;
	}
}



/** Preprocessor used during Classification
 *
 * This class converts the labels int Integral labels which are used by the 
 * standard split functor to address memory in the node objects.
 */
template<class LabelType, class T1, class C1, class T2, class C2>
class Processor<ClassificationTag, LabelType, T1, C1, T2, C2>
{
    public:
	typedef Int32 LabelInt;
	typedef	MultiArrayView<2, T1, C1> Feature_t;
	typedef MultiArray<2, T1> FeatureWithMemory_t;
	typedef MultiArrayView<2,LabelInt> Label_t;
	MultiArrayView<2, T1, C1>const & 	features_;
	MultiArray<2, LabelInt> 			intLabels_;
	MultiArrayView<2, LabelInt> 		strata_;

	template<class T>
    Processor(MultiArrayView<2, T1, C1>const & features,   
			  MultiArrayView<2, T2, C2>const & response,
              RandomForestOptions &options,         
			  ProblemSpec<T> &ext_param)
    :
		features_( features) // do not touch the features. 
    {
		vigra_precondition(!detail::contains_nan(features), "Processor(): Feature Matrix "
						   								   "Contains NaNs");
		vigra_precondition(!detail::contains_nan(response), "Processor(): Response "
						   								   "Contains NaNs");
		vigra_precondition(!detail::contains_inf(features), "Processor(): Feature Matrix "
						   								   "Contains inf");
		vigra_precondition(!detail::contains_inf(response), "Processor(): Response "
						   								   "Contains inf");
		// set some of the problem specific parameters 
        ext_param.column_count_  = features.shape(1);
        ext_param.row_count_     = features.shape(0);
        ext_param.problem_type_  = CLASSIFICATION;
        ext_param.used_          = true;
        intLabels_.reshape(response.shape());

        //get the class labels
        if(ext_param.class_count_ == 0)
        {
			// fill up a map with the current labels and then create the 
			// integral labels.
            std::set<T2>             		labelToInt;
            for(MultiArrayIndex k = 0; k < features.shape(0); ++k)
                labelToInt.insert(response(k,0));
            std::vector<T2> tmp_(labelToInt.begin(), labelToInt.end());
            ext_param.classes_(tmp_.begin(), tmp_.end());
        }
        for(MultiArrayIndex k = 0; k < features.shape(0); ++k)
        {
            if(std::find(ext_param.classes.begin(), ext_param.classes.end(), response(k,0)) == ext_param.classes.end())
            {
                throw std::runtime_error("unknown label type");
            }
            else
                intLabels_(k, 0) = std::find(ext_param.classes.begin(), ext_param.classes.end(), response(k,0))
                                    - ext_param.classes.begin();
        }
        // set class weights
        if(ext_param.class_weights_.size() == 0)
        {
            ArrayVector<T2> 
				tmp((std::size_t)ext_param.class_count_, 
					NumericTraits<T2>::one());
            ext_param.class_weights(tmp.begin(), tmp.end());
        }

        // set mtry and msample
        detail::fill_external_parameters(options, ext_param);

        // set strata
        strata_ = intLabels_;

    }

	/** Access the processed features
	 */
	MultiArrayView<2, T1, C1>const & features()
	{
		return features_;
	}

	/** Access processed labels
	 */
	MultiArrayView<2, LabelInt>& response()
	{
		return intLabels_;
	}

	/** Access processed strata
	 */
	ArrayVectorView < LabelInt>  strata()
	{
		return ArrayVectorView<LabelInt>(intLabels_.size(), intLabels_.data());
	}

	/** Access strata fraction sized - not used currently
	 */
	ArrayVectorView< double> strata_prob()
	{
		return ArrayVectorView< double>();
	}
};



/** Regression Preprocessor - This basically does not do anything with the
 * data.
 */
template<class LabelType, class T1, class C1, class T2, class C2>
class Processor<RegressionTag,LabelType, T1, C1, T2, C2>
{
public:
	// only views are created - no data copied.
	MultiArrayView<2, T1, C1> 	features_;
	MultiArrayView<2, T2, C2> 	response_;
	RandomForestOptions const & options_;
	ProblemSpec<LabelType> const &
								ext_param_;
	// will only be filled if needed
	MultiArray<2, int> 	 	strata_;
	bool strata_filled;

	// copy the views.
	template<class T>
	Processor(	MultiArrayView<2, T1, C1> 	features,
				MultiArrayView<2, T2, C2> 	response,
				RandomForestOptions const &	options,
				ProblemSpec<T>	ext_param)
	:
		features_(features),
		response_(response),
		options_(options),
		ext_param_(ext_param)
	{
		detail::fill_external_parameters(options, ext_param);
		vigra_precondition(!detail::contains_nan(features), "Processor(): Feature Matrix "
						   								   "Contains NaNs");
		vigra_precondition(!detail::contains_nan(response), "Processor(): Response "
						   								   "Contains NaNs");
		vigra_precondition(!detail::contains_inf(features), "Processor(): Feature Matrix "
						   								   "Contains inf");
		vigra_precondition(!detail::contains_inf(response), "Processor(): Response "
						   								   "Contains inf");
		strata_ = MultiArray<2, int> (MultiArrayShape<2>::type(response_.shape(0), 1));
	}

	/** access preprocessed features
	 */
	MultiArrayView<2, T1, C1> & features()
	{
		return features_;
	}

	/** access preprocessed response
	 */
	MultiArrayView<2, T2, C2> & response()
	{
		return response_;
	}

	/** access strata - this is not used currently
	 */
	MultiArrayView<2, int> & strata()
	{
		return strata_;
	}
};
}
#endif //VIGRA_RF_PREPROCESSING_HXX



