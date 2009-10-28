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

#ifndef VIGRA_RF_PREPROCESSING_HXX
#define VIGRA_RF_PREPROCESSING_HXX


#include "rf_common.hxx"

namespace vigra
{


template<class Tag, class T1, class C1, class T2, class C2>
class Processor
{};

namespace detail
{
	void fill_external_parameters(RF_Traits::Options_t & options,
								  RF_Traits::ProblemSpec_t & ext_param)
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
                    int(1+(std::log(double(ext_param.column_count_))/std::log(2.0)));
                break;
            case RF_FUNCTION:
                ext_param.actual_mtry_ =
                    options.mtry_func_(ext_param.column_count_);
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
}

template<class T1, class C1, class T2, class C2>
class Processor<ClassificationTag, T1, C1, T2, C2>
{
    public:
	typedef Int32 LabelInt;
	MultiArrayView<2, T1, C1>const & 	features_;
	MultiArray<2, LabelInt> 	intLabels_;
	MultiArrayView<2, LabelInt> 	strata_;

    Processor(MultiArrayView<2, T1, C1>const & features,   MultiArrayView<2, T2, C2>const & response,
              RF_Traits::Options_t &options,         RF_Traits::ProblemSpec_t &ext_param)
    :features_( features)
    {
        ext_param.column_count_  = features.shape(1);
        ext_param.row_count_     = features.shape(0);
        ext_param.problem_type_  = RF_Traits::ProblemSpec_t::CLASSIFICATION;
        ext_param.used_          = true;
        intLabels_.reshape(response.shape());

        //get the class labels
        if(ext_param.class_count_ == 0)
        {
            std::map<T2, LabelInt > labelToInt;
            ArrayVector<T2>         classes_;
            for(MultiArrayIndex k = 0; k < features.shape(0); ++k)
            {
                typename std::map<T2, LabelInt >::iterator lit = labelToInt.find(response(k, 0));
                if(lit == labelToInt.end())
                {
                    intLabels_(k, 0) = labelToInt[response(k,0)] = classes_.size();
                    classes_.push_back(response(k,0));
                }
                else
                {
                    intLabels_(k, 0) = lit->second;
                }
            }
            ext_param.classes_(classes_.begin(), classes_.end());
        }

        // set class weights
        if(ext_param.class_weights_.size() == 0)
        {
            ArrayVector<T2> tmp((std::size_t)ext_param.class_count_, NumericTraits<T2>::one());
            ext_param.class_weights(tmp.begin(), tmp.end());
        }

        // set mtry and msample
        detail::fill_external_parameters(options, ext_param);

        // set strata
        strata_ = intLabels_;

    }

	MultiArrayView<2, T1, C1>const & features()
	{
		return features_;
	}

	MultiArrayView<2, LabelInt>& response()
	{
		return intLabels_;
	}

	MultiArrayView<2, LabelInt>&  strata()
	{
		return intLabels_;
	}

	ArrayVectorView< double> strata_prob()
	{
		return ArrayVectorView< double>();
	}
};




template<class T1, class C1, class T2, class C2>
class Processor<RegressionTag, T1, C1, T2, C2>
{
public:
	typedef T1 		Feature_t;
	typedef T2 		Response_t;

	MultiArrayView<2, T1, C1> 	features_;
	MultiArrayView<2, T2, C2> 	response_;
	RF_Traits::Options_t const & options_;
	RF_Traits::ProblemSpec_t const &
								ext_param_;
	// will only be filled if needed
	MultiArray<2, int> 	 	strata_;
	bool strata_filled;

	Processor(	MultiArrayView<2, T1, C1> 	feats,
				MultiArrayView<2, T2, C2> 	response,
				RF_Traits::Options_t			options,
				RF_Traits::ProblemSpec_t	ext_param)
	:
		features_(feats),
		response_(response),
		options_(options),
		ext_param_(ext_param)
	{
		detail::fill_external_parameters(options, ext_param);
	}

	MultiArrayView<2, T1, C1> & features()
	{
		return features_;
	}

	MultiArrayView<2, T2, C2> & response()
	{
		return response_;
	}

	MultiArrayView<2, int> & strata()
	{
		if(strata_filled)
		{
			return strata_;
		}
		else
		{
			//
		}
	}
};
}
#endif //VIGRA_RF_PREPROCESSING_HXX



