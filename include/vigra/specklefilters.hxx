/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2014 by Benjamin Seppke                 */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_SPECKLEFILTER_HXX
#define VIGRA_SPECKLEFILTER_HXX

#include "basicimage.hxx"
#include "inspectimage.hxx"

#include "applywindowfunction.hxx"

namespace vigra {

namespace detail 
{

/**
 * Helper to store distances in lookuptable (LUT)    
 */
vigra::FImage distanceLUT(vigra::Diff2D const & window_shape)
{
    vigra::FImage res(window_shape);
    
    int y, x;
    double  w_half = window_shape.x/2.0,
            h_half = window_shape.y/2.0,
            x_diff, y_diff;
    
    for(y=0; y != window_shape.y; y++)
    {	
        for(x=0; x != window_shape.x; x++)
        {
            x_diff = x-w_half;
            y_diff = y-h_half;
            res(x,y) = sqrt(x_diff*x_diff + y_diff*y_diff);
        }
    }
    
    return res;
}

} //end namespace detail

/*********************************************************************
 *                                                                   *
 * The (Basic) Frost Filter                                          *
 *                                                                   *
 *     Parameters:		window_shape   The size of the filter        *
 *                      k              The damping factor (0,...,1)  *
 *********************************************************************/
 
/**  
	This function tries to reduce the speckle noise of an image by means of applying the 
	basic Frost filter using a window of given size and a damping factor k. The implementation
	is according to the article by 
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
*/
//@{

/** \brief 	This function tries to reduce the speckle noise of an image by applying the basic Frost filter.

	The user has to provide a window size and a damping factor k. The implementation is
	according to the article by  
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
    
    All restrictions of the called functions \ref applyWindowFunction apply.
    
    <b> Preconditions:</b>
    \code  
    0.0 < k <= 1.0
    \endcode
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        frostFilter(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> dest,
                    Diff2D window_shape, float k,
                    BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);

    }
    \endcode

    \deprecatedAPI{frostFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void frostFilter(SrcIterator supperleft,
                         SrcIterator slowerright, SrcAccessor sa,
                         DestIterator dupperleft, DestAccessor da,
                     	 Diff2D window_shape, float k,
                         BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        frostFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    Diff2D window_shape, float k,
                    BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/specklefilters.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int w=1000, h=1000;
    MultiArray<2, float> src(w,h), dest(w,h);
    ...
	
    // apply a (basic) frost filter with a window size of 5x5 and a damping factor of 0.5
    frostFilter(src, dest, Diff2D(5,5), 0.5);
    \endcode
*/

doxygen_overloaded_function(template <...> void frostFilter)
    
template<typename VALUETYPE>
class FrostFunctor
{
public:
	FrostFunctor(Diff2D window_shape, float k)
	: m_window_shape(window_shape),
	  m_k(k),
      m_dist(detail::distanceLUT(window_shape))
	{
		using namespace vigra;		
		vigra_precondition( k>0 && k<=1 , "vigra::FrostFunctor(): Damping factor k has to be: 0 < k <= 1!");
	}
	
	template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
	void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc)
	{		
		using namespace vigra;
		
		SrcIterator s_ul = s - m_window_shape/2,
					s_lr = s_ul + m_window_shape;
		
		FindAverageAndVariance<VALUETYPE> averageAndVariance;   // init functor
		
		inspectImage(s_ul, s_lr, s_acc, averageAndVariance);
		
		/*As defined in: Lopez & Touzi & Nezry: Adaptive speckle filters and scene heterogenity*/
		VALUETYPE	C_I2 = averageAndVariance.variance() / (averageAndVariance.average() * averageAndVariance.average()),
					sum_m = 0.0,
					sum_pm = 0.0,
					m = 0.0,
					dist = 0.0,
					p = 0.0;
		
		SrcIterator ys = s_ul;
		SrcIterator xs = ys;
        
        FImage::const_traverser ydist = m_dist.upperLeft();
        FImage::const_traverser xdist = ydist;
        FImage::Accessor   dist_acc = m_dist.accessor();
		
		//convolve mask with each impulse response to compute the result of the frost filter
		int y, x;
		for(y=0 ; ys.y!=s_lr.y; ys.y++, ydist.y++, y++)
		{	
			for(xs=ys, xdist=ydist, x=0; xs.x!=s_lr.x; xs.x++, xdist.x++, x++)
			{
				p = s_acc(xs);
				
				//impuls response of the frost filter
				m = exp(-1 * m_k * C_I2 * dist_acc(xdist));
				
				//convolve
				sum_pm += m * p;
				sum_m += m;
			}
		}
		//normalize
		d_acc.set(sum_pm/sum_m, d);
	}
	
	Diff2D windowShape() const
	{
		return m_window_shape;
	}
	
private:
	Diff2D m_window_shape;
	float m_k;	
    FImage m_dist;
};

template <class SrcIterator, class SrcAccessor, 
		  class DestIterator, class DestAccessor>
inline void frostFilter(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                 		DestIterator d_ul, DestAccessor d_acc, 
                 		Diff2D window_shape, float k,
                 		BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    FrostFunctor<typename SrcIterator::value_type> func(window_shape, k);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class SrcIterator, class SrcAccessor, 
class DestIterator, class DestAccessor>
inline void frostFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                 		pair<DestIterator, DestAccessor> d, 
                 		Diff2D window_shape, float k,
                 		BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    frostFilter(s.first, s.second, s.third,
                d.first, d.second, 
                window_shape, k,
                border);
}


template <class T1, class S1, 
          class T2, class S2>
inline void frostFilter(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, T2, S2> dest,
                 		Diff2D window_shape, float k,
                 		BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
	vigra_precondition(src.shape() == dest.shape(),
						"vigra::frostFilter(): Shape mismatch between input and output.");
	frostFilter(srcImageRange(src),
				destImage(dest),  
                window_shape, k,
                border);
}


/*********************************************************************
 *                                                                   *
 * The Enhanced Frost Filter                                         *
 *                                                                   *
 *     Parameters:		window_shape   The size of the filter        *
 *                      k              The damping factor (0,...,1)  *
 *                      enl            Eq. Num. Looks for comp. of   *
 *                                     the thresholds C_u and C_max  *
 *********************************************************************/
 
/**  
	This function tries to reduce the speckle noise of an image by means of applying the 
	enhanced Frost filter using a window of given size, a damping factor k, and the equivalent
	numbers of look (enl). The implementation is according to the article by 
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
*/

/** \brief 	This function tries to reduce the speckle noise of an image by applying the Enhanced Frost filter.

	The user has to provide a window size, a damping factor k, and the equivalent
	numbers of look (enl). The implementation is according to the article by  
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
    
    All restrictions of the called functions \ref applyWindowFunction apply.
    
    <b> Preconditions:</b>
    \code  
	1. 0.0 < k <= 1.0
	2. enl > 0
    \endcode
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        enhancedFrostFilter(MultiArrayView<2, T1, S1> const & src,
                   			 MultiArrayView<2, T2, S2> dest,
                   			 Diff2D window_shape, float k, int enl,
                    		 BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);

    }
    \endcode

    \deprecatedAPI{enhancedFrostFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void enhancedFrostFilter(SrcIterator supperleft,
                       			  SrcIterator slowerright, SrcAccessor sa,
                       			  DestIterator dupperleft, DestAccessor da,
                     			  Diff2D window_shape, float k, int enl,
                       			  BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        enhancedFrostFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   			 pair<DestIterator, DestAccessor> dest,
                   			 Diff2D window_shape, float k, int enl,
                   			 BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/specklefilters.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int w=1000, h=1000;
    MultiArray<2, float> src(w,h), dest(w,h);
    ...
	
    // apply an enhanced frost filter with a window size of 5x5 and a damping factor of 0.5, where
    // the image was composed by 3 equivalent looks:
    enhancedFrostFilter(src, dest, Diff2D(5,5), 0.5, 3);
    \endcode
*/

template<typename VALUETYPE>
class EnhancedFrostFunctor
{
public:
	EnhancedFrostFunctor(Diff2D window_shape, float k, int enl)
	: m_window_shape(window_shape),
	  m_k(k),
    m_enl(enl),
    m_dist(detail::distanceLUT(window_shape))

	{
		using namespace vigra;
		vigra_precondition( k>0 && k<=1 , "vigra::EnhancedFrostFunctor(): Damping factor k has to be: 0 < k <= 1!");
		vigra_precondition( enl>0, "vigra::EnhancedFrostFunctor(): Equivalent number of looks (enl) must be larger than zero!");
	}
	
	template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
	void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc)
	{
		using namespace vigra;
		
		SrcIterator s_ul = s - m_window_shape/2,
                    s_lr = s_ul + m_window_shape;
		
		FindAverageAndVariance<VALUETYPE> averageAndVariance;   // init functor
		
		inspectImage(s_ul, s_lr, s_acc, averageAndVariance);
		
		/*As defined in: Lopez & Touzi & Nezry: Adaptive speckle filters and scene heterogenity*/
		/* With ENL -> C_u and ENL -> C_max from ENVI: online_help/Using_Adaptive_Filters.html */
		VALUETYPE		C_u    = 0.523/sqrt((double)m_enl),
						C_max  = sqrt(1+2.0/m_enl),
						C_I = sqrt(averageAndVariance.variance()) / averageAndVariance.average(),
						sum_m = 0.0,
						sum_pm = 0.0,
						m = 0.0,
						dist = 0.0,
						p = 0.0;
		
		SrcIterator ys = s_ul;
		SrcIterator xs = ys;
        
        FImage::const_traverser ydist = m_dist.upperLeft();
        FImage::const_traverser xdist = ydist;
        FImage::Accessor   dist_acc = m_dist.accessor();
		
		//convolve mask with each impulse response to compute the result of the frost filter
		int y, x;
		for(y=0 ; ys.y!=s_lr.y; ys.y++, ydist.y++, y++)
		{	
			for(xs=ys, xdist=ydist, x=0; xs.x!=s_lr.x; xs.x++, xdist.x++, x++)
			{
				p = s_acc(xs);
				
				//impuls response of the frost filter
				m = exp(-m_k * func(C_I, C_max, C_u) * dist_acc(xdist));
				
				//convolve
				sum_pm += m * p;
				sum_m += m;
			}
		}
		//normalize
		d_acc.set(sum_pm/sum_m, d);
	}
	
	Diff2D windowShape() const
	{
		return m_window_shape;
	}
	
private:
	//The penalisier function:
	//As defined in: Shi & Fung: A comparison of Digital Speckle Filters
	inline double func(double C_I, double C_max, double C_u) const
	{
		if(C_I < C_u)
		{
			return 0;
		}
		else if (C_I <= C_max)
		{
			return (C_I - C_u)/(C_max - C_I);
		}
		else
		{
			return 1.0e100;
		}
	}
	
	Diff2D m_window_shape;
	float m_k;	
	int m_enl;
    FImage m_dist;
};
    
template <class SrcIterator, class SrcAccessor, 
		  class DestIterator, class DestAccessor>
inline void enhancedFrostFilter(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                 				DestIterator d_ul, DestAccessor d_acc, 
                 				Diff2D window_shape, float k, int enl,
                 				BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    EnhancedFrostFunctor<typename SrcIterator::value_type> func(window_shape, k, enl);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class SrcIterator, class SrcAccessor, 
		  class DestIterator, class DestAccessor>
inline void enhancedFrostFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                         		pair<DestIterator, DestAccessor> d, 
                         		Diff2D window_shape, float k, int enl,
                        		BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    enhancedFrostFilter(s.first, s.second, s.third,
                        d.first, d.second, 
                        window_shape, k, enl,
                        border);
}


template <class T1, class S1, 
          class T2, class S2>
inline void enhancedFrostFilter(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, T2, S2> dest,
                 		Diff2D window_shape, float k, int enl,
                 		BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
	vigra_precondition(src.shape() == dest.shape(),
						"vigra::enhancedFrostFilter(): Shape mismatch between input and output.");
	enhancedFrostFilter(srcImageRange(src),
						destImage(dest),  
                		window_shape, k, enl,
                		border);
}



/*********************************************************************
 *                                                                   *
 * The Gamma Maximum A Posteriori (MAP) Filter                       *
 *                                                                   *
 *     Parameters:		window_shape   The size of the filter        *
 *               		enl            Eq. Num. Looks for comp. of   *
 *                                     the thresholds C_u and C_max  *
 *********************************************************************/

/**  
	This function tries to reduce the speckle noise of an image by means of applying the 
	Gamma Maximum A Posteriori (MAP) filter using a window of given size, and the equivalent
	numbers of look (enl). The implementation is according to the article by 
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
*/

/** \brief 	This function tries to reduce the speckle noise of an image by applying the Gamma Maximum A Posteriori (MAP) filter.

	The user has to provide a window size and the equivalent numbers of look (enl).
	The implementation is according to the article by  
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
    
    All restrictions of the called functions \ref applyWindowFunction apply.
    
    <b> Preconditions:</b>
    \code  
	enl > 0
    \endcode
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        gammaMAPFilter(MultiArrayView<2, T1, S1> const & src,
                   			 MultiArrayView<2, T2, S2> dest,
                   			 Diff2D window_shape, int enl,
                    		 BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);

    }
    \endcode

    \deprecatedAPI{gammaMAPFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void gammaMAPFilter(SrcIterator supperleft,
                       			  SrcIterator slowerright, SrcAccessor sa,
                       			  DestIterator dupperleft, DestAccessor da,
                     			  Diff2D window_shape, int enl,
                       			  BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gammaMAPFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   			 pair<DestIterator, DestAccessor> dest,
                   			 Diff2D window_shape, int enl,
                   			 BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/specklefilters.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int w=1000, h=1000;
    MultiArray<2, float> src(w,h), dest(w,h);
    ...
	
    // apply a Gamma MAP filter with a window size of 5x5, where
    // the image was composed by 3 equivalent looks:
    gammaMAPFilter(src, dest, Diff2D(5,5), 3);
    \endcode
*/

template<typename VALUETYPE>
class GammaMAPFunctor
{
public:
	GammaMAPFunctor(Diff2D window_shape, int enl)
	: m_window_shape(window_shape),
	  m_enl(enl)
	{
		using namespace vigra;		
		vigra_precondition( enl>0, "vigra::GamaMAPFunctor(): Equivalent number of looks (enl) must be larger than zero!");
	}
	
	template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
	void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc)
	{
		using namespace vigra;
		
		SrcIterator s_ul = s - m_window_shape/2,
                    s_lr = s_ul + m_window_shape;
		
		FindAverageAndVariance<VALUETYPE> averageAndVariance;   // init functor
		inspectImage(s_ul, s_lr, s_acc, averageAndVariance);
		
		//As defined in: Shi & Fung: A comparison of Digital Speckle Filters
		/* With ENL -> C_u and ENL -> C_max from ENVI: online_help/Using_Adaptive_Filters.html */
		VALUETYPE		C_u    = 0.523/sqrt((double)m_enl),
						C_max  = sqrt(1+2.0/m_enl),
						I_mean = averageAndVariance.average(),
						C_I = sqrt(averageAndVariance.variance()) / I_mean;
		
		if(C_I <= C_u)
		{			
			d_acc.set(averageAndVariance.average(), d);
		}
		else if(C_I < C_max)
		{
			double	alpha = (1 + C_u*C_u) / (C_I*C_I - C_u*C_u),
					aL1   = alpha - m_enl - 1,
					result =		(aL1 * I_mean + sqrt(I_mean*I_mean * aL1*aL1 + 4*alpha*m_enl*I_mean)) 
								/	(2 * alpha); 
			d_acc.set(result, d);
		}
		else {
			d_acc.set(s_acc(s), d);
		}

	}
	
	Diff2D windowShape() const
	{
		return m_window_shape;
	}
	
private:
	Diff2D m_window_shape;
	int m_enl;
};

    
template <class SrcIterator, class SrcAccessor, 
		  class DestIterator, class DestAccessor>
inline void gammaMAPFilter(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                    	   DestIterator d_ul, DestAccessor d_acc, 
                    	   Diff2D window_shape, int enl,
                    	   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    GammaMAPFunctor<typename SrcIterator::value_type> func(window_shape, enl);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class SrcIterator, class SrcAccessor, 
		  class DestIterator, class DestAccessor>
inline void gammaMAPFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
                   		   pair<DestIterator, DestAccessor> d, 
                    	   Diff2D window_shape, int enl,
                    	   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    gammaMAPFilter(s.first, s.second, s.third,
                        d.first, d.second, 
                        window_shape, enl,
                        border);
}

template <class T1, class S1, 
          class T2, class S2>
inline void gammaMAPFilter(MultiArrayView<2, T1, S1> const & src,
                        	MultiArrayView<2, T2, S2> dest,
                 			Diff2D window_shape, int enl,
                 			BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
	vigra_precondition(src.shape() == dest.shape(),
						"vigra::gammaMAPFilter(): Shape mismatch between input and output.");
	gammaMAPFilter(srcImageRange(src),
						destImage(dest),  
                		window_shape, enl,
                		border);
}



/*********************************************************************
 *                                                                   *
 * The Kuan Filter (with parameter window_shape)                     *
 *                                                                   *
 *     Parameters:		window_shape   The size of the filter        *
 *                      enl            Eq. Num. Looks                *
 *********************************************************************/

/**  
	This function tries to reduce the speckle noise of an image by means of applying the 
	Kuan filter using a window of given size, and the equivalent
	numbers of look (enl). The implementation is according to the article by 
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
*/

/** \brief 	This function tries to reduce the speckle noise of an image by applying the Kuan filter.

	The user has to provide a window size and the equivalent numbers of look (enl).
	The implementation is according to the article by  
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
    
    All restrictions of the called functions \ref applyWindowFunction apply.
    
    <b> Preconditions:</b>
    \code  
	enl > 0
    \endcode
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        kuanFilter(MultiArrayView<2, T1, S1> const & src,
              	   MultiArrayView<2, T2, S2> dest,
                   Diff2D window_shape, int enl,
                   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);

    }
    \endcode

    \deprecatedAPI{kuanFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void kuanFilter(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIterator dupperleft, DestAccessor da,
                     	Diff2D window_shape, int enl,
                       	BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        kuanFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest,
                   Diff2D window_shape, int enl,
                   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/specklefilters.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int w=1000, h=1000;
    MultiArray<2, float> src(w,h), dest(w,h);
    ...
	
    // apply a Kuan filter with a window size of 5x5, where
    // the image was composed by 3 equivalent looks:
    kuanFilter(src, dest, Diff2D(5,5), 3);
    \endcode
*/
template<typename VALUETYPE>
class KuanFunctor
{
public:
	KuanFunctor(Diff2D window_shape, int enl)
	: m_window_shape(window_shape),
	  m_enl(enl)
	{
		using namespace vigra;		
		vigra_precondition( enl>0, "vigra::KuanFunctor(): Equivalent number of looks (enl) must be larger than zero!");
	}

	template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
	void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc)
	{
		using namespace vigra;
		
		SrcIterator s_ul = s - m_window_shape/2,
                    s_lr = s_ul + m_window_shape;
		
		FindAverageAndVariance<VALUETYPE> averageAndVariance;   // init functor
		inspectImage(s_ul, s_lr, s_acc, averageAndVariance);
		
		/*As defined in: Lopez & Touzi & Nezry: Adaptive speckle filters and scene heterogenity*/
		VALUETYPE	/*C_u2 = m_var_u/(m_mean_u*m_mean_u),*/
					C_u2    = (0.523*0.523)/m_enl,					 
					C_I2 = averageAndVariance.variance() / (averageAndVariance.average()*averageAndVariance.average()), 
					W    = (1 - C_u2/C_I2)/(1 + C_u2),
					I    = s_acc(s),
					R    = I * W + averageAndVariance.average() * (1 - W);
		
		d_acc.set(R, d);
	}
	
	Diff2D windowShape() const
	{
		return m_window_shape;
	}
	
private:
	Diff2D m_window_shape;
	int m_enl;
};
    
template <class SrcIterator, class SrcAccessor, 
class DestIterator, class DestAccessor>
inline void kuanFilter(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
                	   DestIterator d_ul, DestAccessor d_acc, 
                	   Diff2D window_shape, int enl,
                	   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    KuanFunctor<typename SrcIterator::value_type> func(window_shape, enl);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class SrcIterator, class SrcAccessor, 
class DestIterator, class DestAccessor>
inline void kuanFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
       			       pair<DestIterator, DestAccessor> d, 
                	   Diff2D window_shape, int enl,
                	   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    kuanFilter(s.first, s.second, s.third,
               d.first, d.second, 
               window_shape, enl,
               border);
}

template <class T1, class S1, 
          class T2, class S2>
inline void kuanFilter(MultiArrayView<2, T1, S1> const & src,
                       MultiArrayView<2, T2, S2> dest,
                 	   Diff2D window_shape, int enl,
                 	   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
	vigra_precondition(src.shape() == dest.shape(),
						"vigra::kuanFilter(): Shape mismatch between input and output.");
	kuanFilter(srcImageRange(src),
			   destImage(dest),  
               window_shape, enl,
               border);
}


/*********************************************************************
 *                                                                   *
 * The (Basic) Lee Filter                                            *
 *                                                                   *
 *     Parameters:		window_shape   The size of the filter        *
 *                      enl            Eq. Num. Looks                *
 *********************************************************************/

/**  
	This function tries to reduce the speckle noise of an image by means of applying the 
	(basic) Lee filter using a window of given size, and the equivalent
	numbers of look (enl). The implementation is according to the article by 
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
*/

/** \brief 	This function tries to reduce the speckle noise of an image by applying the basic Lee filter.

	The user has to provide a window size and the equivalent numbers of look (enl).
	The implementation is according to the article by  
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
    
    All restrictions of the called functions \ref applyWindowFunction apply.
    
    <b> Preconditions:</b>
    \code  
	enl > 0
    \endcode
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        leeFilter(MultiArrayView<2, T1, S1> const & src,
              	   MultiArrayView<2, T2, S2> dest,
                   Diff2D window_shape, int enl,
                   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);

    }
    \endcode

    \deprecatedAPI{leeFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void leeFilter(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIterator dupperleft, DestAccessor da,
                     	Diff2D window_shape, int enl,
                       	BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        leeFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest,
                   Diff2D window_shape, int enl,
                   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/specklefilters.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int w=1000, h=1000;
    MultiArray<2, float> src(w,h), dest(w,h);
    ...
	
    // apply a basic Lee filter with a window size of 5x5, where
    // the image was composed by 3 equivalent looks:
    leeFilter(src, dest, Diff2D(5,5), 3);
    \endcode
*/
template<typename VALUETYPE = float>
class LeeFunctor
{
public:
	LeeFunctor(Diff2D window_shape, int enl)
	: m_window_shape(window_shape),
	  m_enl(enl)
	{
		using namespace vigra;		
		vigra_precondition( enl>0, "vigra::LeeFunctor(): Equivalent number of looks (enl) must be larger than zero!");
	}
	
	
	template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
	void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc)
	{
		using namespace vigra;
		
		SrcIterator s_ul = s - m_window_shape/2,
                    s_lr = s_ul + m_window_shape;
		
		FindAverageAndVariance<VALUETYPE> averageAndVariance;   // init functor
		inspectImage(s_ul, s_lr, s_acc, averageAndVariance);
		
		/*As defined in: Lopez & Touzi & Nezry: Adaptive speckle filters and scene heterogenity*/
		VALUETYPE	C_u2    = (0.523*0.523)/m_enl,
					C_I2 = averageAndVariance.variance() / (averageAndVariance.average()*averageAndVariance.average()), 
					W    = (1.0 - C_u2/C_I2),
					I    = s_acc(s),
					R    = I * W + averageAndVariance.average() * (1 - W);
		
		d_acc.set(R, d);
	}
	
	Diff2D windowShape() const
	{
		return m_window_shape;
	}
	
private:
	Diff2D m_window_shape;
	int m_enl;
};

template <class SrcIterator, class SrcAccessor, 
class DestIterator, class DestAccessor>
void leeFilter(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
               DestIterator d_ul, DestAccessor d_acc, 
               Diff2D window_shape, int enl,
               BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    LeeFunctor<typename SrcIterator::value_type> func(window_shape, enl);
    applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class SrcIterator, class SrcAccessor, 
class DestIterator, class DestAccessor>
void leeFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
               pair<DestIterator, DestAccessor> d, 
               Diff2D window_shape, int enl,
               BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
    leeFilter(s.first, s.second, s.third,
              d.first, d.second, 
              window_shape, enl,
              border);
}

template <class T1, class S1, 
          class T2, class S2>
inline void leeFilter(MultiArrayView<2, T1, S1> const & src,
                       MultiArrayView<2, T2, S2> dest,
                 	   Diff2D window_shape, int enl,
                 	   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
	vigra_precondition(src.shape() == dest.shape(),
						"vigra::leeFilter(): Shape mismatch between input and output.");
	leeFilter(srcImageRange(src),
			  destImage(dest),  
              window_shape, enl,
              border);
}


/*********************************************************************
 *                                                                   *
 * The Enhanced Lee Filter                                           *
 *                                                                   *
 *     Parameters:		window_shape   The size of the filter        *
 *                      k              The damping factor (0,...,1)  *
 *                      enl            Eq. Num. Looks for comp. of   *
 *                                     the thresholds C_u and C_max  *
 *********************************************************************/
 
/**  
	This function tries to reduce the speckle noise of an image by means of applying the 
	enhanced Lee filter using a window of given size, a damping factor k, and the equivalent
	numbers of look (enl). The implementation is according to the article by 
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
*/

/** \brief 	This function tries to reduce the speckle noise of an image by applying the Enhanced Lee filter.

	The user has to provide a window size, a damping factor k, and the equivalent
	numbers of look (enl). The implementation is according to the article by  
	Lopez & Touzi & Nezry (1990): Adaptive speckle filters and scene heterogenity.
    
    All restrictions of the called functions \ref applyWindowFunction apply.
    
    <b> Preconditions:</b>
    \code  
	1. 0.0 < k <= 1.0
	2. enl > 0
    \endcode
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        enhancedLeeFilter(MultiArrayView<2, T1, S1> const & src,
                   			 MultiArrayView<2, T2, S2> dest,
                   			 Diff2D window_shape, float k, int enl,
                    		 BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);

    }
    \endcode

    \deprecatedAPI{enhancedLeeFilter}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void enhancedLeeFilter(SrcIterator supperleft,
                       			  SrcIterator slowerright, SrcAccessor sa,
                       			  DestIterator dupperleft, DestAccessor da,
                     			  Diff2D window_shape, float k, int enl,
                       			  BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        enhancedLeeFilter(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   			 pair<DestIterator, DestAccessor> dest,
                   			 Diff2D window_shape, float k, int enl,
                   			 BorderTreatmentMode border = BORDER_TREATMENT_REPEAT);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/specklefilters.hxx\><br/>
    Namespace: vigra

    \code
    unsigned int w=1000, h=1000;
    MultiArray<2, float> src(w,h), dest(w,h);
    ...
	
    // apply an enhanced Lee filter with a window size of 5x5 and a damping factor of 0.5, where
    // the image was composed by 3 equivalent looks:
    enhancedLeeFilter(src, dest, Diff2D(5,5), 0.5, 3);
    \endcode
*/

template<typename VALUETYPE>
class EnhancedLeeFunctor
{
public:
	EnhancedLeeFunctor(Diff2D window_shape, float k, int enl)
	: m_window_shape(window_shape),
      m_k(k),
      m_enl(enl)
	{
		using namespace vigra;		
		vigra_precondition( k>0 && k<=1 , "vigra::EnhancedLeeFunctor(): Damping factor k has to be: 0 < k <= 1!");
		vigra_precondition( enl>0, "vigra::EnhancedLeeFunctor(): Equivalent number of looks (enl) must be larger than zero!");
	}
	
	template <class SrcIterator,  class SrcAccessor, class DestIterator,  class DestAccessor>
	void operator()(SrcIterator s, SrcAccessor s_acc, DestIterator d, DestAccessor d_acc)
	{	
		using namespace vigra;
		
		SrcIterator s_ul = s - m_window_shape/2,
                    s_lr = s_ul + m_window_shape;
		
		
		FindAverageAndVariance<VALUETYPE> averageAndVariance;   // init functor
		inspectImage(s_ul, s_lr, s_acc, averageAndVariance);
		
		/*As defined in: Lopez & Touzi & Nezry: Adaptive speckle filters and scene heterogenity*/
		/* With ENL -> C_u and ENL -> C_max from ENVI: online_help/Using_Adaptive_Filters.html */
		VALUETYPE	C_u    = 0.523/sqrt((double)m_enl),
					C_max  = sqrt(1+2.0/m_enl),
					C_A    = sqrt(averageAndVariance.variance()) / averageAndVariance.average(), 
					W      = exp(-m_k * (C_A - C_u)/(C_max - C_A)),
					I      = s_acc(s);
		
		if( C_A <= C_u )
		{
			d_acc.set(averageAndVariance.average(), d);
		}
		else if(C_A < C_max)
		{
			d_acc.set(I * W + averageAndVariance.average() * (1 - W), d);
		}
		else {
			d_acc.set(I, d);
		}
	}
	
	Diff2D windowShape() const
	{
		return m_window_shape;
	}
	
private:
	
	Diff2D m_window_shape;
	float m_k;
	int m_enl;
};
    
template <class SrcIterator, class SrcAccessor, 
		  class DestIterator, class DestAccessor>
void enhancedLeeFilter(SrcIterator s_ul,  SrcIterator s_lr,   SrcAccessor s_acc,
					   DestIterator d_ul, DestAccessor d_acc, 
					   Diff2D window_shape, float k, int enl,
					   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
	EnhancedLeeFunctor<typename SrcIterator::value_type> func(window_shape, k, enl);
	applyWindowFunction(s_ul, s_lr, s_acc, d_ul, d_acc, func, border);
}

template <class SrcIterator, class SrcAccessor, 
		  class DestIterator, class DestAccessor>
void enhancedLeeFilter(triple<SrcIterator, SrcIterator, SrcAccessor> s,
					   pair<DestIterator, DestAccessor> d, 
					   Diff2D window_shape, float k, int enl,
					   BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
	enhancedLeeFilter(s.first, s.second, s.third,
					  d.first, d.second, 
					  window_shape, k, enl,
					  border);
}

template <class T1, class S1, 
          class T2, class S2>
inline void enhancedLeeFilter(MultiArrayView<2, T1, S1> const & src,
                        	  MultiArrayView<2, T2, S2> dest,
                 			  Diff2D window_shape, float k, int enl,
                 			  BorderTreatmentMode border = BORDER_TREATMENT_REPEAT)
{
	vigra_precondition(src.shape() == dest.shape(),
						"vigra::enhancedLeeFilter(): Shape mismatch between input and output.");
	enhancedLeeFilter(srcImageRange(src),
						destImage(dest),  
                		window_shape, k, enl,
                		border);
}

//@}

} //end of namespace vigra

#endif //VIGRA_SPECKLEFILTERS_HXX