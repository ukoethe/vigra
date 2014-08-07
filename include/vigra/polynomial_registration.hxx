/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2013 by Benjamin Seppke                 */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/************************************************************************/
#ifndef VIGRA_POLYNOMIAL_REGISTRATION_HXX
#define VIGRA_POLYNOMIAL_REGISTRATION_HXX

#include "mathutil.hxx"
#include "matrix.hxx"
#include "linear_solve.hxx"
#include "tinyvector.hxx"
#include "splineimageview.hxx"

namespace vigra
{
    
/** \addtogroup Registration Image Registration
 */
//@{

namespace detail
{    
/**
 * Iterative function for determinination of the polynom weights:
 *
 * Example: order=2, x, y
 *   -----> 
 *          [1, 
 *              x, y,
 *                    x^2, x*y, y^2]
 *
 * This function is needed, because the polynomial transformation Matrix
 * has the the same number of rows. the target position is then determined
 * by multiplying each x- and y-transformation result value with the 
 * corresponding weight for the current x- and y-coordinate, given by this
 * function.
 */
std::vector<double> polynomWeights(double x, double y, unsigned int polynom_order)
{
    unsigned int poly_count = (polynom_order+1)*(polynom_order+2)/2;
    
    std::vector<double> weights(poly_count);
    
    unsigned int weight_idx=0;
    
    for (unsigned int order=0; order<=polynom_order; order++)
    {
        for(unsigned int i=0; i<=order; i++, weight_idx++) 
        {
			weights[weight_idx] = pow(x,(double)order-i)*pow(y,(double)i);
		}
    }
    return weights;
}
	  
} //namespace detail

/********************************************************/
/*                                                      */
/*     polynomialMatrix2DFromCorrespondingPoints        */
/*                                                      */
/********************************************************/

/** \brief Create polynomial matrix of a certain degree that maps corresponding points onto each other.

    For use with \ref polynomialWarpImage() of same degree. 
    
    Since polynoms are usually non-linear functions, a special semantics is embedded to define 
    a matrix here. Each matrix consist of two rows, containing x- and y-factors of the polynom.
    
    The meaning of the matrix is explained at the example of a polynom of 2nd order:
    
    First  Row = [a_x b_x c_x d_x e_x f_x]
    Second Row = [a_y b_y c_y d_y e_y f_y]
    
    The transformed coordinate p'=[x' y'] of a position p=[x y] is then:
    
    x' = a_x + b_x*x + c_x*y + d_x*x^2 + e_x*x*y + f_x*y^2
    y' = a_y + b_y*x + c_y*y + d_y*x^2 + e_y*x*y + f_y*y^2 
    
    Note that the order of the polynom's factors is directly influenced by the
    \ref detail::polynomWeights() function and follows the intuitive scheme.
*/
template <int PolynomOrder, 
		  class SrcPointIterator, 
		  class DestPointIterator>
linalg::TemporaryMatrix<double>
polynomialMatrix2DFromCorrespondingPoints(SrcPointIterator s, SrcPointIterator s_end, 
                                          DestPointIterator d)
{
    int point_count = s_end - s;
    int poly_count = (PolynomOrder+1)*(PolynomOrder+2)/2;
    
	vigra::Matrix<double> A(point_count,poly_count), b1(point_count,1), res1(poly_count,1), b2(point_count,1), res2(poly_count,1);
    std::vector<double> weights;
	
    for (int i =0; i<point_count; ++i, ++s, ++d)
    {   
        weights = detail::polynomWeights((*d)[0], (*d)[1], PolynomOrder);
		
        for(unsigned int c=0; c<poly_count; c++)
        {
            A(i,c) = weights[c];
        }
						   
		b1(i,0)=(*s)[0];b2(i,0)=(*s)[1];		
    }
		
    if(!vigra::linearSolve(  A, b1, res1 ) || !vigra::linearSolve(  A, b2, res2 ))
		vigra_fail("polynomialMatrix2DFromCorrespondingPoints(): singular solution matrix.");
				
    vigra::Matrix<double> res(poly_count,2);
    
    for(int c=0; c<poly_count; c++)
    {
        res(c,0) = res1(c,0);
        res(c,1) = res2(c,0);
    }
    
    return res;
}
    
    
/********************************************************/
/*                                                      */
/*                polynomialWarpImage                   */
/*                                                      */
/********************************************************/

/** \brief Warp an image according to an polynomial transformation.
 	
 	To get more information about the structure of the matrix, see \ref polynomialMatrix2DFromCorrespondingPoints()
 
    <b>\#include</b> \<vigra/polynomial_registration.hxx\><br>
	Namespace: vigra

	pass 2D array views:
	\code
	namespace vigra {
        template <int ORDER, class T,
                  class T2, class S2, 
                  class C>
		void 
		polynomialWarpImage(SplineImageView<ORDER, T> const & src,
							MultiArrayView<2, T2, S2> dest,
                            MultiArrayView<2, double, C> const & polynomialMatrix);
	}
	\endcode

	\deprecatedAPI{projectiveWarpImage}
	
    pass arguments explicitly:
    \code
    namespace vigra {
        template <int ORDER, class T, 
                  class DestIterator, class DestAccessor,
                  class C>
        void polynomialWarpImage(SplineImageView<ORDER, T> const & src,
                                 DestIterator dul, DestIterator dlr, DestAccessor dest, 
                                 MultiArrayView<2, double, C> const & polynomialMatrix);
    }
    \endcode
 
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <int ORDER, class T, 
                  class DestIterator, class DestAccessor,
                  class C>
        void polynomialWarpImage(SplineImageView<ORDER, T> const & src,
                                 triple<DestIterator, DestIterator, DestAccessor> dest, 
                                 MultiArrayView<2, double, C> const & polynomialMatrix);
    }
    \endcode
	\deprecatedEnd
 */
doxygen_overloaded_function(template <...> void polynomialWarpImage)
    
template <int PolynomOrder, 
          int ORDER, class T, 
          class DestIterator, class DestAccessor,
          class C>
void polynomialWarpImage(SplineImageView<ORDER, T> const & src,
                         DestIterator dul, DestIterator dlr, DestAccessor dest, 
                         MultiArrayView<2, double, C> const & polynomialMatrix)
{
    int poly_count = (PolynomOrder+1)*(PolynomOrder+2)/2;
    
    vigra_precondition(rowCount(polynomialMatrix) == poly_count && columnCount(polynomialMatrix) == 2,
                           "polynomialWarpImage(): matrix doesn't represent a polynomial transformation of given degreee in 2D coordinates.");
	
    double w = dlr.x - dul.x;
    double h = dlr.y - dul.y;
    
    std::vector<double> weights(poly_count);
	
    for(double y = 0.0; y < h; ++y, ++dul.y)
    {
        typename DestIterator::row_iterator rd = dul.rowIterator();
        for(double x=0.0; x < w; ++x, ++rd)
        {
            weights = detail::polynomWeights(x,y, PolynomOrder);
            
            double sx=0;
            double sy=0;
            
			for(int c=0; c<poly_count; c++)
            {
				sx += weights[c]*polynomialMatrix(c,0);
				sy += weights[c]*polynomialMatrix(c,1);
			}
            
            if(src.isInside(sx, sy))
                dest.set(src(sx, sy), rd);
        }
    }
}
    
template <int PolynomOrder, 
          int ORDER, class T, 
          class DestIterator, class DestAccessor,
          class C>
inline
void polynomialWarpImage(SplineImageView<ORDER, T> const & src,
                         triple<DestIterator, DestIterator, DestAccessor> dest, 
                         MultiArrayView<2, double, C> const & polynomialMatrix)
{
    polynomialWarpImage<PolynomOrder>(src, dest.first, dest.second, dest.third, polynomialMatrix);
}
  

template <int PolynomOrder, 
          int ORDER, class T, 
          class T2, class S2,
          class C>
inline
void polynomialWarpImage(SplineImageView<ORDER, T> const & src,
                     	 MultiArrayView<2, T2, S2> dest, 
                     	 MultiArrayView<2, double, C> const & polynomialMatrix)
{
    polynomialWarpImage<PolynomOrder>(src, destImageRange(dest), polynomialMatrix);
}


//@}
    
} // namespace vigra


#endif /* VIGRA_POLYNOMIAL_REGISTRATION_HXX */
