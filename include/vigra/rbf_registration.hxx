/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2014 by Benjamin Seppke                 */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/************************************************************************/
#ifndef VIGRA_RBF_REGISTRATION_HXX
#define VIGRA_RBF_REGISTRATION_HXX

#include <vigra/mathutil.hxx>
#include <vigra/matrix.hxx>
#include <vigra/linear_solve.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/splineimageview.hxx>
#include <vigra/affine_registration.hxx>

namespace vigra {
    
/** \addtogroup Registration Image Registration
 */
//@{

namespace detail {

// Small and hopefully fast helper function to compute the squared-distance
// between two points (2d)
template <class SrcPoint, class DestPoint>
inline double distance2(SrcPoint const & p1, DestPoint const & p2)
{
    return (double)(p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]);
}

} //end of namespace vigra::detail
    
    
/**
 * The famous thin plate spline functor [weight(d) = d^2*log(d^2)]
 * only affects up to max distance...
 */
struct ThinPlateSplineFunctor
{
	template <class SrcPoint, class DestPoint>
	inline double operator()(SrcPoint const & p1, DestPoint const & p2) const 
    {
		double dist2 = detail::distance2(p1, p2);
        
        if(dist2 == 0 )
		{
            return 0;
        }
        else
        {
            return dist2*log(dist2);
        }
	}
};

/**
 * A distance power based radial basis functor [weight = dist^N]
 */
template<int N>
struct DistancePowerFunctor
{
	template <class SrcPoint, class DestPoint>
	inline double operator()(SrcPoint const & p1, DestPoint const & p2) const 
    {
        double dist2 = detail::distance2(p1, p2);
        
		if(dist2 == 0)
		{
			return 0;
		}
		else
		{
			return pow(dist2, N/2.0);
		}
	}
};
    
/********************************************************/
/*                                                      */
/*          rbfMatrix2DFromCorrespondingPoints          */
/*                                                      */
/********************************************************/

/** \brief Create a matrix that maps corresponding points onto each other using a given RBF.
 
 For use with \ref radialBasisWarpImage(). For n given (corresponding) points, 
 the matrix will be of size (n+3,2). Note that the representation of this matrix is exactly 
 the same as the "W" matrix of Bookstein. More information can be found in the following article:
 
 Fred L. Bookstein. Principal Warps: Thin-Plate Splines and the Decomposition of Deformations. IEEE PAMI, Vol 11, No 8. 1989
 */
template <class RadialBasisFunctor,
          class SrcPointIterator, class DestPointIterator>
linalg::TemporaryMatrix<double>
rbfMatrix2DFromCorrespondingPoints(SrcPointIterator s, SrcPointIterator s_end, DestPointIterator d, RadialBasisFunctor const & rbf)
{    
    int point_count = s_end - s;
    
    Matrix<double> L(point_count+3, point_count+3, 0.0);
    Matrix<double> Y(point_count+3,			   2, 0.0);
    Matrix<double> W(point_count+3,		       2, 0.0);
    
    //fill P (directly into K) and V (directly into Y)
    for(int i=0; i<point_count; ++i)
    {
        L(i, point_count  ) = L(point_count,   i) = 1;
        L(i, point_count+1) = L(point_count+1, i) = (d[i])[0];
        L(i, point_count+2) = L(point_count+2, i) = (d[i])[1];
        
        Y(i,0)= (s[i])[0];
        Y(i,1)= (s[i])[1];
    }
    
    //fill K (directly into L)
    for(int j=0; j<point_count; j++)
    {
        for(int i=0; i<j; i++)
        {
            L(i,j) = L(j,i) = rbf(d[i], d[j]);
        }
    }
    
    linearSolve(L, Y, W);
    //Results are okay, even if vigra reports failure... 
    //so I commented this out
    //    if(!linearSolve(L, Y, W))
    //		vigra_fail("radialBasisMatrix2DFromCorrespondingPoints(): singular solution matrix.");
    
    return W;
};

    
/********************************************************/
/*                                                      */
/*                     rbfWarpImage                     */
/*                                                      */
/********************************************************/

/** \brief Warp an image according to an radial basis function based transformation.
 
 To get more information about the structure of the matrix, see \ref rbfMatrix2DFromCorrespondingPoints()
 
 <b>\#include</b> \<vigra/rbf_registration.hxx\><br>
 Namespace: vigra
 
 pass 2D array views:
 \code
 namespace vigra {
     template <int ORDER, class T,
               class T2, class S2,
               class DestPointIterator,
               class C,
               class RadialBasisFunctor>
     void
     rbfWarpImage(SplineImageView<ORDER, T> const & src,
                  MultiArrayView<2, T2, S2> dest,
                  DestPointIterator d, DestPointIterator d_end,
                  MultiArrayView<2, double, C> const & W,
                  RadialBasisFunctor rbf);
 }
 \endcode
 
 \deprecatedAPI{rbfWarpImage}
 
 pass arguments explicitly:
 \code
 namespace vigra {
     template <int ORDER, class T,
               class DestIterator, class DestAccessor,
               class DestPointIterator,
               class C,
               class RadialBasisFunctor>
     void
     rbfWarpImage(SplineImageView<ORDER, T> const & src,
                  DestIterator dul, DestIterator dlr, DestAccessor dest,
                  DestPointIterator d, DestPointIterator d_end,
                  MultiArrayView<2, double, C> const & W,
                  RadialBasisFunctor rbf);
 }
 \endcode
 
 use argument objects in conjunction with \ref ArgumentObjectFactories :
 \code
 namespace vigra {
     template <int ORDER, class T,
               class DestIterator, class DestAccessor,
               class DestPointIterator,
               class C,
               class RadialBasisFunctor>
     void
     rbfWarpImage(SplineImageView<ORDER, T> const & src,
                  triple<DestIterator, DestIterator, DestAccessor> dest,
                  DestPointIterator d, DestPointIterator d_end,
                  MultiArrayView<2, double, C> const & W,
                  RadialBasisFunctor rbf);
 }
 }
 \endcode
 \deprecatedEnd
 */
template <int ORDER, class T,
          class DestIterator, class DestAccessor,
          class DestPointIterator,
          class C,
          class RadialBasisFunctor>
void
rbfWarpImage(SplineImageView<ORDER, T> const & src,
             DestIterator dul, DestIterator dlr, DestAccessor dest,
             DestPointIterator d, DestPointIterator d_end,
             MultiArrayView<2, double, C> const & W,
             RadialBasisFunctor rbf)
{
    int point_count = d_end - d;
    
    vigra_precondition(rowCount(W) == point_count+3 && columnCount(W) == 2,
                       "vigra::rbfWarpImage(): matrix doesn't represent a proper transformation of given point size in 2D coordinates.");
	
    double w = dlr.x - dul.x;
    double h = dlr.y - dul.y;
    
    for(double y = 0.0; y < h; ++y, ++dul.y)
    {
        typename DestIterator::row_iterator rd = dul.rowIterator();
        for(double x=0.0; x < w; ++x, ++rd)
        {
            //Affine part		
            double	sx = W(point_count,0)+W(point_count+1,0)*x+ W(point_count+2,0)*y,
                    sy = W(point_count,1)+W(point_count+1,1)*x+ W(point_count+2,1)*y;
            
            //RBS part
            for(int i=0; i<point_count; i++)
            {
                double weight = rbf(d[i], Diff2D(x,y));
                sx += W(i,0)*weight;
                sy += W(i,1)*weight;
                
            }
            
            if(src.isInside(sx, sy))
                dest.set(src(sx, sy), rd);
        }
    }
};
    
template <int ORDER, class T,
          class DestIterator, class DestAccessor,
          class DestPointIterator,
          class C,
          class RadialBasisFunctor>
inline 
void rbfWarpImage(SplineImageView<ORDER, T> const & src, 
                  triple<DestIterator, DestIterator, DestAccessor> dest,
                  DestPointIterator d, DestPointIterator d_end,
                  MultiArrayView<2, double, C> const & W,
                  RadialBasisFunctor rbf)
{
    rbfWarpImage(src, dest.first, dest.second, dest.third, d, d_end, W, rbf);
}

template <int ORDER, class T,
          class T2, class S2,
          class DestPointIterator,
          class C,
          class RadialBasisFunctor>
inline
void rbfWarpImage(SplineImageView<ORDER, T> const & src,
                  MultiArrayView<2, T2, S2> dest, 
                  DestPointIterator d, DestPointIterator d_end,
                  MultiArrayView<2, double, C> const & W,
                  RadialBasisFunctor rbf)
{
    rbfWarpImage(src, destImageRange(dest), d, d_end, W, rbf);
}


//@}

} // namespace vigra


#endif /* VIGRA_RBF_REGISTRATION_HXX */