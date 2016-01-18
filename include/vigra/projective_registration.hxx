/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2014 by Benjamin Seppke                 */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/************************************************************************/
#ifndef VIGRA_PROJECTIVE_REGISTRATION_HXX
#define VIGRA_PROJECTIVE_REGISTRATION_HXX

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

/********************************************************/
/*                                                      */
/*     projectiveMatrix2DFromCorrespondingPoints        */
/*                                                      */
/********************************************************/

/** \brief Create homogeneous matrix that maps corresponding points onto each other.
 
    For use with \ref projectiveWarpImage(). Since four corresponding points are needed to be given,
    the matrix will compute a full projective transform.
 */
template <class SrcPointIterator, class DestPointIterator>
linalg::TemporaryMatrix<double> 
projectiveMatrix2DFromCorrespondingPoints(SrcPointIterator s, SrcPointIterator send, DestPointIterator d)
{
	//Calculate the matrix using least squares of all points of points like the result is:
	// ( x2 )   ( s_x		r1		t_x )   ( x1 )
	// ( y2 ) = (  r2		s_y		t_y ) * ( y1 ) 
	// (  1 )   (  p1		p2		1   )   (  1 )
	int size = send - s;
	
    vigra_assert(size >= 4, 
                 "projectiveMatrix2DFromCorrespondingPoints(): need at least four corresponding points.");
    
	vigra::Matrix<double> A(2*size,8, 0.0), b(2*size,1),  res(8,1);
	for (int i =0; i<size; ++i, ++s, ++d)
    {
		//m_00				m_01				m_02				m_10					m_11					m_12				m_20								m_21						
		//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		A(i,0)=(*d)[0];		A(i,1)=(*d)[1];		A(i,2)=1;			A(i,3)=0;				A(i,4)=0;				A(i,5)=0;			A(i,6)=-1*((*d)[0])*((*s)[0]);			A(i,7)=-1*((*d)[1])*((*s)[0]);	
		b(i,0)=(*s)[0];
		
		A(size+i,0)=0;		A(size+i,1)=0;		A(size+i,2)=0;		A(size+i,3)=(*d)[0];	A(size+i,4)=(*d)[1];	A(size+i,5)=1;		A(size+i,6)=-1*((*d)[0])*((*s)[1]);		A(size+i,7)=-1*((*d)[1])*((*s)[1]);	
		b(size+i,0)=(*s)[1];
		
	}
	
	vigra_assert(linearSolve(A, b, res), 
				"projectiveMatrix2DFromCorrespondingPoints(): singular solution matrix.");
    
    linalg::TemporaryMatrix<double> projectiveMat(3,3);
	projectiveMat(0,0) = res(0,0);		projectiveMat(0,1) = res(1,0);		projectiveMat(0,2) = res(2,0);
	projectiveMat(1,0) = res(3,0);		projectiveMat(1,1) = res(4,0);		projectiveMat(1,2) = res(5,0);
	projectiveMat(2,0) = res(6,0);		projectiveMat(2,1) = res(7,0);		projectiveMat(2,2) = 1;
	
	return projectiveMat;
}
    
/********************************************************/
/*                                                      */
/*                projectiveWarpImage                   */
/*                                                      */
/********************************************************/

/** \brief Warp an image according to an projective transformation.

    Sorry, no \ref detailedDocumentation() available yet.

	<b> Declarations:</b>

	<b>\#include</b> \<vigra/projective_registration.hxx\><br>
	Namespace: vigra

	pass 2D array views:
	\code
	namespace vigra {
        template <int ORDER, class T,
                  class T2, class S2, 
                  class C>
		void 
		projectiveWarpImage(SplineImageView<ORDER, T> const & src,
							MultiArrayView<2, T2, S2> dest,
                            MultiArrayView<2, double, C> const & projectiveMatrix);
	}
	\endcode

	\deprecatedAPI{projectiveWarpImage}
	pass \ref ImageIterators and \ref DataAccessors :
	
    pass arguments explicitly:
    \code
    namespace vigra {
        template <int ORDER, class T, 
                class DestIterator, class DestAccessor,
                class C>
        void projectiveWarpImage(SplineImageView<ORDER, T> const & src,
                            	 DestIterator dul, DestIterator dlr, DestAccessor dest, 
                            	 MultiArrayView<2, double, C> const & projectiveMatrix);
    }
    \endcode
    
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <int ORDER, class T, 
                class DestIterator, class DestAccessor,
                class C>
        void projectiveWarpImage(SplineImageView<ORDER, T> const & src,
                            	 triple<DestIterator, DestIterator, DestAccessor> dest, 
                            	 MultiArrayView<2, double, C> const & projectiveMatrix);
    }
    \endcode
	\deprecatedEnd
    
    The algorithm applies the given \a projectiveMatrix to the <i>destination coordinates</i> and copies
    the image value from the resulting source coordinates, using the given SplineImageView \a src for interpolation. 
    If the resulting coordinate is outside the source image, nothing will be written at that destination point.
    
    \code
        for all dest pixels:
            currentSrcCoordinate = projectiveMatrix * currentDestCoordinate;
            if src.isInside(currentSrcCoordinate):
                dest[currentDestCoordinate] = src[currentSrcCoordinate]; // copy an interpolated value
    \endcode
    
    The matrix represent a 2-dimensional projective transform by means of homogeneous coordinates,
    i.e. it must be a 3x3 matrix whose last row is (p1,p2,1).
    
    <b> Required Interface:</b>
    
    \code
    DestImageIterator dest_upperleft;
    
    double x = ..., y = ...;
    
    if (spline.isInside(x,y))
        dest_accessor.set(spline(x, y), dest_upperleft);

    \endcode
    
    <b>See also:</b> Functions to specify projective transformation: \ref translationMatrix2D(), \ref scalingMatrix2D(), 
                    \ref shearMatrix2D(), \ref rotationMatrix2DRadians(), \ref rotationMatrix2DDegrees() and \ref projectiveMatrix2DFromCorrespondingPoints()
*/
doxygen_overloaded_function(template <...> void projectiveWarpImage)

template <int ORDER, class T, 
          class DestIterator, class DestAccessor,
          class C>
void projectiveWarpImage(SplineImageView<ORDER, T> const & src,
                     DestIterator dul, DestIterator dlr, DestAccessor dest, 
                     MultiArrayView<2, double, C> const & projectiveMatrix)
{
    vigra_precondition(rowCount(projectiveMatrix) == 3 && columnCount(projectiveMatrix) == 3 && projectiveMatrix(2,2) == 1.0,
        "projectiveWarpImage(): matrix doesn't represent an projective transformation with homogeneous 2D coordinates.");
         
    
    double w = dlr.x - dul.x;
    double h = dlr.y - dul.y;
    
    for(double y = 0.0; y < h; ++y, ++dul.y)
    {
        typename DestIterator::row_iterator rd = dul.rowIterator();
        for(double x=0.0; x < w; ++x, ++rd)
        {
        	double fac = 1.0/(x*projectiveMatrix(2,0) + y*projectiveMatrix(2,1) + 1);
            double sx = (x*projectiveMatrix(0,0) + y*projectiveMatrix(0,1) + projectiveMatrix(0,2)) * fac;
            double sy = (x*projectiveMatrix(1,0) + y*projectiveMatrix(1,1) + projectiveMatrix(1,2)) * fac;
            if(src.isInside(sx, sy))
                dest.set(src(sx, sy), rd);
        }
    }
}

template <int ORDER, class T, 
          class DestIterator, class DestAccessor,
          class C>
inline
void projectiveWarpImage(SplineImageView<ORDER, T> const & src,
                     triple<DestIterator, DestIterator, DestAccessor> dest, 
                     MultiArrayView<2, double, C> const & projectiveMatrix)
{
    projectiveWarpImage(src, dest.first, dest.second, dest.third, projectiveMatrix);
}


template <int ORDER, class T, 
          class T2, class S2,
          class C>
inline
void projectiveWarpImage(SplineImageView<ORDER, T> const & src,
                     	 MultiArrayView<2, T2, S2> dest, 
                     	 MultiArrayView<2, double, C> const & projectiveMatrix)
{
    projectiveWarpImage(src, destImageRange(dest), projectiveMatrix);
}


//@}
    
} // namespace vigra


#endif /* VIGRA_PROJECTIVE_REGISTRATION_HXX */