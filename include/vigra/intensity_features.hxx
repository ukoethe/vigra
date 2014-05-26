#ifndef INTENSITY_FEATURES_HXX
#define INTENSITY_FEATURES_HXX

#include "multi_array.hxx"
#include "multi_math.hxx"


namespace vigra {

class CenterOfMassFunctor {

public:

	CenterOfMassFunctor() :
			x(0.0), y(0.0), size(0) {
	}

	void operator()(Diff2D const& diff) {
		++size;
		x += diff.x;
		y += diff.y;
	}

	float xCenter() const {
		return x / size;
	}

	float yCenter() const {
		return y / size;
	}

	Diff2D center() const {
		return Diff2D(x / size, y / size);
	}

	float x;
	float y;
	int size;
};

class Distance2PointFunctor {

public:

	Distance2PointFunctor(Diff2D point) :
			point(point) {
	}

	double operator()(Diff2D const & s) const {
		return (s - point).magnitude();
	}

	Diff2D point;
};

/********************************************************/
/*                                                      */
/*                    avg                               */
/*                                                      */
/********************************************************/

/** \brief  Find the average intensity in an image or ROI.

 <b> Usage:</b>

 \code
 vigra::MultiArray<2, double> imgArray(3, 3);

 // Do something with imgArray

 double avgIntensity = avg(imgArray);
 std::cout << avgIntensity << std::endl;


 \endcode

 */
template<class VALUETYPE>
double avg(MultiArrayView<2, VALUETYPE> inputArray) {
	FindAverage<VALUETYPE> average;
	vigra::inspectImage(inputArray, average);
	return average.average();
}

/********************************************************/
/*                                                      */
/*                    dwavg                             */
/*                                                      */
/********************************************************/

/** \brief  Find the average intensity weighted by the distance to a point in an image or ROI.

 <b> Usage:</b>

 \code
 vigra::MultiArray<2, double> imgArray(3, 3);

 // Do something with imgArray

 double avgIntensity = dwavg(imgArray, Diff2D(1, 1));
 std::cout << avgIntensity << std::endl;


 \endcode

 */
template<class VALUETYPE>
double dwavg(MultiArrayView<2, VALUETYPE> inputArray, Diff2D point) {
	Distance2PointFunctor d2p(point);
	FindAverage<VALUETYPE> average;
	MultiArray<2, double> weights(inputArray.shape());
	transformImage(
			srcIterRange(Diff2D(),
					Diff2D(inputArray.shape()[0], inputArray.shape()[1])),
			destImage(weights), d2p);
	inspectTwoImages(inputArray, weights, average);
	return average.average();
}

/********************************************************/
/*                                                      */
/*                    idwavg                            */
/*                                                      */
/********************************************************/

/** \brief  Find the average intensity weighted by the inverse of the distance to a point in an image or ROI.

 <b> Usage:</b>

 \code
 vigra::MultiArray<2, double> imgArray(3, 3);

 // Do something with imgArray

 double avgIntensity = idwavg(imgArray, Diff2D(1, 1));
 std::cout << avgIntensity << std::endl;


 \endcode

 */
template<class VALUETYPE>
double idwavg(MultiArrayView<2, VALUETYPE> inputArray, Diff2D point) {
	using namespace vigra::multi_math;
	Distance2PointFunctor d2p(point);
	FindAverage<VALUETYPE> average;
	MultiArray<2, double> weights(inputArray.shape());
	transformImage(
			srcIterRange(Diff2D(),
					Diff2D(inputArray.shape()[0], inputArray.shape()[1])),
			destImage(weights), d2p);
	weights = 1 / (weights + 1);
	inspectTwoImages(inputArray, weights, average);
	return average.average();
}

/********************************************************/
/*                                                      */
/*                    centerOfMass                      */
/*                                                      */
/********************************************************/

/** \brief  Find the center of mass in an image or ROI.

 <b> Usage:</b>

 \code
 vigra::MultiArray<2, double> imgArray(3, 3);

 // Do something with imgArray

 Diff2D center = centerOfMass(imgArray);

 std::cout << "X: " << center.x << " " << "Y: " << center.y << std::endl;


 \endcode

 */
template<class VALUETYPE>
Diff2D centerOfMass(MultiArrayView<2, VALUETYPE> inputArray) {
	CenterOfMassFunctor com;
	inspectImageIf(
			srcIterRange(Diff2D(),
					Diff2D(inputArray.shape()[0], inputArray.shape()[1])),
			srcImage(inputArray), com);
	return com.center();
}

/********************************************************/
/*                                                      */
/*                    avgDistance                       */
/*                                                      */
/********************************************************/

/** \brief  Find the average distance to a point in an image or ROI.

 If the input image will be treated like a mask ignoring pixels with value 0.

 <b> Usage:</b>

 \code
 vigra::MultiArray<2, double> imgArray(3, 3);

 // Do something with imgArray

 double avgDistance = avgDistance(imgArray, Diff2D(1, 1));
 std::cout << avgDistance << std::endl;


 \endcode

 */
template<class VALUETYPE>
double avgDistance(MultiArrayView<2, VALUETYPE> inputArray, Diff2D point) {
	Distance2PointFunctor d2p(point);
	FindAverage<VALUETYPE> average;
	MultiArray<2, double> distances(inputArray.shape());
	transformImageIf(
			srcIterRange(Diff2D(),
					Diff2D(inputArray.shape()[0], inputArray.shape()[1])),
			srcImage(inputArray), destImage(distances), d2p);
	inspectImageIf(distances, inputArray, average);
	return average.average();
}

}

#endif
