/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyconvolution_PyArray_API
#define NO_IMPORT_ARRAY

#include <iostream>
#include <Python.h>
#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/separableconvolution.hxx>
#include <vigra/stdconvolution.hxx>
#include <set>
#include "tuples.hpp"
#include "vigranumpykernel.hxx"

#include <cmath>

namespace python = boost::python;
namespace vigra
{
	template<class T>
	void pythonInitSetExplicitlyKernel1D(Kernel1D<T> *k,int left,int right,NumpyArray<1,T> contents)
	{
		vigra_precondition(right-left+1==contents.shape(0),"The distance between left and right must match the length of contents");
		k->initExplicitly(left,right);
		for(int ii=left;ii<=right;++ii)
		{
			(*k)[ii]=contents[ii-left];
		}
	}
	template<class T>
		T pythonGetItemKernel1D(Kernel1D<T> *k,int i)
		{
			return (*k)[i];
		}
	template<class T>
		void pythonSetItemKernel1D(Kernel1D<T> *k,int i,T v)
		{
			(*k)[i]=v;
		}
	template<class T>
		void pythonInitSetExplicitlyKernel2D(Kernel2D<T> *k,Point2D left,Point2D right,NumpyArray<2,T> contents)
		{
			vigra_precondition(right.x-left.x+1==contents.shape(0),"The distance between left and right must match the length of contents");
			vigra_precondition(right.y-left.y+1==contents.shape(1),"The distance between left and right must match the length of contents");
			k->initExplicitly(left,right);
			int xx,yy;
			for(xx=left.x;xx<=right.x;++xx)
			{
				for(yy=left.y;yy<=right.y;++yy)
				{
					(*k)(xx,yy)=contents(xx-left.x,yy-left.y);
				}
			}
		}
	template<class T>
		void pythonInitWithFactoryKernel2D(Kernel2D<T> *k, boost::tuples::tuple< Point2D, Point2D, NumpyArray<2,T> > args)
		{
			//pythonInitSetExplicitlyKernel2D(k, args.get<0>(), args.get<1>(), args.get<2>() );
			pythonInitSetExplicitlyKernel2D(k, boost::tuples::get<0>(args), boost::tuples::get<1>(args), boost::tuples::get<2>(args) );
		}
	template<class T>
		static boost::tuples::tuple< Point2D, Point2D, NumpyArray<2,T> > kernelGaussian2D(int radius)
		{
			Kernel2D<T> k;
			k.initDisk(radius);
			Point2D ul = k.upperLeft();
			Point2D lr = k.lowerRight();
			NumpyArray<2,T> arr( MultiArrayShape<2>::type(k.width(), k.height()) );

			// hack: copy kernel array
			for (int x=ul.x;x<=lr.x;x++)
				for (int y=ul.y;y<=ul.y;y++)
					arr(x,y)=k(x,y);
			//return triple<Point2D, Point2D, NumpyArray<2,T> > (ul, lr, arr);
			return boost::tuples::tuple< Point2D, Point2D, NumpyArray<2,T> > (ul, lr, arr);
		}
	template<class T>
		static boost::tuples::tuple< Point2D, Point2D, NumpyArray<2,T> > kernelDisk2D(int radius)
		{
			Kernel2D<T> k;
			k.initDisk(radius);
			//return k
			Point2D ul = k.upperLeft();
			Point2D lr = k.lowerRight();
			NumpyArray<2,T> arr( MultiArrayShape<2>::type(k.width(), k.height()) );

			// hack: copy kernel array
			for (int x=ul.x;x<=lr.x;x++)
				for (int y=ul.y;y<=ul.y;y++)
					arr(x,y)=k(x,y);
			//return triple<Point2D, Point2D, NumpyArray<2,T> > (ul, lr, arr);
			return boost::tuples::tuple< Point2D, Point2D, NumpyArray<2,T> > (ul, lr, arr);
		}
	template<class T>
		T pythonGetItemKernel2D(Kernel2D<T> *k,Point2D p)
		{
			return(*k)(p.x,p.y);
		}
	template<class T>
		void pythonSetItemKernel2D(Kernel2D<T> *k,Point2D p,T v)
		{
			(*k)(p.x,p.y)=v;
		}

	template<class T>
	void defineKernels()
	{
		using namespace python;

		typedef boost::tuples::tuple< Point2D, Point2D, NumpyArray<2,T> > ktriplet;
		boost::python::register_tuple< ktriplet >();
		//boost::python::register_tuple< typename kTriplet<T>::Type >;

		enum_<BorderTreatmentMode>("BorderTreatmentMode")
			.value("BORDER_TREATMENT_AVOID",BORDER_TREATMENT_AVOID)
			.value("BORDER_TREATMENT_CLIP",BORDER_TREATMENT_CLIP)
			.value("BORDER_TREATMENT_REPEAT",BORDER_TREATMENT_REPEAT)
			.value("BORDER_TREATMENT_REFLECT",BORDER_TREATMENT_REFLECT)
			.value("BORDER_TREATMENT_WRAP",BORDER_TREATMENT_WRAP);

		class_<Kernel1D<T> > kernel1d("Kernel1D",init<>());
		kernel1d
			/*.def("setConst",
				 (void (Kernel1D::*)(float))&Kernel1D<T>::operator=)*/
			.def("initGaussian",
				 (void (Kernel1D<T>::*)(double,T))&Kernel1D<T>::initGaussian,
				 (arg("std_dev"), arg("norm")=1.0))
			.def("initDiscreteGaussian",
				 (void (Kernel1D<T>::*)(double,T))&Kernel1D<T>::initDiscreteGaussian,
				 (arg("std_dev"),arg("norm")=1.0))
			.def("initGaussianDerivative",
				 (void (Kernel1D<T>::*)(double,int,T))&Kernel1D<T>::initGaussianDerivative,
				 (arg("std_dev"),arg("order"),arg("norm")=1.0))
			.def("initOptimalSmoothing3",
				 &Kernel1D<T>::initOptimalSmoothing3)
			.def("initOptimalFirstDerivativeSmoothing3",
				 &Kernel1D<T>::initOptimalFirstDerivativeSmoothing3)
			.def("initOptimalSecondDerivativeSmoothing3",
				 &Kernel1D<T>::initOptimalSecondDerivativeSmoothing3)
			.def("initOptimalSmoothing5",
				 &Kernel1D<T>::initOptimalSmoothing5)
			.def("initOptimalFirstDerivativeSmoothing5",
				 &Kernel1D<T>::initOptimalFirstDerivativeSmoothing5)
			.def("initOptimalSecondDerivativeSmoothing5",
				 &Kernel1D<T>::initOptimalSecondDerivativeSmoothing5)
			.def("initBurtFilter",
				 &Kernel1D<T>::initBurtFilter,
				 (arg("a")=0.04785))
			.def("initBinomial",
				 (void (Kernel1D<T>::*)(int,T))&Kernel1D<T>::initBinomial,
				 (arg("radius"), arg("norm")=1.0))
			.def("initAveraging",
				 (void (Kernel1D<T>::*)(int,T))&Kernel1D<T>::initAveraging,
				 (arg("radius"),arg("norm")=1.0))
			.def("initSymmetricGradient",
				 (void (Kernel1D<T>::*)(T))&Kernel1D<T>::initSymmetricGradient,
				 (arg("norm")=1.0))
			.def("initSymmetricDifference",
				 (void (Kernel1D<T>::*)(T))&Kernel1D<T>::initSymmetricDifference,
				 (arg("norm")=1.0))
			.def("initSecondDifference3",
				 &Kernel1D<T>::initSecondDifference3)
			.def("initOptimalFirstDerivative5",
				 &Kernel1D<T>::initOptimalFirstDerivative5)
			.def("initOptimalSecondDerivative5",
				 &Kernel1D<T>::initOptimalSecondDerivative5)
			.def("initSetExplicitly",
				 registerConverters(&pythonInitSetExplicitlyKernel1D<T>))
			.def("__getitem__",
				 &pythonGetItemKernel1D<T>)
			.def("__setitem__",
				 &pythonSetItemKernel1D<T>)
			.def("left",
				 &Kernel1D<T>::left)
			.def("right",
				 &Kernel1D<T>::right)
			.def("size",
				 &Kernel1D<T>::size)
			.def("borderTreatment",
				 &Kernel1D<T>::borderTreatment)
			.def("setBorderTreatment",
				 &Kernel1D<T>::setBorderTreatment)
			.def("norm",
				 &Kernel1D<T>::norm)
			.def("normalize",
				 (void (Kernel1D<T>::*)(T,unsigned int,double))&Kernel1D<T>::normalize,
				 (arg("norm")=1.0,arg("derivativeOrder")=0,arg("offset")= 0.0))
			;

		class_<Kernel2D<T> > kernel2d("Kernel2D",init<>());
		kernel2d
			.def("initSetExplicitly",
				 registerConverters(&pythonInitSetExplicitlyKernel2D<T>))
			.def("initWithFactoryKernel",
				 registerConverters(&pythonInitWithFactoryKernel2D<T>))
			.def("kernelGaussian", kernelGaussian2D<T>).staticmethod("kernelGaussian")
			.def("kernelDisk", kernelDisk2D<T>).staticmethod("kernelDisk")
			.def("initSeperable",
				 (void (Kernel2D<T>::*)(Kernel1D<T>&,Kernel1D<T>&))&Kernel2D<T>::initSeparable)
			.def("__setitem__",
				 &pythonSetItemKernel2D<T>)
			.def("__getitem__",
				 &pythonGetItemKernel2D<T>)
			.def("initDisk",
				 &Kernel2D<T>::initDisk)
			.def("width",
				 &Kernel2D<T>::width)
			.def("height",
				 &Kernel2D<T>::height)
			.def("norm",
				 &Kernel2D<T>::norm)
			.def("upperLeft",
				 &Kernel2D<T>::upperLeft)
			.def("lowerRight",
				 &Kernel2D<T>::lowerRight);
	}
	void defineKernels()
	{
		defineKernels<KernelValueType>();
	}

	//void registerNumpyArrayConverters();
} // namespace vigra

