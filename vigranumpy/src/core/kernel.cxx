/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyfilters_PyArray_API
#define NO_IMPORT_ARRAY

#include <Python.h>
#include <iostream>
#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include "vigranumpykernel.hxx"

namespace python = boost::python;

namespace vigra
{

template<class T>
void pythonInitExplicitlyKernel1D(Kernel1D<T> &k, int left, int right, NumpyArray<1,T> contents)
{
	vigra_precondition(contents.size() == 1 || right-left+1 == contents.shape(0),
	          "Kernel1D::initExplicitly(): 'contents' must contain as many elements as the kernel (or just one element).");
	          
	k.initExplicitly(left,right);
	for(int i=left; i<=right; ++i)
	{
		k[i] = (contents.size() == 1)
		             ? contents(0)
		             : contents(i-left);
	}
}

#if 0 // alternative implementation
template<class KernelValueType>
void pythonInitExplicitlyKernel1D(Kernel1D<KernelValueType> & self, int left, int right,
    python::object const & args)
{
    vigra_precondition(left <= 0, "left should be <= 0");
    vigra_precondition(right >= 0, "right should be >= 0");

    if(! PySequence_Check(args.ptr()))
    {
        KernelValueType value = python::extract<KernelValueType>(args);
        self.initExplicitly(left, right) = value;
    }
    else
    {
        KernelValueType value = python::extract<KernelValueType>(args[0]);
        Kernel::InitProxy ip = self.initExplicitly(left, right) = value;
        if(python::len(args) != self.size())
        {
            std::stringstream str;
            str << "Wrong number of init values. The number must be ";
            str << self.size();
            PyErr_SetString(PyExc_ValueError, str.str().c_str());
            python::throw_error_already_set();
        }
        else
        {
            int size = self.size();
            for(int i=1; i<size; ++i)
            {
                ip,(python::extract<KernelValueType>(args[i]));
            }
        }
    }
}
#endif // #if 0

template<class T>
T pythonGetItemKernel1D(Kernel1D<T> const & self, int position)
{
    if(self.left() <= position && self.right() >= position)
    {
        return self[position];
    }
    else
    {
        std::stringstream str;
        str << "Bad position: " << position << "." << std::endl;
        str << self.left() << " <= position <= " << self.right();
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
        return 0;
    }
}

template<class T>
void pythonSetItemKernel1D(Kernel1D<T> & self, int position, T value)
{
    if(self.left() <= position && self.right() >= position)
    {
        self[position] = value;
    }
    else
    {
        std::stringstream str;
        str << "Bad position: " << position << "." << std::endl;
        str << self.left() << " <= position <= " << self.right();
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
    }
}
	
template<class T>
void pythonInitExplicitlyKernel2D(Kernel2D<T> &k,
                                  MultiArrayShape<2>::type upperleft, MultiArrayShape<2>::type lowerright,
                                  NumpyArray<2,T> contents)
{
	vigra_precondition(contents.size() == 1 || 
	                   lowerright - upperleft + MultiArrayShape<2>::type(1,1) == contents.shape(),
	          "Kernel2D::initExplicitly(): 'contents' must contain as many elements as the kernel (or just one element).");

    Point2D ul(upperleft[0], upperleft[1]), lr(lowerright[0], lowerright[1]);
    
    k.initExplicitly(ul, lr);
    for(int y = ul.y; y <= lr.y; ++y)
    {
	    for(int x = ul.x; x <= lr.x; ++x)
	    {
		    k(x,y) = (contents.size() == 1)
		                   ? contents(0)
		                   : contents(x-ul.x, y-ul.y);
	    }
    }
}

#if 0 // alternative implementation
void py2DKernel_initExplicitly(TwoDKernel & self, int upperleftX,
    int upperleftY, int lowerrightX, int lowerrightY,
    python::object const & args)
{
    vigra_precondition(upperleftX <= 0 ,
       "initExplicitly(): upperleftX must be <= 0.");
    vigra_precondition(upperleftY <= 0 ,
       "initExplicitly(): upperleftY must be <= 0.");
    vigra_precondition(lowerrightX >= 0 ,
       "initExplicitly(): lowerrightX must be >= 0.");
    vigra_precondition(lowerrightY >= 0 ,
       "initExplicitly(): lowerrightY must be >= 0.");
    Diff2D upperleft(upperleftX, upperleftY);
    Diff2D lowerright(lowerrightX, lowerrightY);
    if(! PySequence_Check(args.ptr()))
    {
        KernelValueType value = python::extract<KernelValueType>(args);
        self.initExplicitly(upperleft, lowerright) = value;
    }
    else
    {
        KernelValueType value = python::extract<KernelValueType>(args[0]);
        TwoDKernel::InitProxy ip = self.initExplicitly(upperleft, lowerright) =
            value;
        if(python::len(args) != (self.width() * self.height()))
        {
            std::stringstream str;
            str << "Wrong number of init values. The number must be ";
            str << self.width() * self.height();
            PyErr_SetString(PyExc_ValueError, str.str().c_str());
            python::throw_error_already_set();
        }
        else
        {
            int size = self.width() * self.height();
            for(int i=1; i<size; ++i)
            {
                ip,(python::extract<KernelValueType>(args[i]));
            }
        }
    }
}
#endif

template<class T>
T pythonGetItemKernel2D(Kernel2D<T> const & self, MultiArrayShape<2>::type position)
{
    if(self.upperLeft().x <= position[0] && self.lowerRight().x >= position[0] &&
       self.upperLeft().y <= position[1] && self.lowerRight().y >= position[1])
    {
        return self(position[0], position[1]);
    }
    else
    {
        std::stringstream str;
        str << "Bad position: " << position << "." << std::endl;
        str << self.upperLeft() << " <= position <= " << self.lowerRight();
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
        return 0;
    }
}

template<class T>
void pythonSetItemKernel2D(Kernel2D<T> & self, MultiArrayShape<2>::type position, T value)
{
    if(self.upperLeft().x <= position[0] && self.lowerRight().x >= position[0] &&
       self.upperLeft().y <= position[1] && self.lowerRight().y >= position[1])
    {
        self(position[0], position[1]) = value;
    }
    else
    {
        std::stringstream str;
        str << "Bad position: " << position << "." << std::endl;
        str << self.upperLeft() << " <= position <= " << self.lowerRight();
        PyErr_SetString(PyExc_ValueError, str.str().c_str());
        python::throw_error_already_set();
    }
}	

template<class T>
void defineKernels()
{
	using namespace python;
    
    docstring_options doc_options(true, true, false);

	enum_<BorderTreatmentMode>("BorderTreatmentMode")
		.value("BORDER_TREATMENT_AVOID",BORDER_TREATMENT_AVOID)
		.value("BORDER_TREATMENT_CLIP",BORDER_TREATMENT_CLIP)
		.value("BORDER_TREATMENT_REPEAT",BORDER_TREATMENT_REPEAT)
		.value("BORDER_TREATMENT_REFLECT",BORDER_TREATMENT_REFLECT)
		.value("BORDER_TREATMENT_WRAP",BORDER_TREATMENT_WRAP);

	class_<Kernel1D<T> > kernel1d("Kernel1D",
                                "Generic 1 dimensional convolution kernel.\n\n"
                                "This kernel may be used for convolution of 1 dimensional signals or "
                                "for separable convolution of multidimensional signals. "
                                "The kernel's size is given by its left() and right() "
                                "methods. The desired border treatment mode is returned by "
                                "getBorderTreatment(). "
                                "The different init functions create a kernel with the specified "
                                "properties. "
                                "For more details, see Kernel1D_ in the C++ documentation.\n\n",
                                init<>("Standard constructor::\n\n   Kernel1D()\n\nCreates an identity kernel.\n"));
	kernel1d
        .def(init< Kernel1D<T> >(args("kernel"),
            "Copy constructor::\n\n"
            "   Kernel1D(other_kernel)\n\n"))
		.def("initGaussian",
			 (void (Kernel1D<T>::*)(double,T))&Kernel1D<T>::initGaussian,
			 (arg("scale"), arg("norm")=1.0),
                "Init kernel as a sampled Gaussian function. The radius of the kernel is "
                "always 3*std_dev. 'norm' denotes the desired sum of all bins of the "
                "kernel (i.e. the kernel is corrected for the normalization error "
                "introduced by windowing the Gaussian to a finite interval). "
                "However, if norm is 0.0, the kernel is normalized to 1 by the "
                "analytic expression for the Gaussian, and no correction for the "
                "windowing error is performed.\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'gaussianKernel()'.\n\n")
		.def("initDiscreteGaussian",
			 (void (Kernel1D<T>::*)(double,T))&Kernel1D<T>::initDiscreteGaussian,
			 (arg("scale"),arg("norm")=1.0),
                "Init kernel as Lindeberg's discrete analog of the Gaussian function. "
                "The radius of the kernel is always 3*std_dev. 'norm' denotes "
                "the desired sum of all bins of the kernel.\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'discreteGaussianKernel()'.\n\n")
		.def("initGaussianDerivative",
			 (void (Kernel1D<T>::*)(double,int,T))&Kernel1D<T>::initGaussianDerivative,
			 (arg("scale"),arg("order"),arg("norm")=1.0),
                "Init kernel as a Gaussian derivative of order 'order'. The radius of "
                "the kernel is always 3*std_dev + 0.5*order. 'norm' denotes "
                "the norm of the kernel. Thus, the kernel will be corrected for "
                "the error introduced by windowing the Gaussian to a finite "
                "interval. However, if norm is 0.0, the kernel is normalized to 1 "
                "by the analytic expression for the Gaussian derivative, and no "
                "correction for the windowing error is performed.\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'gaussianDerivativeKernel()'.\n\n")
		.def("initBurtFilter",
			 &Kernel1D<T>::initBurtFilter,
			 (arg("a")=0.04785),
                "Init kernel as a 5-tap smoothing filter of the form::\n\n"
                "   [ a, 0.25, 0.5 - 2*a, 0.25, a]\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'burtFilterKernel()'.\n\n")
		.def("initBinomial",
			 (void (Kernel1D<T>::*)(int,T))&Kernel1D<T>::initBinomial,
			 (arg("radius"), arg("norm")=1.0),
                "Init kernel as a binomial filter with given radius (i.e. window size 2*radius+1). "
                "'norm' denotes the sum of all bins of the kernel.\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'binomialKernel()'.\n\n")
		.def("initAveraging",
			 (void (Kernel1D<T>::*)(int,T))&Kernel1D<T>::initAveraging,
			 (arg("radius"),arg("norm")=1.0),
                "Init kernel as an averaging filter with given radius (i.e. window size 2*radius+1). "
                "'norm' denotes the sum of all bins of the kernel.\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'averagingKernel()'.\n\n")
		.def("initSymmetricDifference",
			 (void (Kernel1D<T>::*)(T))&Kernel1D<T>::initSymmetricDifference,
			 (arg("norm")=1.0),
                "Init kernel as a symmetric difference filter of the form::\n\n"
                "   [ 0.5 * norm, 0.0 * norm, -0.5 * norm]\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'symmetricDifferenceKernel()'.\n\n")
		.def("initSecondDifference3",
			 &Kernel1D<T>::initSecondDifference3,
                "Init kernel as a 3-tap second difference filter of the form::\n\n"
                "   [ 1, -2, 1]\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'secondDifference3Kernel()'.\n\n")
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
		.def("initOptimalFirstDerivative5",
			 &Kernel1D<T>::initOptimalFirstDerivative5)
		.def("initOptimalSecondDerivative5",
			 &Kernel1D<T>::initOptimalSecondDerivative5)
		.def("initExplicitly",
			 registerConverters(&pythonInitExplicitlyKernel1D<T>),
			 (arg("left"), arg("right"), arg("contents")),
                "Init the kernel with explicit values from 'contents', which must be a "
                "1D numpy.ndarray. 'left' and 'right' are the boundaries of the kernel "
                "(inclusive). If 'contents' contains the wrong number of values, a "
                "run-time error results. It is, however, possible to give just one "
                "initializer. This creates an averaging filter with the given constant. "
                "The norm is set to the sum of the initializer values. \n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'explicitlyKernel()'.\n\n")
		.def("__getitem__",
			 &pythonGetItemKernel1D<T>)
		.def("__setitem__",
			 &pythonSetItemKernel1D<T>)
		.def("left",
			 &Kernel1D<T>::left,
                "Left border of kernel (inclusive).\n")
		.def("right",
			 &Kernel1D<T>::right,
                "Right border of kernel (inclusive).\n")
		.def("size",
			 &Kernel1D<T>::size,
			    "Number of kernel elements (right() - left() + 1).\n")
		.def("borderTreatment",
			 &Kernel1D<T>::borderTreatment,
                "Return current border treatment mode.\n")
		.def("setBorderTreatment",
			 &Kernel1D<T>::setBorderTreatment,
			 args("borderTreatment"),
                "Set border treatment mode.\n")
		.def("norm",
			 &Kernel1D<T>::norm,
                "Return the norm of kernel.\n")
		.def("normalize",
			 (void (Kernel1D<T>::*)(T,unsigned int,double))&Kernel1D<T>::normalize,
			 (arg("norm")=1.0,arg("derivativeOrder")=0,arg("offset")= 0.0),
                "Set a new norm and normalize kernel, use the normalization "
                "formula for the given derivativeOrder.\n")
		;

	class_<Kernel2D<T> > kernel2d("Kernel2D",
	        "Generic 2 dimensional convolution kernel.\n\n"
            "This kernel may be used for convolution of 2 dimensional signals. "
            "The desired border treatment mode is returned by borderTreatment()."
            "(Note that the 2D convolution functions don't currently support all "
            "modes.) "
            "The different init functions create a kernel with the specified "
            "properties. "
            "For more details, see Kernel2D_ in the C++ documentation.\n\n",
            init<>("Standard constructor::\n\n   Kernel2D()\n\nCreates an identity kernel.\n"));
	kernel2d
        .def(init< Kernel2D<T> >(args("kernel"),
            "Copy constructor::\n\n"
            "   Kernel2D(other_kernel)\n\n"))
		.def("initExplicitly",
			 registerConverters(&pythonInitExplicitlyKernel2D<T>),
			 (arg("upperLeft"), arg("lowerRight"), arg("contents")),
                "Init the kernel with explicit values from 'contents', which must be a "
                "2D numpy.ndarray. 'upperLeft' and 'lowerRight' are the boundaries of the "
                "kernel (inclusive), and  must be 2D tuples. "
                "If 'contents' contains the wrong number of values, a run-time error "
                "results. It is, however, possible to give just one initializer. "
                "This creates an averaging filter with the given constant. "
                "The norm is set to the sum of the initializer values. \n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'explicitlyKernel2D()'.\n\n")
		.def("initSeparable",
			 (void (Kernel2D<T>::*)(Kernel1D<T> const &,Kernel1D<T> const &))&Kernel2D<T>::initSeparable,
			 (arg("kernelX"), arg("kernelY")),
		        "Init the 2D kernel as the cartesian product of two 1D kernels of "
                "type Kernel1D. The norm becomes the product of the two original "
                "norms.\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'separableKernel2D()'.\n\n")
		.def("initGaussian",
			 (void (Kernel2D<T>::*)(double, T))&Kernel2D<T>::initGaussian,
			 (arg("scale"), arg("norm")=1.0),
                "Init kernel as a sampled 2D Gaussian function. The radius of the kernel is "
                "always 3*std_dev. 'norm' denotes the desired sum of all bins of the "
                "kernel (i.e. the kernel is corrected for the normalization error "
                "introduced by windowing the Gaussian to a finite interval). "
                "However, if norm is 0.0, the kernel is normalized to 1 by the "
                "analytic expression for the Gaussian, and no correction for the "
                "windowing error is performed.\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'gaussianKernel2D()'.\n\n")
		.def("initDisk",
			 &Kernel2D<T>::initDisk,
			 args("radius"),
                "Init the 2D kernel as a circular averaging filter. The norm will "
                "be calculated as 1 / (number of non-zero kernel values).\n\n"
                "Precondition::\n\n"
                "   radius > 0\n\n"
                "Kernel construction and initialization can be performed in one step "
                "by calling the factory function 'diskKernel2D()'.\n\n")
		.def("__setitem__",
			 &pythonSetItemKernel2D<T>)
		.def("__getitem__",
			 &pythonGetItemKernel2D<T>)
		.def("width",
			 &Kernel2D<T>::width,
			    "Horizontal kernel size (lowerRight()[0] - upperLeft()[0] + 1).\n")
		.def("height",
			 &Kernel2D<T>::height,
			    "Vertical kernel size (lowerRight()[1] - upperLeft()[1] + 1).\n")
		.def("upperLeft",
			 &Kernel2D<T>::upperLeft,
                "Upper left border of kernel (inclusive).\n")
		.def("lowerRight",
			 &Kernel2D<T>::lowerRight,
                "Lower right border of kernel (inclusive).\n")
		.def("norm",
			 &Kernel2D<T>::norm,
			    "Return the norm of the kernel.\n")
		.def("normalize",
		    (void (Kernel2D<T>::*)(T))&Kernel2D<T>::normalize,
			 (arg("norm")=1.0),
			    "Set the kernel's norm and renormalize the values.\n")
		.def("borderTreatment",
			 &Kernel2D<T>::borderTreatment,
                "Return current border treatment mode.\n")
		.def("setBorderTreatment",
			 &Kernel2D<T>::setBorderTreatment,
			 args("borderTreatment"),
                "Set border treatment mode.\n")
	;
}

void defineKernels()
{
	defineKernels<KernelValueType>();
}

} // namespace vigra

