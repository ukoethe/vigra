/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/
 
 
#ifndef VIGRA_STDCONVOLUTION_HXX
#define VIGRA_STDCONVOLUTION_HXX

#include <math.h>
#include "vigra/stdimage.hxx"
#include "vigra/bordertreatment.hxx"
#include "vigra/separableconvolution.hxx"

/** @name Two-dimensional convolution functions
    
    These generic convolution functions implement
    the standard 2D convolution operation for images that fit 
    into the required interface. Arbitrary ROI's are supported
    by the mask version of the algorithm.
    The functions need a suitable 2D kernel to operate.
    
    @memo 2D non-separable convolution, with and without ROI mask
*/
//@{

/** Performs a 2 dimensional convolution of the source image using the given
    kernel. The KernelIterator must point to the center of the kernel, and
    the kernel's size is given by its upper left (x and y of distance <= 0) and 
    lower right (distance >= 0) corners. The image must always be larger than the 
    kernel. At those positions where the kernel does not completely fit 
    into the image, the specified \Ref{BorderTreatmentMode} is 
    applied. Only BORDER_TREATMENT_CLIP and BORDER_TREATMENT_AVOID are currently 
    supported.
    
    The images's pixel type (SrcAccessor::value_type) must be a
    linear space over the kernel's value_type (KernelAccessor::value_type),
    i.e. addition of source values, multiplication with kernel values,
    and NumericTraits must be defined. 
    The kernel's value_type must be an algebraic field,
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must
    be defined.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class SrcIterator, class SrcAccessor, 
	      class DestIterator, class DestAccessor,
	      class KernelIterator, class KernelAccessor>
    void convolveImage(SrcIterator sul, SrcIterator slr, SrcAccessor as,
		       DestIterator dul, DestAccessor ad,
		       KernelIterator ki, KernelAccessor ak, 
		       Diff2D kul, Diff2D klr, BorderTreatmentMode border)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class SrcIterator, class SrcAccessor, 
	      class DestIterator, class DestAccessor,
	      class KernelIterator, class KernelAccessor>
    void convolveImage(
       triple<SrcIterator, SrcIterator, SrcAccessor> src,
       pair<DestIterator, DestAccessor> dest,
       tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D, 
                                         BorderTreatmentMode> kernel)
    \end{verbatim}
    
    {\bf Usage:}
    
    Include-File:
    \URL[vigra/stdconvolution.hxx]{../include/vigra/stdconvolution.hxx}
    
    
    \begin{verbatim}
    FImage src(w,h), dest(w,h);    
    ...
    
    // define horizontal Sobel filter
    Kernel2D<float> sobel;
    
    sobel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =  // upper left and lower right
                         0.125, 0.0, -0.125,
			 0.25,  0.0, -0.25,
			 0.125, 0.0, -0.125;
        
    convolveImage(srcImageRange(src), destImage(dest), kernel2d(sobel));	     
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    ImageIterator sul, slr;
    ImageIterator dul;
    ImageIterator ik;
    
    SrcAccessor src_accessor;
    DestAccessor dest_accessor;
    KernelAccessor kernel_accessor;
    
    NumericTraits<SrcAccessor::value_type>::RealPromote s = src_accessor(sul);

    s = s + s;
    s = kernel_accessor(ik) * s;
    s -= s;

    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(s), dul);

    NumericTraits<KernelAccessor::value_type>::RealPromote k = kernel_accessor(ik);
    
    k += k;
    k -= k;
    k = k / k;

    \end{verbatim}

    {\bf Preconditions:}
    
    \begin{verbatim}
    kul.x <= 0
    kul.y <= 0
    klr.x >= 0
    klr.y >= 0
    slr.x - sul.x >= klr.x + kul.x + 1
    slr.y - sul.y >= klr.y + kul.y + 1
    border == BORDER_TREATMENT_CLIP || border == BORDER_TREATMENT_AVOID
    \end{verbatim}

    If border == BORDER_TREATMENT_CLIP: Sum of kernel elements must be
    != 0.

    @memo
*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
void convolveImage(SrcIterator sul, SrcIterator slr, SrcAccessor as,
                   DestIterator dul, DestAccessor ad,
		   KernelIterator ki, KernelAccessor ak, 
		   Diff2D kul, Diff2D klr, BorderTreatmentMode border)
{
    precondition((border == BORDER_TREATMENT_CLIP  ||
		  border == BORDER_TREATMENT_AVOID),
		  "convolveImage(): "
		  "Border treatment must be BORDER_TREATMENT_CLIP or BORDER_TREATMENT_AVOID.");

    precondition(kul.x <= 0 && kul.y <= 0,
		 "convolveImage(): coordinates of "
                 "kernel's upper left must be <= 0.");
    precondition(klr.x >= 0 && klr.y >= 0,
		 "convolveImage(): coordinates of "
                 "kernel's lower right must be >= 0.");

    // use traits to determine SumType as to prevent possible overflow
    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote SumType;
    typedef typename 
        NumericTraits<typename KernelAccessor::value_type>::RealPromote KSumType;
    
    // calculate width and height of the image
    int w = slr.x - sul.x;  
    int h = slr.y - sul.y;  
    int wk = klr.x - kul.x + 1;  
    int hk = klr.y - kul.y + 1;  
    
    precondition(w >= wk && h >= hk,
		 "convolveImage(): kernel larger than image.");

    int x,y;
    int ystart = (border == BORDER_TREATMENT_AVOID) ?  klr.y : 0;
    int yend   = (border == BORDER_TREATMENT_AVOID) ? h+kul.y : h;
    int xstart = (border == BORDER_TREATMENT_AVOID) ?  klr.x : 0;
    int xend   = (border == BORDER_TREATMENT_AVOID) ? w+kul.x : w;
    
    // create y iterators
    DestIterator yd = dul + Diff2D(xstart, ystart);
    SrcIterator ys = sul + Diff2D(xstart, ystart);
    
    KSumType norm = ak(ki);
    int xx, yy;
    KernelIterator yk  = ki + klr;
    for(yy=0; yy<hk; ++yy, --yk.y)
    {
	KernelIterator xk  = yk;
	
	for(xx=0; xx<wk; ++xx, --xk.x)
	{
	    norm += ak(xk);
	}
    }
    norm -= ak(ki);
    
    
    for(y=ystart; y < yend; ++y, ++ys.y, ++yd.y)
    {
        // create x iterators 
        DestIterator xd(yd);
        SrcIterator xs(ys);

        for(x=xstart; x < xend; ++x, ++xs.x, ++xd.x)
        {
            // how much of the kernel fits into the image ?
	    bool clip = false;
	    
	    int x0, y0, x1, y1;
	    
	    y0 = (y<klr.y) ? clip = true, -y : -klr.y;
	    y1 = (h-y-1<-kul.y) ? clip = true, h-y-1 : -kul.y;	    
	    x0 = (x<klr.x) ? clip = true, -x : -klr.x;
	    x1 = (w-x-1<-kul.x) ? clip = true, w-x-1 : -kul.x;
	    
	    if(!clip)
	    {
		// init the sum
		SumType sum = as(xs);
		
		SrcIterator yys = xs - klr;
		KernelIterator yk  = ki + klr;
		
		int xx, yy;
		for(yy=0; yy<hk; ++yy, ++yys.y, --yk.y)
		{
		    SrcIterator xxs = yys;
		    KernelIterator xk  = yk;
		    
		    for(xx=0; xx<wk; ++xx, ++xxs.x, --xk.x)
		    {
			sum += ak(xk) * as(xxs);
		    }
		}
		
		sum -= as(xs);
		// store average in destination pixel
		ad.set(sum, xd);
            }
	    else
	    {
		// init the sum
		SumType sum = as(xs);
		KSumType ksum = ak(ki);
		
		SrcIterator yys = xs + Diff2D(x0, y0);
		KernelIterator yk  = ki - Diff2D(x0, y0);
		
		int xx, yy, wk, hk;
		wk = x1 - x0 + 1;
		hk = y1 - y0 + 1;
		for(yy=0; yy<hk; ++yy, ++yys.y, --yk.y)
		{
		    SrcIterator xxs = yys;
		    KernelIterator xk  = yk;
		    
		    for(xx=0; xx<wk; ++xx, ++xxs.x, --xk.x)
		    {
			sum += ak(xk) * as(xxs);
			ksum += ak(xk);
		    }
		}
		
		sum -= as(xs);
		ksum -= ak(ki);
		
		// store average in destination pixel
		ad.set((norm / ksum) * sum, xd);
	    }
	}
    }    
}

   
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class KernelIterator, class KernelAccessor>
inline
void convolveImage(
   triple<SrcIterator, SrcIterator, SrcAccessor> src,
   pair<DestIterator, DestAccessor> dest,
   tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D, 
          BorderTreatmentMode> kernel)
{
    convolveImage(src.first, src.second, src.third,
                  dest.first, dest.second,
		  kernel.first, kernel.second, kernel.third, 
		  kernel.fourth, kernel.fifth);
}


/** Performs a 2 dimensional convolution of the source image within the 
    given ROI mask using the given kernel. The ROI is applied as follows: 
    Only pixel under the ROI are used in the calculations. Whenever a part of the 
    kernel lies outside the ROI, the kernel is renormalized to its original
    norm (analogous to the CLIP \Ref{BorderTreatmentMode}). An convolution result is 
    calculated whenever at the current kernel position {\it at least one pixel of the 
    kernel is within the ROI}. I.e., pixels not under the ROI may nevertheless
    be assigned a value if they are {\it near} the ROI. Thus, this algorithm is also
    useful as an interpolator. To get rid of the results outside the ROI mask, a 
    subsequent \Ref{copyImageIf}() must be performed. 
    
    The KernelIterator must point to the center of the kernel, and
    the kernel's size is given by its upper left (x and y of distance <= 0) and 
    lower right (distance >= 0) corners. The image must always be larger than the 
    kernel. At those positions where the kernel does not completely fit 
    into the image, the specified \Ref{BorderTreatmentMode} is 
    applied. Only BORDER_TREATMENT_CLIP and BORDER_TREATMENT_AVOID are currently 
    supported.
    
    The images's pixel type (SrcAccessor::value_type) must be a
    linear space over the kernel's value_type (KernelAccessor::value_type),
    i.e. addition of source values, multiplication with kernel values,
    and NumericTraits must be defined. 
    The kernel's value_type must be an algebraic field,
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must
    be defined.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class SrcIterator, class SrcAccessor, 
	      class MaskIterator, class MaskAccessor,
	      class DestIterator, class DestAccessor,
	      class KernelIterator, class KernelAccessor>
    void 
    convolveImageWithMask(SrcIterator sul, SrcIterator slr, SrcAccessor as,
			  MaskIterator mul, MaskAccessor am,
			  DestIterator dul, DestAccessor ad,
			  KernelIterator ki, KernelAccessor ak, 
			  Diff2D kul, Diff2D klr, BorderTreatmentMode border)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class SrcIterator, class SrcAccessor, 
	      class MaskIterator, class MaskAccessor,
	      class DestIterator, class DestAccessor,
	      class KernelIterator, class KernelAccessor>
    inline
    void convolveImageWithMask(
       triple<SrcIterator, SrcIterator, SrcAccessor> src,
       pair<MaskIterator, MaskAccessor> mask,
       pair<DestIterator, DestAccessor> dest,
       tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D, 
    	                                  BorderTreatmentMode> kernel)
    \end{verbatim}
    
    {\bf Usage:}
    
    Include-File:
    \URL[vigra/stdconvolution.hxx]{../include/vigra/stdconvolution.hxx}
    
    
    \begin{verbatim}
    FImage src(w,h), dest(w,h); 
    CImage mask(w,h);   
    ...
    
    // define 3x3 binomial filter
    Kernel2D<float> binom;
    
    binom.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =   // upper left and lower right
                         0.0625, 0.125, 0.0625,
			 0.125,  0.25,  0.125,
			 0.0625, 0.125, 0.0625;
        
    convolveImage(srcImageRange(src), maskImage(mask), destImage(dest), kernel2d(binom));	     
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    ImageIterator sul, slr;
    ImageIterator mul;
    ImageIterator dul;
    ImageIterator ik;
    
    SrcAccessor src_accessor;
    MaskAccessor mask_accessor;
    DestAccessor dest_accessor;
    KernelAccessor kernel_accessor;
    
    NumericTraits<SrcAccessor::value_type>::RealPromote s = src_accessor(sul);

    s = s + s;
    s = kernel_accessor(ik) * s;
    s -= s;

    if(mask_accessor(mul)) ...;
    
    dest_accessor.set(
        NumericTraits<DestAccessor::value_type>::fromRealPromote(s), dul);

    NumericTraits<KernelAccessor::value_type>::RealPromote k = kernel_accessor(ik);
    
    k += k;
    k -= k;
    k = k / k;

    \end{verbatim}

    {\bf Preconditions:}
    
    \begin{verbatim}
    kul.x <= 0
    kul.y <= 0
    klr.x >= 0
    klr.y >= 0
    slr.x - sul.x >= klr.x + kul.x + 1
    slr.y - sul.y >= klr.y + kul.y + 1
    border == BORDER_TREATMENT_CLIP || border == BORDER_TREATMENT_AVOID
    \end{verbatim}

    Sum of kernel elements must be != 0.

    @memo
*/
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class MaskIterator, class MaskAccessor,
          class KernelIterator, class KernelAccessor>
void 
convolveImageWithMask(SrcIterator sul, SrcIterator slr, SrcAccessor as,
                      MaskIterator mul, MaskAccessor am,
                      DestIterator dul, DestAccessor ad,
	   	      KernelIterator ki, KernelAccessor ak, 
		      Diff2D kul, Diff2D klr, BorderTreatmentMode border)
{
    precondition((border == BORDER_TREATMENT_CLIP  ||
		  border == BORDER_TREATMENT_AVOID),
		  "convolveImageWithMask(): "
		  "Border treatment must be BORDER_TREATMENT_CLIP or BORDER_TREATMENT_AVOID.");

    precondition(kul.x <= 0 && kul.y <= 0,
		 "convolveImageWithMask(): left borders must be <= 0.");
    precondition(klr.x >= 0 && klr.y >= 0,
		 "convolveImageWithMask(): right borders must be >= 0.");

    // use traits to determine SumType as to prevent possible overflow
    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::RealPromote SumType;
    typedef typename 
        NumericTraits<typename KernelAccessor::value_type>::RealPromote KSumType;
    
    // calculate width and height of the image
    int w = slr.x - sul.x;  
    int h = slr.y - sul.y;  
    int wk = klr.x - kul.x + 1;  
    int hk = klr.y - kul.y + 1;  
    
    int x,y;
    int ystart = (border == BORDER_TREATMENT_AVOID) ?  klr.y : 0;
    int yend   = (border == BORDER_TREATMENT_AVOID) ? h+kul.y : h;
    int xstart = (border == BORDER_TREATMENT_AVOID) ?  klr.x : 0;
    int xend   = (border == BORDER_TREATMENT_AVOID) ? w+kul.x : w;
    
    // create y iterators
    DestIterator yd = dul + Diff2D(xstart, ystart);
    SrcIterator ys = sul + Diff2D(xstart, ystart);
    MaskIterator ym = mul + Diff2D(xstart, ystart);
    
    KSumType norm = ak(ki);
    int xx, yy;
    KernelIterator yk  = ki + klr;
    for(yy=0; yy<hk; ++yy, --yk.y)
    {
	KernelIterator xk  = yk;
	
	for(xx=0; xx<wk; ++xx, --xk.x)
	{
	    norm += ak(xk);
	}
    }
    norm -= ak(ki);
    
    
    for(y=ystart; y < yend; ++y, ++ys.y, ++yd.y, ++ym.y)
    {
        // create x iterators 
        DestIterator xd(yd);
        SrcIterator xs(ys);
        MaskIterator xm(ym);

        for(x=xstart; x < xend; ++x, ++xs.x, ++xd.x, ++xm.x)
        {
            // how much of the kernel fits into the image ?
	    int x0, y0, x1, y1;
	    
	    y0 = (y<klr.y) ? -y : -klr.y;
	    y1 = (h-y-1<-kul.y) ? h-y-1 : -kul.y;	    
	    x0 = (x<klr.x) ? -x : -klr.x;
	    x1 = (w-x-1<-kul.x) ? w-x-1 : -kul.x;
	    	    
	    bool first = true;
	    // init the sum
	    SumType sum;
	    KSumType ksum;
	    
	    SrcIterator yys = xs + Diff2D(x0, y0);
	    MaskIterator yym = xm + Diff2D(x0, y0);
	    KernelIterator yk  = ki - Diff2D(x0, y0);
	    
	    int xx, yy, wk, hk;
	    wk = x1 - x0 + 1;
	    hk = y1 - y0 + 1;
	    for(yy=0; yy<hk; ++yy, ++yys.y, --yk.y, ++yym.y)
	    {
		SrcIterator xxs = yys;
		MaskIterator xxm = yym;
		KernelIterator xk  = yk;
		
		for(xx=0; xx<wk; ++xx, ++xxs.x, --xk.x, ++xxm.x)
		{
		    if(!am(xxm)) continue;
		    
		    if(first)
		    {
			sum = ak(xk) * as(xxs);
			ksum = ak(xk);
			first = false;
		    }
		    else
		    {
			sum += ak(xk) * as(xxs);
			ksum += ak(xk);
		    }
		}
	    }
	    // store average in destination pixel
	    if(!first && 
	        ksum != NumericTraits<KSumType>::zero()) 
	    {
	        ad.set((norm / ksum) * sum, xd);
	    }
	}
    }    
}

   
template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor,
          class MaskIterator, class MaskAccessor,
          class KernelIterator, class KernelAccessor>
inline
void convolveImageWithMask(
   triple<SrcIterator, SrcIterator, SrcAccessor> src,
   pair<MaskIterator, MaskAccessor> mask,
   pair<DestIterator, DestAccessor> dest,
   tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D, 
          BorderTreatmentMode> kernel)
{
    convolveImageWithMask(src.first, src.second, src.third,
                  mask.first, mask.second,
                  dest.first, dest.second,
		  kernel.first, kernel.second, kernel.third, 
		  kernel.fourth, kernel.fifth);
}

//@}

/********************************************************/
/*                                                      */
/*                      Kernel2D                        */
/*                                                      */
/********************************************************/

/** Generic 2 dimensional convolution kernel.
    This kernel may be used for convolution of 2 dimensional signals. 
    
    Convolution functions access the kernel via an ImageIterator
    which they get by calling \Ref{center}(). This iterator
    points to the center of the kernel. The kernel's size is given by its upperLeft() 
    (upperLeft().x <= 0, upperLeft().y <= 0) 
    and lowerRight() (lowerRight().x >= 0, lowerRight().y >= 0) methods. 
    The desired border treatment mode is returned by borderTreatment().
    (Note that the \Ref{Standard and masked convolution functions} don't currently 
    support all modes.)
    
    The different init functions create a kernel with the specified
    properties. The requirements for the kernel's value_type depend 
    on the init function used. At least NumericTraits must be defined.
    
    The kernel defines a factory function \Ref{kernel2d()} to create an argument object
    (see \Ref{Kernel Argument Object Factories}).
    
    {\bf Usage:}
    
    Include-File:
    \URL[vigra/stdconvolution.hxx]{../include/vigra/stdconvolution.hxx}
    
    \begin{verbatim}
    FImage src(w,h), dest(w,h);    
    ...
    
    // define horizontal Sobel filter
    Kernel2D<float> sobel;
    
    sobel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =  // upper left and lower right
                         0.125, 0.0, -0.125,
			 0.25,  0.0, -0.25,
			 0.125, 0.0, -0.125;
        
    convolveImage(srcImageRange(src), destImage(dest), kernel2d(sobel));	     
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    value_type v = NumericTraits<value_type>::one();
    \end{verbatim}

    See also the init functions.
    
*/

template <class ARITHTYPE>
class Kernel2D
{
  public:
        /** the kernel's value type
	    @memo
	*/
    typedef ARITHTYPE value_type;
    
        /** 2D random access iterator over the kernel's values
	    @memo
	*/
    typedef typename BasicImage<value_type>::Iterator Iterator;
    
        /** the kernel's accessor
	    @memo
	*/
    typedef typename BasicImage<value_type>::Accessor Accessor;
    
    struct InitProxy
    {
        typedef typename 
	   BasicImage<value_type>::ScanOrderIterator Iterator;
	   
        InitProxy(Iterator i, int count, value_type & norm)
	: iter_(i), base_(i),
	  count_(count), sum_(count),
	  norm_(norm)
	{}
	
	~InitProxy()
	{
	    precondition(count_ == 1 || count_ == sum_,
	          "Kernel2D::initExplicitly(): "
		  "Too few init values.");
	}
	
	InitProxy & operator,(value_type const & v)
	{
	    if(count_ == sum_)  norm_ = *iter_;
	    
	    --count_;
	    precondition(count_ > 0,
	          "Kernel2D::initExplicitly(): "
		  "Too many init values.");
		  
	    norm_ += v;
	    
	    ++iter_;
	    *iter_ = v;
	    
	    return *this;
	}
	
	Iterator iter_, base_;
	int count_, sum_;
	value_type & norm_;
    };
    
    static value_type one() { return NumericTraits<value_type>::one(); }
    
        /** Default constructor.
	    Creates a kernel of size 1x1 which would copy the signal
	    unchanged.
	    @memo
	*/
    Kernel2D()
    : kernel_(1, 1, Kernel2D<ARITHTYPE>::one()),
      left_(0, 0),
      right_(0, 0),
      border_treatment_(BORDER_TREATMENT_CLIP)
    {}
    
        /** Copy constructor.
	    @memo
	*/
    Kernel2D(Kernel2D const & k)
    : left_(k.left_),
      right_(k.right_),
      norm_(k.norm_),
      kernel_(k.kernel_),
      border_treatment_(k.border_treatment_)
    {}
    
        /** Copy assignment.
	    @memo
	*/
    Kernel2D & operator=(Kernel2D const & k)
    {
        if(this != &k)
	{
	    left_ = k.left_;
	    right_ = k.right_;
	    norm_ = k.norm_;
	    kernel_ = k.kernel_;
	}
	return *this;
    }
    
        /** Initialisation. 
	    This initializes the kernel with the given constant. The norm becomes
	    v*width()*height().
	    
	    Instead of a single value an initializer list of length width()*height() 
	    can be used like this:
	    
	    \begin{verbatim}
	    Kernel2D<float> binom;
	    
	    binom.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =  
				 0.0625, 0.125, 0.0625,
				 0.125,  0.25,  0.125,
				 0.0625, 0.125, 0.0625;
	    \end{verbatim}
	    
	    In this case, the norm will be set to the sum of the init values. 
	    An initializer list of wrong length will result in a run-time error.
	    @memo
	*/
    InitProxy operator=(value_type const & v)
    {
	int size = (right_.x - left_.x + 1) * 
	           (right_.y - left_.y + 1);
	kernel_ = v;
	norm_ = (double)size*v;
	
	return InitProxy(kernel_.begin(), size, norm_);
    }
    
        /** Destructor.
	    @memo
	*/
    ~Kernel2D() 
    {}
    
        /** Init the 2D kernel as the cartesian product of two 1D kernels
	    of type \Ref{Kernel1D}. The norm becomes the product of the two original
	    norms.
	    
	    {\bf Required Interface:}
	    
	    The kernel's value_type must be a linear algebra.
	    
	    \begin{verbatim}
	    Kernel2D<...>::value_type v;
	    v = v * v;
	    \end{verbatim}
	*/
    void initSeparable(Kernel1D<value_type> & kx, 
                       Kernel1D<value_type> & ky)
    {
        left_ = Diff2D(kx.left(), ky.left());
	right_ = Diff2D(kx.right(), ky.right());
	int w = right_.x - left_.x + 1;
	int h = right_.y - left_.y + 1;
	kernel_.resize(w, h);
	
	norm_ = kx.norm() * ky.norm();
	
	typedef typename Kernel1D<value_type>::Iterator KIter;
	typename Kernel1D<value_type>::Accessor ka;
	
	KIter kiy = ky.center() + left_.y;
	Iterator iy = center() + left_;
	
	for(int y=left_.y; y<=right_.y; ++y, ++kiy, ++iy.y)
	{
	    KIter kix = kx.center() + left_.x;
	    Iterator ix = iy;
	    for(int x=left_.x; x<=right_.x; ++x, ++kix, ++ix.x)
	    {
	        *ix = ka(kix) * ka(kiy);
	    }
	}
    }
    
        /** Init the 2D kernel as the cartesian product of two 1D kernels
	    given explicitly by iterators and sizes. The norm becomes the 
	    sum of the resulting kernel values.
		    
	    {\bf Required Interface:}
	    
	    The kernel's value_type must be a linear algebra.
	    
	    \begin{verbatim}
	    Kernel2D<...>::value_type v;
	    v = v * v;
	    v += v;
	    \end{verbatim}
	    
	    {\bf Preconditions:}
	    
	    \begin{verbatim}
	    xleft <= 0;
	    xright >= 0;
	    yleft <= 0;
	    yright >= 0;	    
	    \end{verbatim}
	*/
    template <class KernelIterator>
    void initSeparable(KernelIterator kxcenter, int xleft, int xright,
                       KernelIterator kycenter, int yleft, int yright)
    {
   	precondition(xleft <= 0 && yleft <= 0,
		     "Kernel2D::initSeparable(): left borders must be <= 0.");
	precondition(xright >= 0 && yright >= 0,
		     "Kernel2D::initSeparable(): right borders must be >= 0.");
    
        left_ = Diff2D(xleft, yleft);
	right_ = Diff2D(xright, yright);
	
	int w = right_.x - left_.x + 1;
	int h = right_.y - left_.y + 1;
	kernel_.resize(w, h);
	
	KernelIterator kiy = kycenter + left_.y;
	Iterator iy = center() + left_;
	
	for(int y=left_.y; y<=right_.y; ++y, ++kiy, ++iy.y)
	{
	    KernelIterator kix = kxcenter + left_.x;
	    Iterator ix = iy;
	    for(int x=left_.x; x<=right_.x; ++x, ++kix, ++ix.x)
	    {
	        *ix = *kix * *kiy;
	    }
	}

	typename BasicImage<value_type>::SequentialIterator i = kernel_.begin();
	typename BasicImage<value_type>::SequentialIterator iend = kernel_.end();
	norm_ = *i;
	++i;
	
	for(; i!= iend; ++i)  
	{
	    norm_ += *i;
	}
    }
      
        /** Init the 2D kernel as a circular averaging filter. The norm will be 
	    calculated as 
	    #NumericTraits<value_type>::one() / (number of non-zero kernel values)#.
	    The kernel's value_type must be a linear space.
		
	    {\bf Required Interface:}
	    
	    \begin{verbatim}
	    value_type v = NumericTraits<value_type>::one(); 

	    double d; 
	    v = d * v; 
	    \end{verbatim}
	    
	    {\bf Precondition:}
	    
	    \begin{verbatim}
	    radius > 0;
	    \end{verbatim}	
	*/
    void initDisk(int radius)
    {
   	precondition(radius > 0,
		     "Kernel2D::initDisk(): radius must be > 0.");
    
        left_ = Diff2D(-radius, -radius);
	right_ = Diff2D(radius, radius);
	int w = right_.x - left_.x + 1;
	int h = right_.y - left_.y + 1;
	kernel_.resize(w, h);
	norm_ = NumericTraits<value_type>::one();
	
	kernel_ = NumericTraits<value_type>::zero();
	double count = 0.0;
	
	Iterator k = center();
	double r2 = (double)radius*radius;
	
	int i;
	for(i=0; i<= radius; ++i)
	{
	    double r = (double) i - 0.5;
	    int w = (int)(sqrt(r2 - r*r) + 0.5);
	    for(int j=-w; j<=w; ++j)
	    {
	        k(j, i) = NumericTraits<value_type>::one();
	        k(j, -i) = NumericTraits<value_type>::one();
		count += (i != 0) ? 2.0 : 1.0;
	    }
	}
	
	count = 1.0 / count;

	for(int y=-radius; y<=radius; ++y)
	{
	    for(int x=-radius; x<=radius; ++x)
	    {
	        k(x,y) = count * k(x,y);
	    }
	}
    }
    
        /** Init the kernel by an explicit initializer list.
	    The upper left and lower right corners of the kernel must be passed.
	    A comma-separated initializer list is given after the assignment operator.
	    This function is used like this:
		
	    \begin{verbatim}
	    // define horizontal Sobel filter
	    Kernel2D<float> sobel;
	    
	    sobel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 
				 0.125, 0.0, -0.125,
				 0.25,  0.0, -0.25,
				 0.125, 0.0, -0.125;
	    \end{verbatim}
	    
	    The norm is set to the sum of the initialzer values. If the wrong number of 
	    values is given, a run-time error results. It is, however, possible to give 
	    just one initializer. This creates an averaging filter with the given constant:
	    
	    \begin{verbatim}
	    Kernel2D<float> average3x3;
	    
	    average3x3.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 1.0/9.0;
	    \end{verbatim}
	    
	    Here, the norm is set to value*width()*height().
		
	    {\bf Preconditions:}
	    
	    \begin{verbatim}
	    1. upperleft.x <= 0;
	    2. upperleft.y <= 0;
	    3. lowerright.x >= 0;
	    4. lowerright.y >= 0;
	    5. the number of values in the initializer list 
	       is 1 or equals the size of the kernel.
	    \end{verbatim}	
	*/
    Kernel2D & initExplicitly(Diff2D upperleft, Diff2D lowerright)
    {
   	precondition(upperleft.x <= 0 && upperleft.y <= 0,
		     "Kernel2D::initExplicitly(): left borders must be <= 0.");
	precondition(lowerright.x >= 0 && lowerright.y >= 0,
		     "Kernel2D::initExplicitly(): right borders must be >= 0.");
    
        left_ = upperleft;
	right_ = lowerright;
	
	int w = right_.x - left_.x + 1;
	int h = right_.y - left_.y + 1;
	kernel_.resize(w, h);
	
	return *this;
    }
    
        /** Coordinates of the upper left corner of the kernel.
	    @memo
	*/
    Diff2D upperLeft() const { return left_; }
    
        /** Coordinates of the lower right corner of the kernel.
	    @memo
	*/
    Diff2D lowerRight() const { return right_; }

        /** Width of the kernel.
	    @memo
	*/
    int width() const { return right_.x - left_.x + 1; }

        /** Height of the kernel.
	    @memo
	*/
    int height() const { return right_.y - left_.y + 1; }

        /** ImageIterator that points to the center of the kernel (coordinate (0,0)).
	    @memo
	*/
    Iterator center() { return kernel_.upperLeft() - left_; }
    
        /** Norm of the kernel (i.e. sum of its elements).
	    @memo
	*/
    value_type norm() const { return norm_; }
    
        /** The kernels default accessor.
	    @memo
	*/
    Accessor accessor() const { return Accessor(); }
    
        /** Normalize the kernel to the given value. (The norm is the sum of all kernel
	    elements.) The kernel's value_type must be a division algebra or 
	    algebraic field.
		
	    {\bf Required Interface:}
	    
	    \begin{verbatim}
	    value_type v = NumericTraits<value_type>::one(); // if norm is not 
	                                                     // given explicitly
							     
	    v += v; 
	    v = v * v; 
	    v = v / v; 
	    \end{verbatim}
	*/
    void normalize(value_type norm) 
    { 
	typename BasicImage<value_type>::SequentialIterator i = kernel_.begin();
	typename BasicImage<value_type>::SequentialIterator iend = kernel_.end();
	typename NumericTraits<value_type>::RealPromote sum = *i;
	++i;
	
	for(; i!= iend; ++i)  
	{
	    sum += *i;
	}
	
	sum = norm / sum;
	i = kernel_.begin();
	for(; i != iend; ++i)  
	{
	    *i = *i * sum;
	}
	
	norm_ = norm;
    }
    
        /** Normalize the kernel to norm 1.
	    @memo
	*/
    void normalize() 
    {
        normalize(one()); 
    }
    
        /** current border treatment mode
	    @memo
	*/
    BorderTreatmentMode borderTreatment() const 
    { return border_treatment_; }
    
        /** Set border treatment mode. 
	    Only #BORDER_TREATMENT_CLIP# and #BORDER_TREATMENT_AVOID# are currently
	    allowed.
	    @memo
	*/
    void setBorderTreatment( BorderTreatmentMode new_mode)
    { 
	precondition((new_mode == BORDER_TREATMENT_CLIP  ||
	              new_mode == BORDER_TREATMENT_AVOID),
		      "Kernel2D::setBorderTreatment(): "
		      "Border treatment must be BORDER_TREATMENT_CLIP or BORDER_TREATMENT_AVOID.");
    
        border_treatment_ = new_mode; 
    }
    
    
  private:
    BasicImage<value_type> kernel_;
    Diff2D left_, right_;
    value_type norm_;
    BorderTreatmentMode border_treatment_;
};

/**************************************************************/
/*                                                            */
/*         Argument object factories for Kernel2D             */
/*                                                            */
/*     (documentation: see vigra/convolution.hxx)             */
/*                                                            */
/**************************************************************/

template <class KernelIterator, class KernelAccessor>
inline
tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D, BorderTreatmentMode>
kernel2d(KernelIterator ik, KernelAccessor ak, Diff2D kul, Diff2D klr,
         BorderTreatmentMode border)

{
    return 
    tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D, BorderTreatmentMode> (
                                                           ik, ak, kul, klr, border);
}

template <class T>
inline
tuple5<typename Kernel2D<T>::Iterator, typename Kernel2D<T>::Accessor, 
       Diff2D, Diff2D, BorderTreatmentMode>
kernel2d(Kernel2D<T> & k)

{
    return 
	tuple5<typename Kernel2D<T>::Iterator, typename Kernel2D<T>::Accessor, 
	       Diff2D, Diff2D, BorderTreatmentMode>(
                                     k.center(), 
				     k.accessor(), 
				     k.upperLeft(), k.lowerRight(), 
				     k.borderTreatment());
}

template <class T>
inline
tuple5<typename Kernel2D<T>::Iterator, typename Kernel2D<T>::Accessor, 
       Diff2D, Diff2D, BorderTreatmentMode>
kernel2d(Kernel2D<T> & k, BorderTreatmentMode border)

{
    return 
	tuple5<typename Kernel2D<T>::Iterator, typename Kernel2D<T>::Accessor, 
	       Diff2D, Diff2D, BorderTreatmentMode>(
                                     k.center(), 
				     k.accessor(), 
				     k.upperLeft(), k.lowerRight(), 
				     border);
}



#endif // VIGRA_STDCONVOLUTION_HXX
