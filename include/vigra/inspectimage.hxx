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
 
 
#ifndef VIGRA_INSPECTIMAGE_HXX
#define VIGRA_INSPECTIMAGE_HXX

#include <vector>
#include "vigra/utilities.hxx"
#include "vigra/numerictraits.hxx"

/** @name Functions to Inspect Images
    @memo apply read-only functor to every pixel
*/
//@{

/********************************************************/
/*                                                      */
/*                        inspectImage                  */
/*                                                      */
/********************************************************/

/** Apply read-only functor to every pixel in the image.
    This function can be used to collect statistics of the image etc.
    The results must be stored in the functor, which serves as a return
    value.
    The function uses an accessor to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class ImageIterator, class Accessor, class Functor>
    void
    inspectImage(ImageIterator upperleft, ImageIterator lowerright,
                 Accessor a, Functor & f)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class ImageIterator, class Accessor, class Functor>
    void
    inspectImage(triple<ImageIterator, ImageIterator, Accessor> img,
         Functor & f)
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/inspectimage.hxx]{../include/vigra/inspectimage.hxx}
    
    \begin{verbatim}
    BImage::Iterator upperleft = img.upperLeft();
    
    // init functor
    FindMinMax<BImage::PixelType> minmax();   
    
    inspectImage(srcImageRange(img), minmax);
    
    cout << "Min: " << minmax.min << " Max: " << minmax.max;

    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    ConstImageIterator upperleft, lowerright;
    Accessor accessor;
    Functor functor;
    
    functor(accessor(upperleft));         // return not used
    \end{verbatim}
    
    @memo
*/
template <class ImageIterator, class Accessor, class Functor>
void
inspectImage(ImageIterator upperleft, ImageIterator lowerright,
         Accessor a, Functor & f)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
    
    for(int y=0; y<h; ++y, ++upperleft.y)
    {
        ImageIterator ix(upperleft);
    
        for(int x=0; x<w; ++x, ++ix.x)
        {
            f(a(ix));
        }
    }
}
    
template <class ImageIterator, class Accessor, class Functor>
inline
void
inspectImage(triple<ImageIterator, ImageIterator, Accessor> img,
         Functor & f)
{
    inspectImage(img.first, img.second, img.third, f);
}

/********************************************************/
/*                                                      */
/*                      inspectImageIf                  */
/*                                                      */
/********************************************************/

/** Apply read-only functor to every pixel in the ROI.
    This function can be used to collect statistics of the roi etc.
    The functor is called whenever the return value of the mask's 
    accessor is not zero.
    The results must be stored in the functor, which serves as a return
    value.
    Accessors are used to access the pixel and mask data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class ImageIterator, class Accessor,
              class MaskImageIterator, class MaskAccessor, class Functor>
    void
    inspectImageIf(ImageIterator upperleft, ImageIterator lowerright, 
           MaskImageIterator mask_upperleft, MaskAccessor ma,
           Functor & f)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class ImageIterator, class Accessor,
          class MaskImageIterator, class MaskAccessor, class Functor>
    void
    inspectImageIf(triple<ImageIterator, ImageIterator, Accessor> img,
           pair<MaskImageIterator, MaskAccessor> mask,
           Functor & f)
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/inspectimage.hxx]{../include/vigra/inspectimage.hxx}
    
    \begin{verbatim}
    BImage img(100, 100);
    BImage mask(100, 100);
    BImage::Iterator upperleft = img.upperLeft();
    
    // init functor
    FindMinMax<BImage::PixelType> minmax();    
    
    inspectImageIf(srcImageRange(img), 
                   maskImage(mask), minmax);
    
    cout << "Min: " << minmax.min << " Max: " << minmax.max;

    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    ImageIterator upperleft, lowerright;
    MaskImageIterator mask_upperleft;
    
    Accessor accessor;
    MaskAccessor mask_accessor;
    
    Functor functor;

    if(mask_accessor(mask_upperleft)) functor(accessor(upperleft)); 
    \end{verbatim}
    
    @memo
*/
template <class ImageIterator, class Accessor,
      class MaskImageIterator, class MaskAccessor, class Functor>
void
inspectImageIf(ImageIterator upperleft, 
               ImageIterator lowerright, Accessor a,
           MaskImageIterator mask_upperleft, MaskAccessor ma,
           Functor & f)
{
    int w = lowerright.x - upperleft.x;
    int h = lowerright.y - upperleft.y;
            
    for(int y=0; y<h; ++y, ++upperleft.y, ++mask_upperleft.y)
    {
        ImageIterator ix(upperleft);
        MaskImageIterator mx(mask_upperleft);
    
        for(int x=0; x<w; ++x, ++ix.x, ++mx.x)
    {
        if(ma(mx)) f(a(ix));
    }
    }
}

template <class ImageIterator, class Accessor,
      class MaskImageIterator, class MaskAccessor, class Functor>
inline
void
inspectImageIf(triple<ImageIterator, ImageIterator, Accessor> img,
               pair<MaskImageIterator, MaskAccessor> mask,
               Functor & f)
{
    inspectImageIf(img.first, img.second, img.third, 
                   mask.first, mask.second, f);
}

/********************************************************/
/*                                                      */
/*                  inspectTwoImages                    */
/*                                                      */
/********************************************************/

/** Apply read-only functor to every pixel of both images.
    This function can be used to collect statistics for each region of a 
    labeled image, especially in conjunction with 
    the \Ref{ArrayOfRegionStatistics} functor. The results must be 
    stored in the functor which serves as a return value. 
    Accessors are used to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class ImageIterator1, class Accessor1,
          class ImageIterator2, class Accessor2, 
          class Functor>
    void
    inspectTwoImages(ImageIterator1 upperleft1, ImageIterator1 lowerright1, Accessor1 a1, 
             ImageIterator2 upperleft2, Accessor2 a2,
             Functor & f)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class ImageIterator1, class Accessor1,
          class ImageIterator2, class Accessor2, 
          class Functor>
    void
    inspectTwoImages(triple<ImageIterator1, ImageIterator1, Accessor1> img1,
                     pair<ImageIterator2, Accessor2> img2,
             Functor & f)
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/inspectimage.hxx]{../include/vigra/inspectimage.hxx}
    
    \begin{verbatim}
    BImage image1;
    BImage image2;
    
    SomeStatisticsFunctor stats(...);     // init functor
    
    inspectTwoImages(srcImageRange(image1), srcImage(image2),
                     region_stats);
    
    
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    ImageIterator1 upperleft1, lowerright1;
    ImageIterator2 upperleft2;
    
    Accessor1 accessor1;
    Accessor2 accessor2;
    
    Functor functor;
    functor(accessor1(upperleft1), accessor2(upperleft2));  // return not used
    \end{verbatim}
    
    @memo
*/
template <class ImageIterator1, class Accessor1,
          class ImageIterator2, class Accessor2, 
      class Functor>
void
inspectTwoImages(ImageIterator1 upperleft1, ImageIterator1 lowerright1, Accessor1 a1, 
                 ImageIterator2 upperleft2, Accessor2 a2,
         Functor & f)
{
    int w = lowerright1.x - upperleft1.x;
    int h = lowerright1.y - upperleft1.y;
    
    for(int y=0; y<h; ++y, ++upperleft1.y, ++upperleft2.y)
    {
        ImageIterator1 ix1(upperleft1);
        ImageIterator2 ix2(upperleft2);
    
        for(int x=0; x<w; ++x, ++ix1.x, ++ix2.x)
        {
            f(a1(ix1), a2(ix2));
        }
    }
}

template <class ImageIterator1, class Accessor1,
      class ImageIterator2, class Accessor2, 
      class Functor>
inline
void
inspectTwoImages(triple<ImageIterator1, ImageIterator1, Accessor1> img1,
         pair<ImageIterator2, Accessor2> img2,
         Functor & f)
{
    inspectTwoImages(img1.first, img1.second, img1.third, 
                     img2.first, img2.second, f);
}

/********************************************************/
/*                                                      */
/*                inspectTwoImagesIf                    */
/*                                                      */
/********************************************************/

/** Apply read-only functor to those pixels of both images where
    the mask image is non-zero.
    This function can be used to collect statistics for selected regions of a 
    labeled image, especially in conjunction with 
    the \Ref{ArrayOfRegionStatistics} functor. The results must be 
    stored in the functor which serves as a return value. 
    Accessors are used to access the pixel data.
    
    {\bf Declarations:}
    
    pass arguments explicitly:
    \begin{verbatim}
    template <class ImageIterator1, class Accessor1,
              class ImageIterator2, class Accessor2, 
              class MaskImageIterator, class MaskAccessor, 
          class Functor>
    void
    inspectTwoImagesIf(ImageIterator1 upperleft1, ImageIterator1 lowerright1, Accessor1 a1, 
                     ImageIterator2 upperleft2, Accessor2 a2,
                     MaskImageIterator mupperleft, MaskAccessor mask,
             Functor & f)
    \end{verbatim}
    
    
    use argument objects in conjuction with \Ref{Argument Object Factories}:
    \begin{verbatim}
    template <class ImageIterator1, class Accessor1,
          class ImageIterator2, class Accessor2, 
              class MaskImageIterator, class MaskAccessor, 
          class Functor>
    void
    inspectTwoImagesIf(triple<ImageIterator1, ImageIterator1, Accessor1> img1,
             pair<ImageIterator2, Accessor2> img2,
             pair<MaskImageIterator, MaskAccessor> mimg,
             Functor & f)
    \end{verbatim}
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/inspectimage.hxx]{../include/vigra/inspectimage.hxx}
    
    \begin{verbatim}
    BImage image1;
    BImage image2;
    BImage maskimage;
    
    SomeStatisticsFunctor stats(...);     // init functor
    
    inspectTwoImagesIf(srcImageRange(image1), srcImage(image2),
                       srcImage(maskimage), region_stats);
    
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    ImageIterator1 upperleft1, lowerright1;
    ImageIterator2 upperleft2;
    MaskImageIterator upperleftm;
    
    Accessor1 accessor1;
    Accessor2 accessor2;
    MaskAccessor mask;
    
    Functor functor;
    if(mask(upperleftm))
        functor(accessor1(upperleft1), accessor2(upperleft2)); 
    \end{verbatim}
    
    @memo
*/
template <class ImageIterator1, class Accessor1,
          class ImageIterator2, class Accessor2, 
          class MaskImageIterator, class MaskAccessor, 
      class Functor>
void
inspectTwoImagesIf(ImageIterator1 upperleft1, ImageIterator1 lowerright1, Accessor1 a1, 
                 ImageIterator2 upperleft2, Accessor2 a2,
                 MaskImageIterator mupperleft, MaskAccessor mask,
                 Functor & f)
{
    int w = lowerright1.x - upperleft1.x;
    int h = lowerright1.y - upperleft1.y;
    
    for(int y=0; y<h; ++y, ++upperleft1.y, ++upperleft2.y, ++mupperleft.y)
    {
        ImageIterator1 ix1(upperleft1);
        ImageIterator2 ix2(upperleft2);
        MaskImageIterator ixm(mupperleft);
    
        for(int x=0; x<w; ++x, ++ix1.x, ++ix2.x, ++ixm.x)
        {
            if(mask(ixm)) f(a1(ix1), a2(ix2));
        }
    }
}

template <class ImageIterator1, class Accessor1,
          class ImageIterator2, class Accessor2, 
          class MaskImageIterator, class MaskAccessor, 
          class Functor>
inline
void
inspectTwoImagesIf(triple<ImageIterator1, ImageIterator1, Accessor1> img1,
         pair<ImageIterator2, Accessor2> img2,
         pair<MaskImageIterator, MaskAccessor> m,
         Functor & f)
{
    inspectTwoImagesIf(img1.first, img1.second, img1.third, 
                     img2.first, img2.second, 
                     m.first, m.second, 
                     f);
}

//@}

/** @name Functors To Inspect Images
    @memo functors which report image statistics
*/
//@{

/********************************************************/
/*                                                      */
/*                     FindMinMax                       */
/*                                                      */
/********************************************************/

/** Find the minimum and maximum pixel value in an image or ROI.
    In addition the size of the ROI is calculated.
    This Functor can also be used in conjunction with
    \Ref{ArrayOfRegionStatistics} to find the extremes of all regions in
    a labeled image.
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/inspectimage.hxx]{../include/vigra/inspectimage.hxx}
    
    \begin{verbatim}
    BImage::Iterator img;
    
    FindMinMax<BImage::PixelType> minmax;   // init functor
    
    inspectImage(srcImageRange(img), minmax);
    
    cout << "Min: " << minmax.min << " Max: " << minmax.max;
    
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    VALUETYPE v1, v2(v1);
    
    v1 < v2;
    v1 = v2;
    \end{verbatim}
    
*/
template <class VALUETYPE>
class FindMinMax
{
   public:
   
        /** the functor's value type
        @memo
    */
    typedef VALUETYPE value_type;
   
    /** init min and max
    @memo
    */
    FindMinMax()
    : count(0)
    {}
    
    /** update min and max
        @memo
    */
    void operator()(VALUETYPE const & v) 
    {
        if(count)
        {
            if(v < min) min = v;
            if((this->max) < v) max = v;
        }
        else
        {
            min = v;
            max = v;
        }
        ++count;
    }
    
    /** merge two statistics
        @memo
    */
    void operator()(FindMinMax const & v) 
    {
        if(v.count)
        {
            if(count)
                {
                    if(v.min < min) min = v.min;
                if((this->max) < v.max) max = v.max;
                }
                else
                {
                min = v.min;
                max = v.max;
                }
        }
        count += v.count;
    }
    
    /** the current min
        @memo
    */
    VALUETYPE min;
    
    /** the current max
        @memo
    */
    VALUETYPE max;
    
    /** the number of values processed so far
        @memo
    */
    unsigned int count;
    
}; 

/********************************************************/
/*                                                      */
/*                    FindAverage                       */
/*                                                      */
/********************************************************/

/** Find the average pixel value in an image or ROI.
    In addition the size of the ROI is calculated.
    This Functor can also be used in conjunction with
    \Ref{ArrayOfRegionStatistics} to find the average of all regions in
    a labeled image.
    
    {\bf Usage:}
    
        Include-File:
        \URL[vigra/inspectimage.hxx]{../include/vigra/inspectimage.hxx}
    
    \begin{verbatim}
    BImage::Iterator img;
    
    FindAverage<BImage::PixelType> average;   // init functor
    
    inspectImage(srcImageRange(img), average);
    
    cout << "Average: " << average();
    
    \end{verbatim}

    {\bf Required Interface:}
    
    \begin{verbatim}
    VALUETYPE v1, v2(v1);
    
    v1 < v2;
    v1 = v2;
    \end{verbatim}
    
*/
template <class VALUETYPE>
class FindAverage
{
   public:
   
        /** the functor's value type
        @memo
    */
    typedef typename NumericTraits<VALUETYPE>::RealPromote value_type;
   
    /** init average
    @memo
    */
    FindAverage()
    : count(0), sum(NumericTraits<value_type>::zero())
    {}
    
    /** update average
        @memo
    */
    void operator()(VALUETYPE const & v) 
    {
        sum += v;
        ++count;
    }
    
    /** merge two statistics
        @memo
    */
    void operator()(FindAverage const & v) 
    {
        sum += v.sum;
        count += v.count;
    }
    
    /** return current average
        @memo
    */
    value_type average() const
    {
        return sum / (double)count;
    }
    
    /** return current average
        @memo
    */
    value_type operator()() const
    {
        return sum / (double)count;
    }
    
    value_type sum;
    unsigned int count;
    
}; 

/********************************************************/
/*                                                      */
/*                    FindROISize                       */
/*                                                      */
/********************************************************/

/** Calculate the size of an ROI in an image. 
    This Functor is often used in conjunction with
    \Ref{ArrayOfRegionStatistics} to find the sizes of all regions in
    a labeled image.
    
    {\bf Usage:}
    
    Include-File:
    \URL[vigra/inspectimage.hxx]{../include/vigra/inspectimage.hxx}

    \begin{verbatim}
    BImage img, mask;
    
    FindROISize<BImage::PixelType> roisize;   // init functor
    
    inspectImageIf(srcImageRange(img), srcImage(mask), roisize);
    
    cout << "Size of ROI: " << roisize.count;
    
    \end{verbatim}

*/
template <class VALUETYPE>
class FindROISize
{
   public:
   
        /** the functor's value type
        @memo
    */
    typedef VALUETYPE value_type;
   
    /** init counter to 0
    @memo
    */
    FindROISize()
    : count(0)
    {}
    
    /** update counter
        @memo
    */
    void operator()(VALUETYPE const &) 
    {
        ++count;
    }
    
    /** return current size
        @memo
    */
    int operator()() const 
    {
        return count;
    }
    
    /** return current size
        @memo
    */
    int size() const 
    {
        return count;
    }
    
    /** merge two statistics
        @memo
    */
    void operator()(FindROISize const & o) 
    {
        count += o.count;
    }
    
    /** the current counter
        @memo
    */
    unsigned int count;
    
}; 

/********************************************************/
/*                                                      */
/*                 LastValueFunctor                     */
/*                                                      */
/********************************************************/

/** Stores and returns the last value it has seen. 
    This Functor is best used in conjunction with
    \Ref{ArrayOfRegionStatistics} to realize a look-up table.
    
    {\bf Usage:}
    
    Include-File:
    \URL[vigra/inspectimage.hxx]{../include/vigra/inspectimage.hxx}

    \begin{verbatim}
    BImage img;
    
    ArrayOfRegionStatistics<LastValueFunctor<unsigned char> > lut(255);
    
    for(int i=0; i<256; ++i)
    {
        lut[i] = ...; // init look-up table
    }
    
    transformImage(srcImageRange(img), destImage(img), lut);
    
    \end{verbatim}

*/
template <class VALUETYPE>
class LastValueFunctor
{
   public:
   
        /** the functor's value type
        @memo
    */
    typedef VALUETYPE value_type;
   
    /** default initialization of value
        @memo
    */
    LastValueFunctor()
    {}

        /** replace value
        @memo
    */
    void operator=(VALUETYPE const & v) { value = v; }
        
    /** replace value
        @memo
    */
    void operator()(VALUETYPE const & v) { value = v; }
    
    /** return current value
        @memo
    */
    VALUETYPE operator()() const { return value; }

    /** the current value
        @memo
    */
    VALUETYPE value;
    
}; 

/********************************************************/
/*                                                      */
/*              ArrayOfRegionStatistics                 */
/*                                                      */
/********************************************************/

/** Calculate statistics for all regions of a labeled image. 
    This Functor encapsulates an array of statistics functors, one
    for each label, and selects the one to be updated according to the 
    pixel's label. 
    
    {\bf Usage:}
    
    Include-File:
    \URL[vigra/inspectimage.hxx]{../include/vigra/inspectimage.hxx}

    \begin{verbatim}
    BImage img;
    IImage labels;
    int max_label;
    ...
    
    // init functor as an array of 'max_label' FindMinMax-Functors
    ArrayOfRegionStatistics<FindMinMax<BImage::PixelType> > 
                                                         minmax(max_label);   
        
    inspectTwoImages(srcImageRange(img), srcImage(labels), minmax);
    
    for(int i=0; i<= max_label; ++i)
    {
        cout << "Max gray lavel of region " << i << ": " 
             << minmax.region[i].max << endl;
    }
    
    // init functor as an array of 'max_label' FindAverage-Functors
    ArrayOfRegionStatistics<FindAverage<BImage::PixelType> > 
                                                         average(max_label);   
        
    inspectTwoImages(srcImageRange(img), srcImage(labels), average);
    
    // write back the average of each region into the original image
    transformImage(srcImageRange(labels), destImage(img), average);
    
    \end{verbatim}
    
    {\bf Required Interface:}

    \begin{verbatim}
    RegionStatistics region;
    RegionStatistics::value_type v;
    
    region(v);     // update statistics
    v = region();  // return statistics
    
    \end{verbatim}
*/
template <class RegionStatistics>
class ArrayOfRegionStatistics
{
  public:
         ///
    typedef typename RegionStatistics::value_type value_type;
         ///
    typedef typename std::vector<RegionStatistics>::iterator iterator;

    /** init array of RegionStatistics with index domain
        0...max_region_label.
        @memo
    */
    ArrayOfRegionStatistics(int max_region_label)
    : regions(max_region_label+1)
    {}
    
    
    /** update regions statistics for region 'label'
        @memo
    */
    void operator()(value_type const & v, int label) {
        regions[label](v);
    }
    
    /** merge second region into first
        @memo
    */
    void merge(int label1, int label2) {
        regions[label1](regions[label2]);
    }
    
    /** ask for maximal index (label) allowed
        @memo
    */
    int maxRegionLabel() const { return size() - 1; } 

    /** ask for array size (i.e. maxRegionLabel() + 1)
        @memo
    */
    int size() const { return regions.size(); } 

    /** access the statistics for a region via its label
        @memo
    */
    value_type operator()(int label) { return regions[label](); }
    
    /** access the statistics functor for a region via its label
        @memo
    */
    RegionStatistics & operator[](int label) { return regions[label]; }
    
    /** iterator to the begin of the region array
        @memo
    */
    iterator begin() { return regions.begin(); } 
    
    /** iterator to the end of the region array
        @memo
    */
    iterator end() { return regions.end(); } 
    
  private:
    std::vector<RegionStatistics> regions;
};


//@}



#endif // VIGRA_INSPECTIMAGE_HXX
