/************************************************************************/
/*                                                                      */
/*               Copyright 2011-2012 by Ullrich Koethe                  */
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

#ifndef VIGRA_HISTOGRAM_HXX
#define VIGRA_HISTOGRAM_HXX

#include "config.hxx"
#include "array_vector.hxx"
#include <algorithm>

namespace vigra {

template <class DataType, class BinType>
class HistogramView
{
    BinType * bins_;
    int size_, stride_;
    DataType offset_;
    double scale_, scaleInverse_;
    
  public:
    HistogramView(DataType const & min, DataType const & max, int binCount, 
                  BinType * bins = 0, int stride = 1)
    : bins_(bins),
      size_(binCount),
      stride_(stride),
      offset_(min),
      scale_(double(binCount) / (max - min)),
      scaleInverse_(1.0 / scale_)
    {}
    
    HistogramView & setData(BinType * bins , int stride = 1)
    {
        bins_ = bins;
        stride_ = stride;
        return *this;
    }
    
    HistogramView & reset()
    {
        if(hasData())
            for(int k=0; k<size_; ++k)
                *(bins_ +k*stride_) = BinType();
        return *this;
    }
    
    void getBinCenters(ArrayVector<DataType> * centers) const
    {
        double invScale = 1.0 / scale_;
        for(int k=0; k < size_; ++k)
        {
            (*centers)[k] = mapItemInverse(k + 0.5) ;
        }
    }
    
    int size() const
    {
        return size_;
    }
    
    bool hasData() const
    {
        return bins_ != 0;
    }
    
    BinType const & operator[](int k) const
    {
        return *(bins_ + k*stride_);
    }

	double mapItem(DataType const & d) const
	{
		return scale_ * (d - offset_);
	}
    
	DataType mapItemInverse(double d) const
	{
		return DataType(d * scaleInverse_ + offset_);
	}
    
    void add(DataType const & d, BinType weight = NumericTraits<BinType>::one())
    {
        get(int(mapItem(d))) += weight;
    }

  protected:

	BinType & get(int index)
	{
        if(index < 0)
            index = 0;
        if(index >= size_)
            index = size_ - 1;
		return *(bins_ + index*stride_);
	}
};

template <class T>
class TrapezoidKernel
{
  public:
    typedef T value_type;
    
    T operator[](double f) const
    {
        if(f < -0.5)
            return 0.5*(f + 1.5);
        if(f > 0.5)
            return 0.5*(1.5 - f);
        return 0.5;
    }
    
    double radius() const
    {
        return 1.5;
    }
    
    T findMaximum(double l, double c, double r) const
    {
        double curv = -2.0*c + r + l;
        if(curv == 0.0)
            return T(-0.5);
        double extr = 0.5*(l-r) / curv;
        if(curv < 0.0)
        {
            return zero < -0.5
                       ? T(-0.5)
                       : zero > 0.5
                              ? T(0.5)
                              : T(zero);
        }
        else
        {
            return zero < 0.0
                       ? T(0.5)
                       : T(-0.5);
        }
    }
    
    bool findMode(double l, double c, double r, double * m) const
    {
        double curv = -2.0*c + r + l;
        if(curv >= 0.0)
            return false;
        *m = 0.5*(l-r) / curv;
        if(*m < -0.5 || *m > 0.5)
            return false;
        return true;
    }
};

template <class DataType, class KernelType>
class KernelHistogramView
: public HistogramView<DataType, typename KernelType::value_type>
{
    KernelType kernel_;
    int radius_;
    
  public:
  
    typedef typename KernelType::value_type BinType;
    typedef HistogramView<DataType, BinType> BaseType;

    KernelHistogramView(DataType const & min, DataType const & max, int binCount, 
                        BinType * bins = 0, int stride = 1)
    : BaseType(min, max, binCount, bins, stride),
      radius_(kernel_.radius()-0.5) // FIXME: this needs generalization
    {}
    
    void add(DataType const & d, BinType weight = NumericTraits<BinType>::one())
    {
        double mapped = this->mapItem(d);
        double f = mapped - std::floor(mapped) - kernel_.radius();
        int center = int(mapped);
        
        for(int k=center+radius_; k>=center-radius_; --k, f += 1.0)
        {   
            this->get(k) += weight*kernel_[f];
		}
    }
    
    DataType findMode() const
    {
        double mmax = 0, vmax = 0, m;
        
        for(int k=0; k<this->size(); ++k)
        {
            double l = k > 0
                         ? (*this)[k-1]
                         : 0.0;
            double c = (*this)[k];
            double r = k < this->size() - 1
                         ? (*this)[k+1]
                         : 0.0;
            if(kernel_.findMode(l, c, r, &m))
            {
                double v = l*kernel_[m+1.0] + c*kernel_[m] + r*kernel_[m-1.0];
                if(vmax < v)
                {
                    mmax = m + k + 0.5;
                    vmax = v;
                }
            }
        }
        return mapItemInverse(mmax);
    }
    
    template <class Array>
    void findModes(Array * modes)
    {
        double m;
        for(int k=0; k<this->size(); ++k)
        {
            double l = k > 0
                         ? (*this)[k-1]
                         : 0.0;
            double c = (*this)[k];
            double r = k < this->size() - 1
                         ? (*this)[k+1]
                         : 0.0;
            if(kernel_.findMode(l, c, r, &m))
            {
                double v = l*kernel_[m+1.0] + c*kernel_[m] + r*kernel_[m-1.0];
                modes->push_back(std::make_pair(this->mapItemInverse(m + k + 0.5), v));
            }
        }
    }
};

template <class DataType, class BinType>
class Histogram
: public HistogramView<DataType, BinType>
{
  public:
    typedef HistogramView<DataType, BinType> BaseType;
    ArrayVector<BinType> data_;
    
  public:
    Histogram(DataType const & min, DataType const & max, int binCount, 
                  BinType * bins = 0, int stride = 1)
    : HistogramView(min, max, binCount),
      data_(binCount)
    {
        this->setData(&data_[0]);
    }
    
    Histogram const & reset()
    {
        this->setData(&data_[0]);
        BaseType::reset();
        return *this;
    }
 };

} // namespace vigra

#endif // VIGRA_HISTOGRAM_HXX
