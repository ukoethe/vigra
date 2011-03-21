/************************************************************************/
/*                                                                      */
/*       Copyright 2009 by Ullrich Koethe and Hans Meine                */
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

#ifndef VIGRA_NUMPY_ARRAY_TAGGEDSHAPE_HXX
#define VIGRA_NUMPY_ARRAY_TAGGEDSHAPE_HXX

#include <string>
#include "numpy_array_utilities.hxx"
#include "axistags.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*                     PyAxisTags                       */
/*                                                      */
/********************************************************/

// FIXME: right now, we implement this class using the standard
//        Python C-API only. It would be easier and more efficient 
//        to use boost::python here, but it would cause NumpyArray
//        to depend on boost, making it more difficult to use
//        NumpyArray in connection with other glue code generators.
class PyAxisTags
{
  public:
    typedef PyObject * pointer;
    
    python_ptr axistags;
    
    // FIXME: is it OK to always create a copy of the tags?
    PyAxisTags(python_ptr tags)
    {
        if(!tags)
            return;
        // FIXME: do a more elaborate type check here?
        if(!PySequence_Check(tags))
        {
            PyErr_SetString(PyExc_TypeError, 
                           "PyAxisTags(tags): tags argument must have type 'AxisTags'.");
            pythonToCppException(false);
        }
        python_ptr func(PyString_FromString("__copy__"), python_ptr::keep_count);
        axistags = python_ptr(PyObject_CallMethodObjArgs(tags, func.get(), NULL), 
                              python_ptr::keep_count);
    }
    
    PyAxisTags(int ndim, std::string const & order = "")
    {
        axistags = detail::defaultAxistags(ndim, order);
    }
    
    long size() const
    {
        return axistags
                   ? PySequence_Length(axistags)
                   : 0;
    }
    
    long channelIndex(long defaultVal) const
    {
        python_ptr key(PyString_FromString("channelIndex"), python_ptr::keep_count);
        return detail::getAttrLong(axistags, key, defaultVal);
    }

    long channelIndex() const
    {
        return channelIndex(size());
    }

    long majorNonchannelIndex(long defaultVal) const
    {
        static python_ptr key(PyString_FromString("majorNonchannelIndex"), python_ptr::keep_count);
        return detail::getAttrLong(axistags, key, defaultVal);
    }

    long majorNonchannelIndex() const
    {
        return majorNonchannelIndex(size());
    }

    void setChannelDescription(std::string const & description)
    {
        if(!axistags)
            return;
        python_ptr d(PyString_FromString(description.c_str()), python_ptr::keep_count);
        python_ptr func(PyString_FromString("setChannelDescription"), python_ptr::keep_count);
        python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), d.get(), NULL), 
                       python_ptr::keep_count);
        pythonToCppException(res);
    }

    double resolution(long index)
    {
        if(!axistags)
            return 0.0;
        python_ptr func(PyString_FromString("resolution"), python_ptr::keep_count);
        python_ptr i(PyInt_FromLong(index), python_ptr::keep_count);
        python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), i.get(), NULL), 
                       python_ptr::keep_count);
        pythonToCppException(res);
        if(!PyFloat_Check(res))
        {
            PyErr_SetString(PyExc_TypeError, "AxisTags.resolution() did not return float.");
            pythonToCppException(false);
        }
        return PyFloat_AsDouble(res);
    }
 
    void setResolution(long index, double resolution)
    {
        if(!axistags)
            return;
        python_ptr func(PyString_FromString("setResolution"), python_ptr::keep_count);
        python_ptr i(PyInt_FromLong(index), python_ptr::keep_count);
        python_ptr r(PyFloat_FromDouble(resolution), python_ptr::keep_count);
        python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), i.get(), r.get(), NULL), 
                       python_ptr::keep_count);
        pythonToCppException(res);
    }
 
    void scaleAxisResolution(long index, double factor)
    {
        if(!axistags)
            return;
        python_ptr func(PyString_FromString("scaleAxisResolution"), python_ptr::keep_count);
        python_ptr i(PyInt_FromLong(index), python_ptr::keep_count);
        python_ptr f(PyFloat_FromDouble(factor), python_ptr::keep_count);
        python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), i.get(), f.get(), NULL), 
                       python_ptr::keep_count);
        pythonToCppException(res);
    }
 
    void toFrequencyDomain(long index, int size, int sign = 1)
    {
        if(!axistags)
            return;
        python_ptr func(sign == 1
                           ? PyString_FromString("toFrequencyDomain")
                           : PyString_FromString("fromFrequencyDomain"), 
                        python_ptr::keep_count);
        python_ptr i(PyInt_FromLong(index), python_ptr::keep_count);
        python_ptr s(PyInt_FromLong(size), python_ptr::keep_count);
        python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), i.get(), s.get(), NULL), 
                       python_ptr::keep_count);
        pythonToCppException(res);
    }
 
    void fromFrequencyDomain(long index, int size)
    {
        toFrequencyDomain(index, size, -1);
    }
 
    ArrayVector<npy_intp> 
    permutationToNormalOrder(bool ignoreErrors = false) const
    {
        ArrayVector<npy_intp> res;
        getAxisPermutationImpl(res, "permutationToNormalOrder", ignoreErrors);
        return res;
    }

    ArrayVector<npy_intp> 
    permutationFromNormalOrder(bool ignoreErrors = false) const
    {
        ArrayVector<npy_intp> res;
        getAxisPermutationImpl(res, "permutationFromNormalOrder", ignoreErrors);
        return res;
    }
    
    void dropChannelAxis()
    {
        if(!axistags)
            return;
        python_ptr func(PyString_FromString("dropChannelAxis"), 
                               python_ptr::keep_count);
        python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                       python_ptr::keep_count);
        pythonToCppException(res);
    }
    
    void insertChannelAxis()
    {
        if(!axistags)
            return;
        python_ptr func(PyString_FromString("insertChannelAxis"), 
                               python_ptr::keep_count);
        python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                       python_ptr::keep_count);
        pythonToCppException(res);
    }
    
    operator pointer()
    {
        return axistags.get();
    }

    bool operator!() const
    {
        return !axistags;
    }

  private:

    void 
    getAxisPermutationImpl(ArrayVector<npy_intp> & res, const char * name, bool ignoreErrors) const
    {
        if(!axistags)
            return;
        
        python_ptr func(PyString_FromString(name), python_ptr::keep_count);
        python_ptr permutation(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                               python_ptr::keep_count);
        if(!permutation && ignoreErrors)
        {
            PyErr_Clear();
            return;
        }
        pythonToCppException(permutation);
        
        if(!PySequence_Check(permutation) || PySequence_Length(permutation) != size())
        {
            if(ignoreErrors)
                return;
            std::string message = std::string(name) + "() did not return a sequence, or a sequence of wrong length.";
            PyErr_SetString(PyExc_ValueError, message.c_str());
            pythonToCppException(false);
        }
            
        res.resize(size());
        for(int k=0; k<size(); ++k)
        {
            python_ptr i(PySequence_GetItem(permutation, k), python_ptr::keep_count);
            if(!PyInt_Check(i))
            {
                if(ignoreErrors)
                {
                    res.clear();
                    return;
                }
                std::string message = std::string(name) + "() did not return a sequence of int.";
                PyErr_SetString(PyExc_ValueError, message.c_str());
                pythonToCppException(false);
            }
            res[k] = PyInt_AsLong(i);
        }
    }
};

/********************************************************/
/*                                                      */
/*                     TaggedShape                      */
/*                                                      */
/********************************************************/

class TaggedShape
{
  public:
    enum ChannelAxis { first, last, none };
    
    ArrayVector<npy_intp> shape, original_shape;
    PyAxisTags axistags;
    ChannelAxis channelAxis;
    std::string channelDescription;
    
    TaggedShape(MultiArrayIndex size)
    : shape(size),
      axistags(size),
      channelAxis(none)
    {}
    
    template <class U, int N>
    TaggedShape(TinyVector<U, N> const & sh, PyAxisTags tags)
    : shape(sh.begin(), sh.end()),
      original_shape(sh.begin(), sh.end()),
      axistags(tags),
      channelAxis(none)
    {}
    
    template <class T>
    TaggedShape(ArrayVector<T> const & sh, PyAxisTags tags)
    : shape(sh.begin(), sh.end()),
      original_shape(sh.begin(), sh.end()),
      axistags(tags),
      channelAxis(none)
    {}
    
    template <class U, int N>
    TaggedShape(TinyVector<U, N> const & sh)
    : shape(sh.begin(), sh.end()),
      original_shape(sh.begin(), sh.end()),
      axistags(sh.size()),
      channelAxis(none)
    {}
    
    template <class T>
    TaggedShape(ArrayVector<T> const & sh)
    : shape(sh.begin(), sh.end()),
      original_shape(sh.begin(), sh.end()),
      axistags(sh.size()),
      channelAxis(none)
    {}
    
    template <class U, int N>
    TaggedShape & resize(TinyVector<U, N> const & sh)
    {
        int start = channelAxis == first
                        ? 1
                        : 0, 
            stop = channelAxis == last
                        ? (int)size()-1
                        : (int)size();
                        
        vigra_precondition(N == stop - start || size() == 0,
             "TaggedShape.resize(): size mismatch.");
             
        if(size() == 0)
            shape.resize(N);
        
        for(int k=0; k<N; ++k)
            shape[k+start] = sh[k];
            
        return *this;
    }
    
    TaggedShape & resize(MultiArrayIndex v1)
    {
        return resize(TinyVector<MultiArrayIndex, 1>(v1));
    }
    
    TaggedShape & resize(MultiArrayIndex v1, MultiArrayIndex v2)
    {
        return resize(TinyVector<MultiArrayIndex, 2>(v1, v2));
    }
    
    TaggedShape & resize(MultiArrayIndex v1, MultiArrayIndex v2, MultiArrayIndex v3)
    {
        return resize(TinyVector<MultiArrayIndex, 3>(v1, v2, v3));
    }
    
    TaggedShape & resize(MultiArrayIndex v1, MultiArrayIndex v2, 
                         MultiArrayIndex v3, MultiArrayIndex v4)
    {
        return resize(TinyVector<MultiArrayIndex, 4>(v1, v2, v3, v4));
    }
    
    npy_intp & operator[](int i)
    {
        return shape[i];
    }
    
    npy_intp operator[](int i) const
    {
        return shape[i];
    }
    
    unsigned int size() const
    {
        return shape.size();
    }
    
    TaggedShape & operator+=(int v)
    {
        int start = channelAxis == first
                        ? 1
                        : 0, 
            stop = channelAxis == last
                        ? (int)size()-1
                        : (int)size();
        for(int k=start; k<stop; ++k)
            shape[k] += v;
            
        return *this;
    }
    
    TaggedShape & operator-=(int v)
    {
        return operator+=(-v);
    }
    
    TaggedShape & operator*=(int factor)
    {
        int start = channelAxis == first
                        ? 1
                        : 0, 
            stop = channelAxis == last
                        ? (int)size()-1
                        : (int)size();
        for(int k=start; k<stop; ++k)
            shape[k] *= factor;
            
        return *this;
    }
    
    void rotateToNormalOrder()
    {
        if(channelAxis == last)
        {
            int ndim = (int)size();
            
            npy_intp channelCount = shape[ndim-1];            
            for(int k=ndim-1; k>0; --k)
                shape[k] = shape[k-1];
            shape[0] = channelCount;
            
            channelCount = original_shape[ndim-1];            
            for(int k=ndim-1; k>0; --k)
                original_shape[k] = original_shape[k-1];
            original_shape[0] = channelCount;
            
            channelAxis = first;
        }
    }
    
    TaggedShape & setChannelDescription(std::string const & description)
    {
        // we only remember the description here, and will actually set
        // it in the finalize function
        channelDescription = description;
        return *this;
    }
    
    TaggedShape & setChannelIndexLast()
    {
        // FIXME: add some checks?
        channelAxis = last;
        return *this;
    }
    
    // transposeShape() means: only shape and resolution are transposed, not the axis keys
    template <class U, int N>
    TaggedShape & transposeShape(TinyVector<U, N> const & p)
    {
        int ntags = axistags.size();
        ArrayVector<npy_intp> permute = axistags.permutationToNormalOrder();
        
        int tstart = (axistags.channelIndex(ntags) < ntags)
                        ? 1
                        : 0;
        int sstart = (channelAxis == first)
                        ? 1
                        : 0;
        int ndim = ntags - tstart;

        vigra_precondition(N == ndim,
             "TaggedShape.transposeShape(): size mismatch.");
             
        PyAxisTags newAxistags(axistags.axistags); // force copy
        for(int k=0; k<ndim; ++k)
        {
            original_shape[k+sstart] = shape[p[k]+sstart];
            newAxistags.setResolution(permute[k+tstart], axistags.resolution(permute[p[k]+tstart]));
        }
        shape = original_shape;
        axistags = newAxistags;
        
        return *this;
    }

    TaggedShape & toFrequencyDomain(int sign = 1)
    {
        int ntags = axistags.size();
        
        ArrayVector<npy_intp> permute = axistags.permutationToNormalOrder();
        
        int tstart = (axistags.channelIndex(ntags) < ntags)
                        ? 1
                        : 0;
        int sstart = (channelAxis == first)
                        ? 1
                        : 0;
        int size = (int)this->size() - sstart;
        
        for(int k=0; k<size; ++k)
        {
            axistags.toFrequencyDomain(permute[k+tstart], shape[k+sstart], sign);
        }
        
        return *this;
    }

    TaggedShape & fromFrequencyDomain()
    {
        return toFrequencyDomain(-1);
    }
    
    TaggedShape & setChannelCount(int count)
    {
        switch(channelAxis)
        {
          case first:
            shape[0] = count;
            break;
          case last:
            shape[size()-1] = count;
            break;
          case none:
            shape.push_back(count);
            original_shape.push_back(count);
            channelAxis = last;
            break;
        }
        return *this;
    }
};

// inline 
// void scaleAxisResolution(TaggedShape & tagged_shape)
// {
    // if(tagged_shape.size() != tagged_shape.original_shape.size())
        // return;
    
    // int ntags = PySequence_Length(tagged_shape.axistags);
    
    // ArrayVector<npy_intp> permute = detail::permutationToNormalOrder(tagged_shape.axistags);
    
    // int tstart = (detail::channelIndex(tagged_shape.axistags, ntags) < ntags)
                    // ? 1
                    // : 0;
    // int sstart = (tagged_shape.channelAxis == TaggedShape::first)
                    // ? 1
                    // : 0;
    // int size = (int)tagged_shape.size() - sstart;
    
    // for(int k=0; k<size; ++k)
    // {
        // int sk = k + sstart;
        // if(tagged_shape.shape[sk] == tagged_shape.original_shape[sk])
            // continue;
        // double factor = (tagged_shape.original_shape[sk] - 1.0) / (tagged_shape.shape[sk] - 1.0);
        // detail::scaleAxisResolution(tagged_shape.axistags, permute[k+tstart], factor);
    // }
// }

inline 
void scaleAxisResolution(TaggedShape & tagged_shape)
{
    if(tagged_shape.size() != tagged_shape.original_shape.size())
        return;
    
    int ntags = tagged_shape.axistags.size();
    
    ArrayVector<npy_intp> permute = tagged_shape.axistags.permutationToNormalOrder();
    
    int tstart = (tagged_shape.axistags.channelIndex(ntags) < ntags)
                    ? 1
                    : 0;
    int sstart = (tagged_shape.channelAxis == TaggedShape::first)
                    ? 1
                    : 0;
    int size = (int)tagged_shape.size() - sstart;
    
    for(int k=0; k<size; ++k)
    {
        int sk = k + sstart;
        if(tagged_shape.shape[sk] == tagged_shape.original_shape[sk])
            continue;
        double factor = (tagged_shape.original_shape[sk] - 1.0) / (tagged_shape.shape[sk] - 1.0);
        tagged_shape.axistags.scaleAxisResolution(permute[k+tstart], factor);
    }
}


// inline 
// ArrayVector<npy_intp> unifyTaggedShapeSize(TaggedShape & tagged_shape)
// {
    // python_ptr axistags = tagged_shape.axistags;
    // ArrayVector<npy_intp> shape = tagged_shape.shape;

    // if(!PySequence_Check(axistags))
    // {
        // PyErr_SetString(PyExc_TypeError, "constructArray(): axistags have wrong type.");
        // pythonToCppException(false);
    // }
    
    // int ndim = (int)shape.size();
    // int ntags = PySequence_Length(axistags);
    
    // long channelIndex = detail::channelIndex(axistags, ntags);

// #if 0 // debug only
    // std::cerr << "ndim: " << ndim << ", ntags: " << ntags << ", channelIndex: " << channelIndex << "\n";
    // static python_ptr func(PyString_FromString("__repr__"), 
                           // python_ptr::keep_count);
    // python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                   // python_ptr::keep_count);
    // pythonToCppException(res);
    // std::cerr << "axistags: " << PyString_AsString(res) << "\n";
// #endif

    // if(tagged_shape.channelAxis == TaggedShape::none)
    // {
        // // shape has no channel axis
        // if(channelIndex == ntags)
        // {
            // // axistags have no channel axis either => sizes should match
            // vigra_precondition(ndim == ntags,
                 // "constructArray(): size mismatch between shape and axistags.");
        // }
        // else
        // {
            // if(ndim+1 == ntags)
            // {
                // // axistags have have one additional element => drop the channel tag
                // // FIXME: would it be cleaner to make this an error ?
                // static python_ptr func(PyString_FromString("dropChannelAxis"), 
                                       // python_ptr::keep_count);
                // python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                               // python_ptr::keep_count);
                // pythonToCppException(res);
            // }
            // else
                // vigra_precondition(ndim == ntags,
                     // "constructArray(): size mismatch between shape and axistags.");
            
        // }
    // }
    // else
    // {
        // // shape has a channel axis
        // if(channelIndex == ntags)
        // {
            // // axistags have no channel axis => should be one element shorter
            // vigra_precondition(ndim == ntags+1,
                 // "constructArray(): size mismatch between shape and axistags.");
                 
            // if(shape[0] == 1)
            // {
                // // we have a singleband image => drop the channel axis
                // shape.erase(shape.begin());
                // ndim -= 1;
            // }
            // else
            // {
                // // we have a multiband image => add a channel tag
                // static python_ptr func(PyString_FromString("insertChannelAxis"), 
                                       // python_ptr::keep_count);
                // python_ptr res(PyObject_CallMethodObjArgs(axistags, func.get(), NULL), 
                               // python_ptr::keep_count);
                // pythonToCppException(res);
            // }
        // }
        // else
        // {
            // // axistags have channel axis => sizes should match
            // vigra_precondition(ndim == ntags,
                 // "constructArray(): size mismatch between shape and axistags.");
        // }
    // }
    
    // return shape;
// }

inline 
ArrayVector<npy_intp> unifyTaggedShapeSize(TaggedShape & tagged_shape)
{
    PyAxisTags axistags = tagged_shape.axistags;
    ArrayVector<npy_intp> shape = tagged_shape.shape;

    int ndim = (int)shape.size();
    int ntags = axistags.size();
    
    long channelIndex = axistags.channelIndex();

    if(tagged_shape.channelAxis == TaggedShape::none)
    {
        // shape has no channel axis
        if(channelIndex == ntags)
        {
            // axistags have no channel axis either => sizes should match
            vigra_precondition(ndim == ntags,
                 "constructArray(): size mismatch between shape and axistags.");
        }
        else
        {
            if(ndim+1 == ntags)
                // axistags have have one additional element => drop the channel tag
                // FIXME: would it be cleaner to make this an error ?
                axistags.dropChannelAxis();
            else
                vigra_precondition(ndim == ntags,
                     "constructArray(): size mismatch between shape and axistags.");
            
        }
    }
    else
    {
        // shape has a channel axis
        if(channelIndex == ntags)
        {
            // axistags have no channel axis => should be one element shorter
            vigra_precondition(ndim == ntags+1,
                 "constructArray(): size mismatch between shape and axistags.");
                 
            if(shape[0] == 1)
            {
                // we have a singleband image => drop the channel axis
                shape.erase(shape.begin());
                ndim -= 1;
            }
            else
            {
                // we have a multiband image => add a channel tag
                axistags.insertChannelAxis();
            }
        }
        else
        {
            // axistags have channel axis => sizes should match
            vigra_precondition(ndim == ntags,
                 "constructArray(): size mismatch between shape and axistags.");
        }
    }
    
    return shape;
}

inline // FIXME
ArrayVector<npy_intp> finalizeTaggedShape(TaggedShape & tagged_shape)
{
    tagged_shape.rotateToNormalOrder();
    
    if(tagged_shape.axistags)
    {
        // we assume here that the axistag object belongs to the array to be created
        // so that we can freely edit it
        scaleAxisResolution(tagged_shape);
                
        if(tagged_shape.channelDescription != "")
            tagged_shape.axistags.setChannelDescription(tagged_shape.channelDescription);
            
        // this must be last, as it may destroy snyc between shape and original_shape
        return unifyTaggedShapeSize(tagged_shape);
    }
    else
    {
        return tagged_shape.shape;
    }
}

} // namespace vigra

#endif // VIGRA_NUMPY_ARRAY_TAGGEDSHAPE_HXX
