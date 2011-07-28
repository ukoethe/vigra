/************************************************************************/
/*                                                                      */
/*               Copyright 2010-2011 by Ullrich Koethe                  */
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

#ifndef VIGRA_AXISTAGS_HXX
#define VIGRA_AXISTAGS_HXX

#include "utilities.hxx"
#include "array_vector.hxx"
#include "algorithm.hxx"
#include "error.hxx"
#include "functorexpression.hxx"
#include <string>
#include <sstream>
#include <iomanip>

namespace vigra {

class AxisInfo
{
  public:
  
    // this particular assignment of bits to types is crucial for 
    // canonical axis ordering
    enum AxisType { Channels = 1, 
                    Space = 2, 
                    Angle = 4, 
                    Time = 8, 
                    Frequency = 16, 
                    UnknownAxisType = 32, 
                    NonChannel = Space | Angle | Time | Frequency | UnknownAxisType,
                    AllAxes = 2*UnknownAxisType-1 };

    AxisInfo(std::string key = "?", AxisType typeFlags = UnknownAxisType, 
             double resolution = 0.0, std::string description = "")
    : key_(key),
      description_(description),
      resolution_(resolution),
      flags_(typeFlags)
    {}
    
    std::string key() const
    {
        return key_;
    }
    
    std::string description() const
    {
        return description_;
    }
    
    void setDescription(std::string const & description)
    {
        description_ = description;
    }
    
    double resolution() const
    {
        return resolution_;
    }
    
    void setResolution(double resolution)
    {
        resolution_ = resolution;
    }
    
    AxisType typeFlags() const
    {
        return flags_ == 0
                  ? UnknownAxisType
                  : flags_;
    }
    
    bool isUnknown() const
    {
        return isType(UnknownAxisType);
    }
    
    bool isSpatial() const
    {
        return isType(Space);
    }
    
    bool isTemporal() const
    {
        return isType(Time);
    }
    
    bool isChannel() const
    {
        return isType(Channels);
    }
    
    bool isFrequency() const
    {
        return isType(Frequency);
    }
    
    bool isAngular() const
    {
        return isType(Angle);
    }
    
    bool isType(AxisType type) const
    {
        return (typeFlags() & type) != 0;
    }
    
    std::string repr() const
    {
        std::string res("AxisInfo: '");
        res += key_ + "' (type:";
        if(isUnknown())
        {
            res += " none";
        }
        else
        {
            if(isChannel())
                res += " Channels";
            if(isSpatial())
                res += " Space";
            if(isTemporal())
                res += " Time";
            if(isAngular())
                res += " Angle";
            if(isFrequency())
                res += " Frequency";
        }
        if(resolution_ > 0.0)
        {
            res += ", resolution=";
            res += asString(resolution_);
        }
        res += ")";
        if(description_ != "")
        {
            res += " ";
            res += description_;
        }
        return res;
    }
    
    AxisInfo toFrequencyDomain(unsigned int size = 0, int sign = 1) const
    {
        AxisType type;
        if(sign == 1)
        {
            vigra_precondition(!isFrequency(),
                "AxisInfo::toFrequencyDomain(): axis is already in the Fourier domain.");
            type = AxisType(Frequency | flags_);
        }
        else
        {
            vigra_precondition(isFrequency(),
                "AxisInfo::fromFrequencyDomain(): axis is not in the Fourier domain.");
            type = AxisType(~Frequency & flags_);
        }
        AxisInfo res(key(), type, 0.0, description_);
        if(resolution_ > 0.0 && size > 0u)
            res.resolution_ = 1.0 / (resolution_ * size);
        return res;
    }
    
    AxisInfo fromFrequencyDomain(unsigned int size = 0) const
    {
        return toFrequencyDomain(size, -1);
    }

    bool compatible(AxisInfo const & other) const
    {
        return isUnknown() || other.isUnknown() || 
               ((typeFlags() & ~Frequency) == (other.typeFlags() & ~Frequency) &&
                 key() == other.key());
    }

    bool operator==(AxisInfo const & other) const
    {
        return typeFlags() == other.typeFlags() && key() == other.key();
    }

    bool operator!=(AxisInfo const & other) const
    {
        return !operator==(other);
    }
    
    // the primary ordering is according to axis type:
    //     Channels < Space < Angle < Time < Frequency < Unknown
    // the secondary ordering is the lexicographic ordering of the keys
    //     "x" < "y" < "z"
    bool operator<(AxisInfo const & other) const
    {
        return (typeFlags() < other.typeFlags()) || 
                (typeFlags() == other.typeFlags() && key() < other.key());
    }

    bool operator<=(AxisInfo const & other) const
    {
        return !(other < *this);
    }
    
    bool operator>(AxisInfo const & other) const
    {
        return other < *this;
    }
    
    bool operator>=(AxisInfo const & other) const
    {
        return !(*this < other);
    }
    
    // factory functions for standard tags
    static AxisInfo x(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("x", Space, resolution, description);
    }
    
    static AxisInfo y(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("y", Space, resolution, description);
    }
    
    static AxisInfo z(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("z", Space, resolution, description);
    }
    
    static AxisInfo t(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("t", Time, resolution, description);
    }
    
    static AxisInfo fx(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("x", AxisType(Space | Frequency), resolution, description);
    }
    
    static AxisInfo fy(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("y", AxisType(Space | Frequency), resolution, description);
    }
    
    static AxisInfo fz(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("z", AxisType(Space | Frequency), resolution, description);
    }
    
    static AxisInfo ft(double resolution = 0.0, std::string const & description = "")
    {
        return AxisInfo("t", AxisType(Time | Frequency), resolution, description);
    }
    
    static AxisInfo c(std::string const & description = "")
    {
        return AxisInfo("c", Channels, 0.0, description);
    }
    
    std::string key_, description_;
    double resolution_;
    AxisType flags_;
};

class AxisTags
{
  public:
   
    AxisTags()
    {}
    
    AxisTags(int size)
    : axes_(size)
    {}
    
    AxisTags(AxisInfo const & i1)
    {
        push_back(i1);
    }
    
    AxisTags(AxisInfo const & i1, AxisInfo const & i2)
    {
        push_back(i1);
        push_back(i2);
    }
    
    AxisTags(AxisInfo const & i1, AxisInfo const & i2, AxisInfo const & i3)
    {
        push_back(i1);
        push_back(i2);
        push_back(i3);
    }
    
    AxisTags(AxisInfo const & i1, AxisInfo const & i2,
             AxisInfo const & i3, AxisInfo const & i4)
    {
        push_back(i1);
        push_back(i2);
        push_back(i3);
        push_back(i4);
    }
    
    AxisTags(AxisInfo const & i1, AxisInfo const & i2,
             AxisInfo const & i3, AxisInfo const & i4, AxisInfo const & i5)
    {
        push_back(i1);
        push_back(i2);
        push_back(i3);
        push_back(i4);
        push_back(i5);
    }
    
    // static AxisTags fromJSON(std::string const & repr);

    std::string toJSON() const
    {
        std::stringstream s;
        s << "{\n  \"axes\": [";
        for(unsigned int k=0; k<size(); ++k)
        {
            if(k > 0)
                s << ",";
            s << "\n";
            s << "    {\n";
            s << "      \"key\": \"" << axes_[k].key() << "\",\n";
            s << "      \"typeFlags\": " << (unsigned int)axes_[k].typeFlags() << ",\n";
            s << "      \"resolution\": " << std::setprecision(17) << axes_[k].resolution() << ",\n";
            s << "      \"description\": \"" << axes_[k].description() << "\"\n";
            s << "    }";
        }
        s << "\n  ]\n}";
        return s.str();
    }

    unsigned int size() const
    {
        return axes_.size();
    }
    
    int axisTypeCount(AxisInfo::AxisType type) const
    {
        int res = 0;
        for(unsigned int k=0; k<size(); ++k)
            if(axes_[k].isType(type))
                ++res;
        return res;
    }
    
    std::string repr() const
    {
        std::string res;
        if(size() > 0)
            res += axes_[0].key();
        for(unsigned int k=1; k<size(); ++k)
        {
            res += " ";
            res += axes_[k].key();
        }
        return res;
    }
    
    AxisInfo & get(int k)
    {
        checkIndex(k);
        if(k < 0)
            k += size();
        return axes_[k];
    }
    
    AxisInfo & get(std::string const & key)
    {
        return get(index(key));
    }
    
    AxisInfo const & get(int k) const
    {
        checkIndex(k);
        if(k < 0)
            k += size();
        return axes_[k];
    }
    
    AxisInfo const & get(std::string const & key) const
    {
        return get(index(key));
    }
    
    void set(int k, AxisInfo const & info)
    {
        checkIndex(k);
        if(k < 0)
            k += size();
        checkDuplicates(k, info);
        axes_[k] = info;
    }

    void set(std::string const & key, AxisInfo const & info)
    {
        set(index(key), info);
    }

    void insert(int k, AxisInfo const & i)
    {
        if(k == (int)size())
        {
            push_back(i);
        }
        else
        {
            checkIndex(k);
            if(k < 0)
                k += size();
            checkDuplicates(size(), i);
            axes_.insert(axes_.begin()+k, i);
        }
    }

    void push_back(AxisInfo const & i)
    {
        checkDuplicates(size(), i);
        axes_.push_back(i);
    }
    
    void dropAxis(int k)
    {
        checkIndex(k);
        ArrayVector<AxisInfo>::iterator i = k < 0
                                                 ? axes_.end() + k
                                                 : axes_.begin() + k;
        axes_.erase(i, i+1);
    }
    
    void dropAxis(std::string const & key)
    {
        dropAxis(index(key));
    }
    
    void dropChannelAxis()
    {
        int k = channelIndex();
        if(k < (int)size())
            axes_.erase(axes_.begin() + k, axes_.begin() + k + 1);
    }
    
    int index(std::string const & key) const
    {
        for(unsigned int k=0; k<size(); ++k)
            if(axes_[k].key() == key)
                return k;
        return (int)size();
    }
    
    double resolution(int k) const
    {
        return get(k).resolution_;
    }
    
    double resolution(std::string const & key) const
    {
        return resolution(index(key));
    }
    
    void setResolution(int k, double r) 
    {
        get(k).resolution_ = r;
    }
    
    void setResolution(std::string const & key, double r) 
    {
        setResolution(index(key), r);
    }
    
    void scaleResolution(int k, double factor)
    {
        get(k).resolution_ *= factor;
    }
    
    void scaleResolution(std::string const & key, double factor)
    {
        get(key).resolution_ *= factor;
    }
    
    std::string description(int k) const
    {
        return get(k).description_;
    }
    
    std::string description(std::string const & key) const
    {
        return description(index(key));
    }
    
    void setDescription(int k, std::string const & d) 
    {
        get(k).setDescription(d);
    }
    
    void setDescription(std::string const & key, std::string const & d) 
    {
        setDescription(index(key), d);
    }
    
    void setChannelDescription(std::string const & description)
    {
        int k = channelIndex();
        if(k < (int)size())
            axes_[k].setDescription(description);
    }
    
    void toFrequencyDomain(int k, int size = 0, int sign = 1)
    {
        get(k) = get(k).toFrequencyDomain(size, sign);
    }
    
    void toFrequencyDomain(std::string const & key, int size = 0, int sign = 1)
    {
        toFrequencyDomain(index(key), size, sign);
    }
    
    void fromFrequencyDomain(int k, int size = 0)
    {
        toFrequencyDomain(k, size, -1);
    }
    
    void fromFrequencyDomain(std::string const & key, int size = 0)
    {
        toFrequencyDomain(key, size, -1);
    }
    
    bool hasChannelAxis() const
    {
        return channelIndex() != (int)size();
    }
    
    // FIXME: cache the results of these functions?
    int channelIndex() const
    {
        for(unsigned int k=0; k<size(); ++k)
            if(axes_[k].isChannel())
                return k;
        return (int)size();
    }
    
    int majorNonchannelIndex() const
    {
        int k = 0;
        for(; k<(int)size(); ++k)
            if(!axes_[k].isChannel())
                break;
        for(int i=k+1; i<(int)size(); ++i)
        {
            if(axes_[i].isChannel())
                continue;
            if(axes_[i] < axes_[k])
                k = i;
        }
        return k;
    }
    
    void swapaxes(int i1, int i2)
    {
        checkIndex(i1);
        checkIndex(i2);
        if(i1 < 0)
            i1 += size();
        if(i2 < 0)
            i2 += size();
        std::swap(axes_[i1], axes_[i2]);
    }
    
    template <class T>
    void transpose(ArrayVector<T> const & permutation)
    {
        if(permutation.size() == 0)
        {
            transpose();
        }
        else
        {
            vigra_precondition(permutation.size() == size(),
                "AxisTags::transpose(): Permutation has wrong size.");
            ArrayVector<AxisInfo> newAxes(size());
            applyPermutation(permutation.begin(), permutation.end(), axes_.begin(), newAxes.begin());
            axes_.swap(newAxes);
        }
    }
    
    void transpose()
    {
        std::reverse(axes_.begin(), axes_.end());
    }
    
    template <class T>
    void permutationToNormalOrder(ArrayVector<T> & permutation)
    {
        permutation.resize(size());
        indexSort(axes_.begin(), axes_.end(), permutation.begin());
    }
    
    template <class T>
    void permutationToNormalOrder(ArrayVector<T> & permutation, AxisInfo::AxisType types)
    {
        ArrayVector<AxisInfo> matchingAxes;
        for(int k=0; k<(int)size(); ++k)
            if(axes_[k].isType(types))
                matchingAxes.push_back(axes_[k]);
        permutation.resize(matchingAxes.size());
        indexSort(matchingAxes.begin(), matchingAxes.end(), permutation.begin());
    }
    
    template <class T>
    void permutationFromNormalOrder(ArrayVector<T> & inverse_permutation)
    {
        ArrayVector<T> permutation;
        permutationToNormalOrder(permutation);
        inverse_permutation.resize(permutation.size());
        indexSort(permutation.begin(), permutation.end(), inverse_permutation.begin());
    }   
    
    template <class T>
    void permutationFromNormalOrder(ArrayVector<T> & inverse_permutation, AxisInfo::AxisType types)
    {
        ArrayVector<T> permutation;
        permutationToNormalOrder(permutation, types);
        inverse_permutation.resize(permutation.size());
        indexSort(permutation.begin(), permutation.end(), inverse_permutation.begin());
    }   
    
    template <class T>
    void permutationToNumpyOrder(ArrayVector<T> & permutation)
    {
        permutation.resize(size());
        indexSort(axes_.begin(), axes_.end(), permutation.begin(), std::greater<AxisInfo>());
    }
    
    template <class T>
    void permutationFromNumpyOrder(ArrayVector<T> & inverse_permutation)
    {
        ArrayVector<T> permutation;
        permutationToNumpyOrder(permutation);
        inverse_permutation.resize(permutation.size());
        indexSort(permutation.begin(), permutation.end(), inverse_permutation.begin());
    }   
    
#if 0
    ArrayVector<UInt32> matchOrdering(AxisTags const & other)
    {
        vigra_precondition(size() == other.size(),
            "AxisTags::matchOrdering(): size mismatch.");
        
        ArrayVector<UInt32> permutation(size());
        for(unsigned int k = 0; k<size(); ++k)
        {
            std::string key = other.get(k).key();
            unsigned int l=0;
            for(; l<size(); ++l)
            {
                if(key == get(l).key())
                    break;
            }
            vigra_precondition(l < size(),
                "AxisTags::matchOrdering(): key mismatch.");
            permutation[k] = l;
        }
        return permutation;
    }    
#endif

    bool compatible(AxisTags const & other) const
    {
        if(size() == 0 || other.size() == 0)
            return true;
        if(size() != other.size())
            return false;
        for(unsigned int k=0; k<size(); ++k)
            if(!axes_[k].compatible(other.axes_[k]))
                return false;
        return true;
    }

    bool operator==(AxisTags const & other) const
    {
        if(size() != other.size())
            return false;
        return std::equal(axes_.begin(), axes_.end(), other.axes_.begin());
    }

    bool operator!=(AxisTags const & other) const
    {
        return !operator==(other);
    }
    
  protected:
    
    void checkIndex(int k) const
    {
        vigra_precondition(k < (int)size() && k >= -(int)size(),
            "AxisTags::checkIndex(): index out of range.");
    }

    void checkDuplicates(int i, AxisInfo const & info)
    {
        if(info.isChannel())
        {  
            for(int k=0; k<(int)size(); ++k)
            {
                vigra_precondition(k == i || !axes_[k].isChannel(),
                     "AxisTags::checkDuplicates(): can only have one channel axis.");
            }
        }
        else if(!info.isUnknown())
        {
            for(int k=0; k<(int)size(); ++k)
            {
                vigra_precondition(k == i || axes_[k].key() != info.key(),
                     std::string("AxisTags::checkDuplicates(): axis key '" + 
                                  info.key() + "' already exists."));
            }
        }
    }

    ArrayVector<AxisInfo> axes_;
};

// #if 0
// struct PyGetFunctor
// {
    // AxisInfo const & operator()(python::object const & o) const
    // {
        // return python::extract<AxisInfo const &>(o)();
    // }
// };

// class PyAxisTags
// : public AxisTags<python::object, PyGetFunctor>
// {
    // typedef AxisTags<python::object, PyGetFunctor> BaseType;
  // public:
    // PyAxisTags()
    // {}

    // PyAxisTags(python::object i1, python::object i2,
             // python::object i3, python::object i4, python::object i5)
    // {
        // if(PySequence_Check(i1.ptr()))
        // {
            // int size = len(i1);
            // for(int k=0; k<size; ++k)
                // if(python::extract<AxisInfo const &>(i1[k]).check())
                    // push_back(i1[k]);
        // }
        // else if(PyInt_Check(i1.ptr()))
        // {
            // int size = python::extract<int>(i1)();
            // for(int k=0; k<size; ++k)
                // push_back(python::object(AxisInfo()));
        // }
        // else
        // {
            // if(python::extract<AxisInfo const &>(i1).check())
                // push_back(i1);
            // if(python::extract<AxisInfo const &>(i2).check())
                // push_back(i2);
            // if(python::extract<AxisInfo const &>(i3).check())
                // push_back(i3);
            // if(python::extract<AxisInfo const &>(i4).check())
                // push_back(i4);
            // if(python::extract<AxisInfo const &>(i5).check())
                // push_back(i5);
        // }
    // }
    
    // python::object getitem(int k)
    // {
        // if(!checkIndex(k))
        // {
            // PyErr_SetString(PyExc_IndexError, "AxisInfo::getitem(): Invalid index or key.");
            // python::throw_error_already_set();
        // }
        // if(k < 0)
            // k += this->size();
        // return this->axes_[k];
    // }
    
    // python::object getitem(std::string const & key)
    // {
        // return getitem(this->findKey(key));
    // }

    // void setitem(int k, python::object i)
    // {
        // if(!this->checkIndex(k))
        // {
            // PyErr_SetString(PyExc_IndexError, "AxisInfo::setitem(): Invalid index or key.");
            // python::throw_error_already_set();
        // }
        // if(!python::extract<AxisInfo const &>(i).check())
        // {
            // PyErr_SetString(PyExc_TypeError, "AxisInfo::setitem(): Item type must be AxisInfo.");
            // python::throw_error_already_set();
        // }
        
        // if(k < 0)
            // k += this->size();
        // this->axes_[k] = i;
    // }
    
    // void setitem(std::string const & key, python::object i)
    // {
        // setitem(this->findKey(key), i);
    // }

    // void append(python::object i)
    // {
        // insert(size(), i);
    // }

    // void insert(int k, python::object i)
    // {
        // if(k < 0)
            // k += this->size();
        // if(k < 0)
            // k = 0;
        // if(k > (int)this->size())
            // k = this->size();
        // if(!python::extract<AxisInfo const &>(i).check())
        // {
            // PyErr_SetString(PyExc_TypeError, "AxisInfo::insert(): Item type must be AxisInfo.");
            // python::throw_error_already_set();
        // }
        // this->axes_.insert(this->axes_.begin()+k, i);
    // }
    
    // void insert(std::string const & key, python::object i)
    // {
        // insert(this->findKey(key), i);
    // }

    // python::list axesByFlag(AxisType typeFlags) const
    // {
        // python::list res;
        // for(unsigned int k=0; k<this->size(); ++k)
            // if(this->get(k).typeFlags() == typeFlags)
                // res.append(k);
        // return res;
    // }

    // python::list spatialAxes() const
    // {
        // python::list res;
        // for(unsigned int k=0; k<this->size(); ++k)
            // if(this->get(k).isSpatial())
                // res.append(k);
        // return res;
    // }

    // python::list temporalAxes() const
    // {
        // python::list res;
        // for(unsigned int k=0; k<this->size(); ++k)
            // if(this->get(k).isTemporal())
                // res.append(k);
        // return res;
    // }

    // python::list channelAxes() const
    // {
        // python::list res;
        // for(unsigned int k=0; k<this->size(); ++k)
            // if(this->get(k).isChannel())
                // res.append(k);
        // return res;
    // }

    // python::list frequencyAxes() const
    // {
        // python::list res;
        // for(unsigned int k=0; k<this->size(); ++k)
            // if(this->get(k).isFrequency())
                // res.append(k);
        // return res;
    // }

    // python::list angularAxes() const
    // {
        // python::list res;
        // for(unsigned int k=0; k<this->size(); ++k)
            // if(this->get(k).isAngular())
                // res.append(k);
        // return res;
    // }

    // python::list untaggedAxes() const
    // {
        // python::list res;
        // for(unsigned int k=0; k<this->size(); ++k)
            // if(this->get(k).isUnknown())
                // res.append(k);
        // return res;
    // }

    // template <class U>
    // python::list vectorToPython(ArrayVector<U> const & v) const
    // {
        // python::list res;
        // for(unsigned int k=0; k<v.size(); ++k)
            // res.append(v[k]);
        // return res;
    // }

    // python::list canonicalOrdering()
    // {
        // return vectorToPython(BaseType::canonicalOrdering());
    // }

    // python::list matchOrdering(PyAxisTags const & other)
    // {
        // return vectorToPython(BaseType::matchOrdering(other));
    // }

    // void transpose(python::object const & o)
    // {
        // unsigned int osize = len(o);
        // ArrayVector<UInt32> permutation(osize);
        
        // for(unsigned int k=0; k<this->size(); ++k)
            // permutation[k] = python::extract<UInt32>(o[k])();
            
        // BaseType::transpose(permutation);
    // }

    // void transpose()
    // {
        // BaseType::transpose();
    // }
// };

// class TaggedShape
// {
  // public:

    // ArrayVector<npy_intp> shape;
    // python_ptr axistags;
    // npy_intp channelCount;
    // std::string channelDescription;
    
    // TaggedShape(MultiArrayIndex size)
    // : shape(size)
    // {}
    
    // template <int N>
    // TaggedShape(typename MultiArrayShape<N>::type const & sh)
    // : shape(sh.begin(), sh.end())
    // {}
    
    // npy_intp & operator[](int i)
    // {
        // // rotate indices so that channels are located at index 0
        // return shape[(i+1) % shape.size()];
    // }
    
    // npy_intp operator[](int i) const
    // {
        // return shape[(i+1) % shape.size()];
    // }
    
    // unsigned int size() const
    // {
        // return shape.size();
    // }
    
    // // void setChannelDescription(std::string const & description)
    // // {
        // // if(axistags)
        // // {
            // // python_ptr func(PyString_FromString("setChannelDescription"), 
                                         // // python_ptr::keep_count);
            // // pythonToCppException(res);
            
            // // python_ptr d(PyString_FromString(d.c_str()), python_ptr::keep_count);
            // // pythonToCppException(d);
            
            // // python_ptr res(PyObject_CallMethodObjArgs(axistags, func, d.get(), NULL),
                           // // python_ptr::keep_count);
            // // pythonToCppException(res);
        // // }
    // // }
    
    // // void setChannelCount(int channelCount)
    // // {
        // // shape[0] = channelCount;
    // // }
    
    // void setChannelDescription(std::string const & description)
    // {
        // channelDescription = description;
    // }
    
    // void setChannelCount(int count)
    // {
        // channelCount = count;
    // }
    
    // void setChannelConfig(int channelCount, std::string const & description)
    // {
        // setChannelCount(channelCount);
        // setChannelDescription(description);
    // }
// };
// #endif

} // namespace vigra

#endif /* VIGRA_AXISTAGS_HXX */
