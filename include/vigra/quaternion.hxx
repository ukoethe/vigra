/************************************************************************/
/*                                                                      */
/*        Copyright 2004-2010 by Hans Meine und Ullrich Koethe          */
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

#ifndef VIGRA_QUATERNION_HXX
#define VIGRA_QUATERNION_HXX

#include "config.hxx"
#include "numerictraits.hxx"
#include "tinyvector.hxx"
#include "matrix.hxx"
#include <iosfwd>   // ostream


namespace vigra {

template<class ValueType>
class Quaternion {
  public:
    typedef TinyVector<ValueType, 3> Vector;
    
        /** the quaternion's valuetype
        */
    typedef ValueType value_type;

        /** reference (return of operator[]).
        */
    typedef ValueType & reference;

        /** const reference (return of operator[] const).
        */
    typedef ValueType const & const_reference;

        /** the quaternion's squared norm type
        */
    typedef typename NormTraits<ValueType>::SquaredNormType SquaredNormType;

        /** the quaternion's norm type
        */
    typedef typename SquareRootTraits<SquaredNormType>::SquareRootResult NormType;


    Quaternion(ValueType w = 0, ValueType x = 0, ValueType y = 0, ValueType z = 0)
    : w_(w), v_(x, y, z)
    {}
    
    Quaternion(ValueType w, const Vector &v)
    : w_(w), v_(v)
    {}

    Quaternion(const Quaternion &q)
    : w_(q.w_), v_(q.v_)
    {}
    
    Quaternion & operator=(Quaternion const & other)
    {
        w_ = other.w_;
        v_ = other.v_;
        return *this;
    }
    
    Quaternion & operator=(ValueType other)
    {
        w_ = other;
        v_.init(0);
        return *this;
    }

        /**
         * Creates a Quaternion which represents the operation of
         * rotating around the given axis by the given angle.
         *
         * The angle should be in the range -pi..3*pi for sensible
         * results.
         */
    static Quaternion
    createRotation(double angle, const Vector &rotationAxis)
    {
        // the natural range would be -pi..pi, but the reflective
        // behavior around pi is too unexpected:
        if(angle > M_PI)
            angle -= 2.0*M_PI;
        double t(VIGRA_CSTD::sin(angle/2));
		double norm(rotationAxis.magnitude());
        return Quaternion(VIGRA_CSTD::sqrt(1-t*t), t*rotationAxis/norm);
    }

        // rename to s
    ValueType w() const { return w_; }
    ValueType &w() { return w_; }
    void setW(ValueType w) { w_ = w; }

    const Vector &v() const { return v_; }
    Vector &v() { return v_; }
    void setV(const Vector & v) { v_ = v; }
    void setV(ValueType x, ValueType y, ValueType z)
    {
        v_[0] = x;
        v_[1] = y;
        v_[2] = z;
    }

    ValueType x() const { return v_[0]; }
    ValueType y() const { return v_[1]; }
    ValueType z() const { return v_[2]; }
    ValueType &x() { return v_[0]; }
    ValueType &y() { return v_[1]; }
    ValueType &z() { return v_[2]; }
    void setX(ValueType x) { v_[0] = x; }
    void setY(ValueType y) { v_[1] = y; }
    void setZ(ValueType z) { v_[2] = z; }
    
    value_type & operator[](int index)
    {
        return index == 0
                   ? w_
                   : v_[index - 1];
    }
    
    value_type operator[](int index) const
    {
        return index == 0
                   ? w_
                   : v_[index - 1];
    }
    
    NormType magnitude() const
    {
        return VIGRA_CSTD::sqrt((NormType)squaredMagnitude());
    }

    SquaredNormType squaredMagnitude() const
    {
        return w_*w_ + v_.squaredMagnitude();
    }

    Quaternion &operator+=(value_type const &w)
    {
        w_ += w;
        return *this;
    }

    Quaternion &operator+=(Quaternion const &other)
    {
        w_ += other.w_;
        v_ += other.v_;
        return *this;
    }

    Quaternion &operator-=(value_type const &w)
    {
        w_ -= w;
        return *this;
    }

    Quaternion &operator-=(Quaternion const &other)
    {
        w_ -= other.w_;
        v_ -= other.v_;
        return *this;
    }

    Quaternion operator-() const
    {
        return Quaternion(-w_, -v_);
    }

    Quaternion &operator*=(Quaternion const &other)
    {
        value_type newW = w_*other.w_ - dot(v_, other.v_);
        v_              = w_ * other.v_ + other.w_ * v + cross(v_, other.v_);
        w_              = newW;
        return *this;
    }

    Quaternion &operator*=(double scale)
    {
        w_ *= scale;
        v_ *= scale;
        return *this;
    }

    Quaternion &operator/=(Quaternion const &other)
    {
        (*this) *= conj(other) / squaredNorm(other);
        return *this;
    }

    Quaternion &operator/=(double scale)
    {
        w_ /= scale;
        v_ /= scale;
        return *this;
    }
    
    Quaternion operator/(double scale) const
    {
        Quaternion result(*this);
        result /= scale;
        return result;
    }

    bool operator==(Quaternion const &other) const
    {
      return (w_ == other.w_) && (v_ == other.v_);
    }

    bool operator!=(Quaternion const &other) const
    {
      return (w_ != other.w_) || (v_ != other.v_);
    }

        /**
         * Fill the first 3x3 elements of the given matrix with a
         * rotation matrix performing the same 3D rotation as this
         * quaternion.  If matrix is in column-major format, it should
         * be pre-multiplied with the vectors to be rotated, i.e.
         * matrix[0][0-3] will be the rotated X axis.
         */
    template<class MatrixType>
    void fillRotationMatrix(MatrixType &matrix) const
    {
        // scale by 2 / norm
        typename NumericTraits<ValueType>::RealPromote s =
            2 / (typename NumericTraits<ValueType>::RealPromote)squaredMagnitude();

        Vector
            vs = v_ * s,
            wv = w_ * vs,
            vv = vs * v_;
        value_type
            xy = vs[0] * v_[1],
            xz = vs[0] * v_[2],
            yz = vs[1] * v_[2];

        matrix[0][0] = 1 - (vv[1] + vv[2]);
        matrix[0][1] =     ( xy   - wv[2]);
        matrix[0][2] =     ( xz   + wv[1]);

        matrix[1][0] =     ( xy   + wv[2]);
        matrix[1][1] = 1 - (vv[0] + vv[2]);
        matrix[1][2] =     ( yz   - wv[0]);

        matrix[2][0] =     ( xz   - wv[1]);
        matrix[2][1] =     ( yz   + wv[0]);
        matrix[2][2] = 1 - (vv[0] + vv[1]);
    }

    void fillRotationMatrix(Matrix<value_type> &matrix) const
    {
        // scale by 2 / norm
        typename NumericTraits<ValueType>::RealPromote s =
            2 / (typename NumericTraits<ValueType>::RealPromote)squaredMagnitude();

        Vector
            vs = v_ * s,
            wv = w_ * vs,
            vv = vs * v_;
        value_type
            xy = vs[0] * v_[1],
            xz = vs[0] * v_[2],
            yz = vs[1] * v_[2];

        matrix(0, 0) = 1 - (vv[1] + vv[2]);
        matrix(0, 1) =     ( xy   - wv[2]);
        matrix(0, 2) =     ( xz   + wv[1]);

        matrix(1, 0) =     ( xy   + wv[2]);
        matrix(1, 1) = 1 - (vv[0] + vv[2]);
        matrix(1, 2) =     ( yz   - wv[0]);

        matrix(2, 0) =     ( xz   - wv[1]);
        matrix(2, 1) =     ( yz   + wv[0]);
        matrix(2, 2) = 1 - (vv[0] + vv[1]);
    }

  protected:
    ValueType w_;
    Vector v_;
};

template<class ValueType>
Quaternion<ValueType> conj(Quaternion<ValueType> const & q)
{
    return Quaternion<ValueType>(q.w(), -q.v());
}

template<typename Type>
inline Quaternion<Type>
operator+(const Quaternion<Type>& t1,
           const Quaternion<Type>& t2) 
{
  return Quaternion<Type>(t1) += t2;
}

template<typename Type>
inline Quaternion<Type>
operator+(const Quaternion<Type>& t1,
           const Type& t2) 
{
  return Quaternion<Type>(t1) += t2;
}

template<typename Type>
inline Quaternion<Type>
operator+(const Type& t1,
           const Quaternion<Type>& t2) 
{
  return Quaternion<Type>(t1) += t2;
}

template<typename Type>
inline Quaternion<Type>
operator-(const Quaternion<Type>& t1,
           const Quaternion<Type>& t2) 
{
  return Quaternion<Type>(t1) -= t2;
}

template<typename Type>
inline Quaternion<Type>
operator-(const Quaternion<Type>& t1,
           const Type& t2) 
{
  return Quaternion<Type>(t1) -= t2;
}

template<typename Type>
inline Quaternion<Type>
operator-(const Type& t1,
           const Quaternion<Type>& t2) 
{
  return Quaternion<Type>(t1) -= t2;
}

template<typename Type>
inline Quaternion<Type>
operator*(const Quaternion<Type>& t1,
           const Quaternion<Type>& t2) 
{
  return Quaternion<Type>(t1) *= t2;
}

template<typename Type>
inline Quaternion<Type>
operator*(const Quaternion<Type>& t1,
           double t2) 
{
  return Quaternion<Type>(t1) *= t2;
}
  
template<typename Type>
inline Quaternion<Type>
operator*(double t1,
           const Quaternion<Type>& t2)
{
  return Quaternion<Type>(t1) *= t2;
}

template<typename Type>
inline Quaternion<Type>
operator/(const Quaternion<Type>& t1,
           const Quaternion<Type>& t2) 
{
  return Quaternion<Type>(t1) /= t2;
}

template<typename Type>
inline Quaternion<Type>
operator/(const Quaternion<Type>& t1,
           double t2) 
{
  return Quaternion<Type>(t1) /= t2;
}
  
template<typename Type>
inline Quaternion<Type>
operator/(double t1,
           const Quaternion<Type>& t2)
{
  return Quaternion<Type>(t1) /= t2;
}

    /// squared norm
template<typename Type>
inline
typename Quaternion<Type>::SquaredNormType
squaredNorm(Quaternion<Type> const & q)
{
    return q.squaredMagnitude();
}

} // namespace vigra

template<class ValueType>
inline
std::ostream & operator<<(std::ostream & os, vigra::Quaternion<ValueType> const & q)
{
    os << q.w() << " " << q.x() << " " << q.y() << " " << q.z();
    return os;
}

template<class ValueType>
inline
std::istream & operator>>(std::istream & is, vigra::Quaternion<ValueType> & q)
{
    ValueType w, x, y, z;
    is >> w >> x >> y >> z;
    q.setW(w);
    q.setX(x);
    q.setY(y);
    q.setZ(z);
    return is;
}

#endif // VIGRA_QUATERNION_HXX
