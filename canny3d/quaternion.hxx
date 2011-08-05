/************************************************************************/
/*                                                                      */
/*                  Copyright 2004-2005 by Hans Meine                   */
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

#ifndef VIGRA_QUATERNION_HXX
#define VIGRA_QUATERNION_HXX

#include <cmath>
#include <vigra/config.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/matrix.hxx>

namespace vigra {

template<class ValueType>
class Quaternion {
  public:
    typedef ValueType value_type;
    typedef TinyVector<ValueType, 3> Vector;

    Quaternion()
    {}
    Quaternion(const Quaternion &q)
    : w_(q.w_), v_(q.v_)
    {}
    Quaternion(ValueType w, const Vector &v)
    : w_(w), v_(v)
    {}
    Quaternion(ValueType w, ValueType x, ValueType y, ValueType z)
    : w_(w), v_(x, y, z)
    {}

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

    typename NumericTraits<ValueType>::RealPromote
    magnitude() const
    {
        return VIGRA_CSTD::sqrt(
            (typename NumericTraits<ValueType>::RealPromote)squaredMagnitude());
    }

    typename NumericTraits<ValueType>::Promote
    squaredMagnitude() const
    {
        return w_*w_ + v_.squaredMagnitude();
    }

    Quaternion &operator+=(Quaternion const &other)
    {
        w_ += other.w_;
        v_ += other.v_;
        return *this;
    }
    Quaternion operator+(Quaternion const &other) const
    {
        Quaternion result(*this);
        result += other;
        return result;
    }

    Quaternion &operator-=(Quaternion const &other)
    {
        w_ -= other.w_;
        v_ -= other.v_;
        return *this;
    }
    Quaternion operator-(Quaternion const &other) const
    {
        Quaternion result(*this);
        result -= other;
        return result;
    }

    Quaternion &operator*=(Quaternion const &other)
    {
        value_type newW(w_*other.w_ - dot(v_, other.v_));
        // if TinyVector had ^ for cross-product:
        // v_ = w_ * other.v_ + other.w_ * v + v_ ^ other.v_;
        v_ = Vector(w_*other.x() + x()*other.w_  + y()*other.z() - z()*other.y(),
                    w_*other.y() - x()*other.z() + y()*other.w_  + z()*other.x(),
                    w_*other.z() + x()*other.y() - y()*other.x() + z()*other.w_ );
        w_ = newW;
        return *this;
    }
    Quaternion operator*(Quaternion const &other) const
    {
        Quaternion result(*this);
        result *= other;
        return result;
    }

    Quaternion &operator*=(double scale)
    {
        w_ *= scale;
        v_ *= scale;
        return *this;
    }
    Quaternion operator*(double scale) const
    {
        Quaternion result(*this);
        result *= scale;
        return result;
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



} // namespace vigra

#include <iostream>

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
