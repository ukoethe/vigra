/************************************************************************/
/*                                                                      */
/*               Copyright 2014-2015 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the VIGRA2 computer vision library.          */
/*    The VIGRA2 Website is                                             */
/*        http://ukoethe.github.io/vigra2                               */
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

#pragma once

#ifndef VIGRA_ANY_HXX
#define VIGRA_ANY_HXX

#include "config.hxx"
#include "error.hxx"
#include <typeinfo>
#include <type_traits>

namespace vigra {

    /** \brief Typesafe storage of arbitrary values.

        Items are always stored by value, but it is of course possible
        to store pointers and smart pointers.

        <b>Usage:</b>

        \code
        Any a(10);  // store integer '10'

        assert(a.is_type<int>());
        assert(a.is_convertible<int>());

        // retrieve the stored value (throws when types don't match)
        assert(a.get<int>() == 10);

        // change the value
        a = 20;
        assert(a.get<int>() == 20);

        // values of arithmetic types can be converted into each other
        // (this is currently not implemented for other types)
        assert(a.is_convewrtible<double>());
        assert(a.cast<double>() == 20.0);

        // delete the stored value
        a.release();
        assert(a.empty());
        assert(a == false);

        // store a shared_ptr
        typedef std::shared_ptr<int> Ptr;
        Any p(Ptr(new int(5))), q = p;
        assert(*(p.get<Ptr>()) == 5);
        // the addresses of the elements in p and q are the same
        assert(p.get<Ptr>().get() == p.get<Ptr>().get());
        \endcode
    */
class Any
{
    struct Handle
    {
        Handle() {}
        virtual ~Handle() {}
        virtual const std::type_info & type() const = 0;
        virtual Handle * clone() const = 0;
        virtual bool equal(Handle const *) const = 0;

      private:
        Handle(Handle const &);
        Handle & operator=(Handle const &);
    };

    template <class T>
    struct TypedHandle
    : public Handle
    {
        T value_;

        TypedHandle(T const & t)
        : value_(t)
        {}

        const std::type_info & type() const
        {
            return typeid(T);
        }

        Handle * clone() const
        {
            return new TypedHandle(value_);
        }

        bool equal(Handle const * h) const
        {
            TypedHandle const * ptr = dynamic_cast<TypedHandle const *>(h);
            return ptr != 0 && value_ == ptr->value_;
        }
    };

    struct ConvertibleHandle
    : public Handle
    {
        template <class T>
        struct TypeTag {};

        ConvertibleHandle() {}

        virtual signed char      cast(TypeTag<signed char>) const = 0;
        virtual signed short     cast(TypeTag<signed short>) const = 0;
        virtual signed int       cast(TypeTag<signed int>) const = 0;
        virtual signed long      cast(TypeTag<signed long>) const = 0;
        virtual signed long long cast(TypeTag<signed long long>) const = 0;

        virtual unsigned char      cast(TypeTag<unsigned char>) const = 0;
        virtual unsigned short     cast(TypeTag<unsigned short>) const = 0;
        virtual unsigned int       cast(TypeTag<unsigned int>) const = 0;
        virtual unsigned long      cast(TypeTag<unsigned long>) const = 0;
        virtual unsigned long long cast(TypeTag<unsigned long long>) const = 0;

        virtual float       cast(TypeTag<float>) const = 0;
        virtual double      cast(TypeTag<double>) const = 0;
        virtual long double cast(TypeTag<long double>) const = 0;
    };

    #define VIGRA_ANY_OF_CONVERTIBLE(TYPE) \
    template <> \
    struct TypedHandle<TYPE> \
    : public ConvertibleHandle \
    { \
        TYPE value_; \
        \
        TypedHandle(TYPE const & t) \
        : value_(t) \
        {} \
        \
        const std::type_info & type() const \
        { \
            return typeid(value_); \
        } \
        \
        Handle * clone() const \
        { \
            return new TypedHandle(value_); \
        } \
        \
        bool equal(Handle const * h) const \
        { \
            TypedHandle const * ptr = dynamic_cast<TypedHandle const *>(h); \
            return ptr != 0 && value_ == ptr->value_; \
        } \
        \
        virtual signed char      cast(TypeTag<signed char>) const \
        { return static_cast<signed char>(value_); } \
        virtual signed short     cast(TypeTag<signed short>) const \
        { return static_cast<signed short>(value_); } \
        virtual signed int       cast(TypeTag<signed int>) const \
        { return static_cast<signed int>(value_); } \
        virtual signed long      cast(TypeTag<signed long>) const \
        { return static_cast<signed long>(value_); } \
        virtual signed long long cast(TypeTag<signed long long>) const \
        { return static_cast<signed long long>(value_); } \
        \
        virtual unsigned char      cast(TypeTag<unsigned char>) const \
        { return static_cast<unsigned char>(value_); } \
        virtual unsigned short     cast(TypeTag<unsigned short>) const \
        { return static_cast<unsigned short>(value_); } \
        virtual unsigned int       cast(TypeTag<unsigned int>) const \
        { return static_cast<unsigned int>(value_); } \
        virtual unsigned long      cast(TypeTag<unsigned long>) const \
        { return static_cast<unsigned long>(value_); } \
        virtual unsigned long long cast(TypeTag<unsigned long long>) const \
        { return static_cast<unsigned long long>(value_); } \
        \
        virtual float       cast(TypeTag<float>) const \
        { return static_cast<float>(value_); } \
        virtual double      cast(TypeTag<double>) const \
        { return static_cast<double>(value_); } \
        virtual long double cast(TypeTag<long double>) const \
        { return static_cast<long double>(value_); } \
    };

    VIGRA_ANY_OF_CONVERTIBLE(signed char     )
    VIGRA_ANY_OF_CONVERTIBLE(signed short    )
    VIGRA_ANY_OF_CONVERTIBLE(signed int      )
    VIGRA_ANY_OF_CONVERTIBLE(signed long     )
    VIGRA_ANY_OF_CONVERTIBLE(signed long long)

    VIGRA_ANY_OF_CONVERTIBLE(unsigned char     )
    VIGRA_ANY_OF_CONVERTIBLE(unsigned short    )
    VIGRA_ANY_OF_CONVERTIBLE(unsigned int      )
    VIGRA_ANY_OF_CONVERTIBLE(unsigned long     )
    VIGRA_ANY_OF_CONVERTIBLE(unsigned long long)

    VIGRA_ANY_OF_CONVERTIBLE(float      )
    VIGRA_ANY_OF_CONVERTIBLE(double     )
    VIGRA_ANY_OF_CONVERTIBLE(long double)

    #undef VIGRA_ANY_OF_CONVERTIBLE

    VIGRA_UNIQUE_PTR<Handle> handle_;

  public:

        /** Construct empty 'Any' object.
        */
    Any()
    : handle_((Handle*)0)
    {}

        /** Construct 'Any' object holding the given value.
        */
    template <class T>
    Any(T const & t)
    : handle_(new TypedHandle<T>(t))
    {}

        /** Construct 'Any' object holding a copy of other's value.
        */
    Any(Any const & other)
    : handle_(other.handle_ ? other.handle_->clone() : (Handle*)0)
    {}

        /** Assign the given value to this 'Any' object
            (overwrites the old value, regardless of types).
        */
    template <class T>
    Any & operator=(T const & t)
    {
        handle_.reset(new TypedHandle<T>(t));
        return *this;
    }

        /** Assign a copy of other's value to this 'Any' object
            (overwrites the old value, regardless of types).
        */
    Any & operator=(Any const & other)
    {
        if(this != &other)
            handle_.reset(other.handle_ ? other.handle_->clone() : (Handle*)0);
        return *this;
    }

        /** Delete the contained object (make this 'Any' object empty).
        */
    void release()
    {
        handle_.release();
    }

        /** Exchange the value of this object with other's.
        */
    void swap(Any & other)
    {
        handle_.swap(other.handle_);
    }

        /** Exchange the value of objects l and r.
        */
    friend void swap(Any & l, Any & r)
    {
        l.swap(r);
    }

        /** Check if this object contains the same type and value as other.
            Also true if both 'Any' objects are empty.
        */
    bool operator==(Any const & other) const
    {
        return (handle_.get() == 0 && other.handle_.get() == 0) ||
               (handle_.get() != 0 && handle_->equal(other.handle_.get()));
    }

        /** Check if this object differs from other by type or value.
        */
    bool operator!=(Any const & other) const
    {
        return !operator==(other);
    }

    bool operator==(bool other) const
    {
        return bool(*this) == other;
    }

    bool operator!=(bool other) const
    {
        return bool(*this) != other;
    }

        /** Convert 'Any' to <tt>false</tt> if this object is empty, <tt>true</tt> otherwise.
        */
    operator bool() const
    {
        return handle_.get() != 0;
    }

        /** Check if this object is empty (holds no value).
        */
    bool empty() const
    {
        return handle_.get() == 0;
    }

        /** Check if this object holds a value of the given type.
        */
    template <class T>
    bool is_type() const
    {
        return bool(handle_) && dynamic_cast<TypedHandle<T> const *>(handle_.get()) != 0;
    }

        /** Check if this object's value is convertible to the given type.
            At present, this is only implemented for arithmetic types.
        */
    template <class T>
    bool is_convertible() const
    {
        return std::is_arithmetic<T>::value && bool(handle_) &&
               dynamic_cast<ConvertibleHandle const *>(handle_.get()) != 0;
    }

        /** Read-write access to the contained value. This throws an exception
            if the types don't match.
        */
    template <class T>
    T & get()
    {
        vigra_precondition(bool(handle_), "Any::get(): object empty.");
        TypedHandle<T> * ptr = dynamic_cast<TypedHandle<T> *>(handle_.get());
        vigra_precondition(ptr != 0, "Any::get(): object is not an instance of the target type.");
        return ptr->value_;
    }

        /** Read-only access to the contained value. This throws an exception
            if the types don't match.
        */
    template <class T>
    T const & get() const
    {
        vigra_precondition(bool(handle_), "Any::get(): object empty.");
        TypedHandle<T> const * ptr = dynamic_cast<TypedHandle<T> const *>(handle_.get());
        vigra_precondition(ptr != 0, "Any::get(): object is not an instance of the target type.");
        return ptr->value_;
    }

        /** Convert the contained value to the given type. This throws an exception
            if <tt>T</tt> or the contained type are not arithmetic types.
        */
    template <class T>
    T cast() const
    {
        vigra_precondition(bool(handle_), "Any::cast(): object empty.");
        ConvertibleHandle const * ptr = dynamic_cast<ConvertibleHandle const *>(handle_.get());
        vigra_precondition(ptr != 0, "Any::cast(): object is not covertible to the target type.");
        return ptr->cast(ConvertibleHandle::TypeTag<T>());
    }
};

} // namespace vigra

#endif // VIGRA_ANY_HXX
