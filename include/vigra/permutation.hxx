#ifndef VIGRA_PERMUTATION_HXX
#define VIGRA_PERMUTATION_HXX

#include "config.hxx"
#include "error.hxx"
#include "array_vector.hxx"
#include "tinyvector.hxx"

namespace vigra {

template <unsigned int N>
class Permutation : public TinyVector<unsigned int, N>
{
  public:
    typedef TinyVector<unsigned int, N>         base_type;
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::value_type      value_type;
    typedef typename base_type::iterator        iterator;
    typedef typename base_type::const_iterator  const_iterator;
    typedef typename base_type::reference       reference;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::pointer         pointer;
    typedef typename base_type::const_pointer   const_pointer;
    typedef int                                 integral_type;

    Permutation() : base_type() {}

    Permutation(const Permutation<N-1> & other, const size_type index)
    : base_type()
    {
        vigra_precondition(
                index < N,
                "Permutation::Permutation(): Invalid index");
        for (size_type n = 0; n < N; n++)
        {
            if (n < index)
            {
                (*this)[n] = other[n];
            }
            else if (n == index)
            {
                (*this)[n] = N - 1;
            }
            else
            {
                (*this)[n] = other[n-1];
            }
        }
        if ((N - 1 - index) % 2 == 0)
        {
            sign_ = (other.sign() == 1);
        }
        else
        {
            sign_ = (other.sign() == -1);
        }
    }

    integral_type sign() const
    {
        if (sign_)
        {
            return +1;
        }
        else
        {
            return -1;
        }
    }

  private:
    bool sign_;

};

template <>
class Permutation<1> : public TinyVector<unsigned int, 1>
{
  public:
    typedef TinyVector<unsigned int, 1> base_type;
    typedef base_type::size_type        size_type;
    typedef base_type::value_type       value_type;
    typedef base_type::iterator         iterator;
    typedef base_type::const_iterator   const_iterator;
    typedef base_type::reference        reference;
    typedef base_type::const_reference  const_reference;
    typedef base_type::pointer          pointer;
    typedef base_type::const_pointer    const_pointer;
    typedef int                         integral_type;

    Permutation() : base_type()
    {
        (*this)[0] = 0;
        (*this).sign_ = true;
    }

    integral_type sign() const
    {
        if (sign_)
        {
            return +1;
        }
        else
        {
            return -1;
        }
    }

  private:
    bool sign_;
};

template <unsigned int N>
class PlainChangesPermutations : public ArrayVector<Permutation<N> >
{
  public:
    typedef ArrayVector<Permutation<N> >        base_type;
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::value_type      value_type;
    typedef typename base_type::iterator        iterator;
    typedef typename base_type::const_iterator  const_iterator;
    typedef typename base_type::reference       reference;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::pointer         pointer;
    typedef typename base_type::const_pointer   const_pointer;

    PlainChangesPermutations() : base_type()
    {
        PlainChangesPermutations<N-1> permutations;
        for (auto permutation : permutations)
        {
            if (permutation.sign() == -1)
            {
                for (unsigned int n = 0; n < N; n++)
                {
                    this->push_back(Permutation<N>(permutation, n));
                }
            }
            else
            {
                for (unsigned int n = N; n > 0; n--)
                {
                    this->push_back(Permutation<N>(permutation, n - 1));
                }
            }
        }
    }
};

template <>
class PlainChangesPermutations<1> : public ArrayVector<Permutation<1> >
{
  public:
    typedef ArrayVector<Permutation<1> >  base_type;
    typedef base_type::size_type          size_type;
    typedef base_type::value_type         value_type;
    typedef base_type::iterator           iterator;
    typedef base_type::const_iterator     const_iterator;
    typedef base_type::reference          reference;
    typedef base_type::const_reference    const_reference;
    typedef base_type::pointer            pointer;
    typedef base_type::const_pointer      const_pointer;

    PlainChangesPermutations() : base_type()
    {
        this->push_back(Permutation<1>());
    }
};

} /* namespace vigra */

#endif
