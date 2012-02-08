#ifndef VIGRA_BIT_ARRAY_HXX
#define VIGRA_BIT_ARRAY_HXX

#include <functional>
#include <ostream>
#include "metaprogramming.hxx"

namespace vigra {

template <class> // undefined class to provoke usable error messages
class vigra_error_BitArray_accepts_only_unsigned_underlying_types_and_no_;

template <class X> // bitwise operators do not necessarily work for bool
struct EnableBitArray
    : public enable_if<HasMetaLog2<X>::value && !IsSameType<X, bool>::value> {};

// BitArray: a minimal subset of std::bitset with the extension of compile-time
// access functions set<unsigned>(), test<unsigned>(), reset<unsigned>(), and
// flip<unsigned>(), plus all relational operators;
// furthermore, there are no range checks.

template <unsigned SIZE, class WORD_TYPE = unsigned, class = void>
class BitArray
    : public
      vigra_error_BitArray_accepts_only_unsigned_underlying_types_and_no_
      <WORD_TYPE>
{};

template <unsigned SIZE, class WORD_TYPE>
class BitArray<SIZE, WORD_TYPE, typename EnableBitArray<WORD_TYPE>::type>
{
    // 'unsigned' will be the most efficent word type for most CPUs,
    // since very long immediates such as a possible 64 bit 'unsigned long'
    // are slower for many typical uses of BitArray
  protected:
    static const unsigned bit_size = SIZE;
    static const unsigned word_len = MetaLog2<WORD_TYPE>::value;
    static const unsigned array_len = (bit_size + word_len - 1) / word_len;
    static const unsigned last_pos = array_len - 1;
    template <unsigned pos>
    struct bit_index
    {
        static const unsigned  word_pos = pos / word_len;
        static const unsigned   bit_pos = pos % word_len;
        static const WORD_TYPE bit_mask = WORD_TYPE(1) << bit_pos;
    };
    typedef bit_index<bit_size> size_index;
    static const WORD_TYPE ones_mask = ~(WORD_TYPE(0));
    static const unsigned border_pos = size_index::bit_pos;
    static const WORD_TYPE last_mask = !border_pos ? 0
                                                   : size_index::bit_mask - 1;
    static const bool does_fit = border_pos == 0;
    unsigned word_pos(unsigned pos) const
    {
        return pos / word_len;
    };
    WORD_TYPE bit_mask(unsigned pos) const
    {
        return WORD_TYPE(1) << (pos % word_len); // the compiler knows as well..
    };

    WORD_TYPE set_bits[array_len];

  public:
    unsigned size()
    {
        return bit_size;
    }
    void clear()
    {
        for (unsigned i = 0; i != array_len; ++i)
            set_bits[i] = 0;
    }
    BitArray()
    {
        clear();
    }
    template <unsigned pos>
    void set()
    {
        typedef bit_index<pos> index;
        set_bits[index::word_pos] |= index::bit_mask;
    }
    template <unsigned pos>
    void reset()
    {
        typedef bit_index<pos> index;
        set_bits[index::word_pos] &= ~index::bit_mask;
    }
    template <unsigned pos>
    void flip()
    {
        typedef bit_index<pos> index;
        set_bits[index::word_pos] ^= index::bit_mask;
    }
    template <unsigned pos>
    bool test() const
    {
        typedef bit_index<pos> index;
        return (set_bits[index::word_pos] & index::bit_mask) != 0;
    }

    BitArray & set(unsigned pos, bool value = true)
    {
        (set_bits[word_pos(pos)] &= ~bit_mask(pos))
                                 |= value ? bit_mask(pos) : 0;
        return *this;
    }
    BitArray & reset(unsigned pos)
    {
        set_bits[word_pos(pos)] &= ~bit_mask(pos);
        return *this;
    }
    BitArray & flip(unsigned pos)
    {
        set_bits[word_pos(pos)] ^= bit_mask(pos);
        return *this;
    }
    bool test(unsigned pos) const
    {
        return set_bits[word_pos(pos)] & bit_mask(pos);
    }
    bool operator[](unsigned pos) const
    {
        return test(pos);
    }

    BitArray & set()
    {
        for (unsigned i = 0; i != last_pos + does_fit; ++i)
            set_bits[i] = ones_mask;
        if (!does_fit)
            set_bits[last_pos] = last_mask;
        return *this;
    }
    BitArray & reset()
    {
        for (unsigned i = 0; i != array_len; ++i)
            set_bits[i] = 0;
        return *this;
    }
    BitArray & flip()
    {
        for (unsigned i = 0; i != last_pos + does_fit; ++i)
            set_bits[i] ^= ones_mask;
        if (!does_fit)
            set_bits[last_pos] ^= last_mask;
        return *this;
    }

    operator bool() const
    {
        for (unsigned i = 0; i != array_len; ++i)
            if (set_bits[i] != 0)
                return true;
        return false;
    }
    bool operator!() const
    {
        return !bool(*this);
    }
    bool any() const
    {
        return *this;
    }
    bool none() const
    {
        return !*this;
    }
    bool all() const
    {
        for (unsigned i = 0; i != last_pos + does_fit; ++i)
            if (set_bits[i] != ones_mask)
                return false;
        if (!does_fit)
            return set_bits[last_pos] == last_mask;
        return true;
    }
    
    BitArray operator~() const
    {
        BitArray x(*this);
        x.flip();
        return x;
    }
   
  protected:
    template <class F>
    bool mutual_compare(const BitArray & t, F f, bool if_equal = false) const
    {
        for (int i = last_pos; i >= 0; i--)
        {
            WORD_TYPE x =   set_bits[i];
            WORD_TYPE y = t.set_bits[i];
            if (f(x, y))
                return true;
            if (f(y, x))
                return false;
        }
        return if_equal;
    }
    typedef std::less<WORD_TYPE>    less;
    typedef std::greater<WORD_TYPE> greater;
    
  public:
    bool operator<(const BitArray & t) const
    {
        return mutual_compare(t, less());
    }
    bool operator>(const BitArray & t) const
    {
        return mutual_compare(t, greater());
    }

    bool operator<=(const BitArray & t) const
    {
        return mutual_compare(t, less(), true);
    }
    bool operator>=(const BitArray & t) const
    {
        return mutual_compare(t, greater(), true);
    }

    bool operator!=(const BitArray & t) const
    {
        for (unsigned i = 0; i != array_len; ++i)
            if (set_bits[i] != t.set_bits[i])
                return true;
        return false;
    }
    bool operator==(const BitArray & t) const
    {
        return !operator!=(t);
    }

  protected:
    struct bit_and_assign
    {
        static void assign(WORD_TYPE & a, WORD_TYPE b) { a &= b; }
    };
    struct exclusive_or_assign
    {
        static void assign(WORD_TYPE & a, WORD_TYPE b) { a ^= b; }
    };
    struct bit_or_assign
    {
        static void assign(WORD_TYPE & a, WORD_TYPE b) { a |= b; }
    };
    template <class A>
    BitArray & assign_operator(const BitArray & x)
    {
        for (unsigned i = 0; i != array_len; ++i)
            A::assign(set_bits[i], x.set_bits[i]);
        return *this;
    }
  public:
    BitArray & operator&=(const BitArray & x)
    {
        return assign_operator<bit_and_assign>(x);
    }
    BitArray & operator^=(const BitArray & x)
    {
        return assign_operator<exclusive_or_assign>(x);
    }
    BitArray & operator|=(const BitArray & x)
    {
        return assign_operator<bit_or_assign>(x);
    }
   
  protected:
    template <class A>
    BitArray & bit_operator(const BitArray & y) const
    {
        BitArray x(*this);
        return x.assign_operator<A>(y);
    }
  public:
    BitArray operator&(const BitArray & y) const
    {
        return bit_operator<bit_and_assign>(y);
    }
    BitArray operator^(const BitArray & y) const
    {
        return bit_operator<exclusive_or_assign>(y);
    }
    BitArray operator|(const BitArray & y) const
    {
        return bit_operator<bit_or_assign>(y);
    }

    bool operator&&(const BitArray & y) const
    {
        return *this && y;
    }
    bool operator||(const BitArray & y) const
    {
        return *this || y;
    }

    friend std::ostream & operator<<(std::ostream & os, const BitArray & z)
    {
        for (int i = bit_size - 1; i >= 0; i--)
            os << (z[i] ? "1" : "0");
        return os;
    }
};

// work around GCC's zero-sized array extension
template <class WORD_TYPE>
class BitArray<0, WORD_TYPE>
{
    bool error[-(long int)sizeof(WORD_TYPE)];
};

} // namespace vigra

#endif // VIGRA_BIT_ARRAY_HXX
