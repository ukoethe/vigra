/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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


#ifndef VIGRA_RANDOM_HXX
#define VIGRA_RANDOM_HXX

#include "mathutil.hxx"
#include "functortraits.hxx"

namespace vigra {

namespace detail {

template<UInt32 StateLength>
struct RandomState;

    /* Tempered twister according to TT800 by M. Matsumoto */
template<>
struct RandomState<25>
{
    static const UInt32 N = 25, M = 7;
    
    mutable UInt32 state_[N];
    mutable UInt32 current_;
                   
    RandomState()
    : current_(0)
    {
        UInt32 seeds[N] = { 
	        0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
	        0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
	        0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
	        0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
	        0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
        };
         
        for(UInt32 i=1; i<N; ++i)
            state_[i] = seeds[i];
    }
    
    void seed(UInt32 theSeed)
    {
        state_[0] = theSeed;
        for(UInt32 i=1; i<N; ++i)
        {
            state_[i] = 1812433253UL * (state_[i-1] ^ (state_[i-1] >> 30)) + i;
        }
    }

  protected:  

    UInt32 get() const
    {
        if(current_ == N)
            generateNumbers();
            
        UInt32 y = state_[current_++];
        y ^= (y << 7) & 0x2b5b2500; 
        y ^= (y << 15) & 0xdb8b0000; 
        y ^= (y >> 16);
        return y;
    }
    
    void generateNumbers() const;
};

void RandomState<25>::generateNumbers() const
{
    UInt32 mag01[2]= { 0x0, 0x8ebfd028 };

    for(UInt32 i=0; i<N-M; ++i)
    {
    	state_[i] = state_[i+M] ^ (state_[i] >> 1) ^ mag01[state_[i] % 2];
    }
    for (UInt32 i=M; i<N; ++i) 
    {
        state_[i] = state_[i+(M-N)] ^ (state_[i] >> 1) ^ mag01[state_[i] % 2];
    }
    current_ = 0;
}

    /* Mersenne twister according to MT19937 by M. Matsumoto */
template<>
struct RandomState<624>
{
    static const UInt32 N = 624, M = 397;
    
    mutable UInt32 state_[N];
    mutable UInt32 current_;
                   
    RandomState()
    : current_(0)
    {
        seed(0xDEADBEEF);
    }
    
    void seed(UInt32 theSeed)
    {
        state_[0] = theSeed;
        for(UInt32 i=1; i<N; ++i)
        {
            state_[i] = 1812433253UL * (state_[i-1] ^ (state_[i-1] >> 30)) + i;
        }
    }

  protected:  

    UInt32 get() const
    {
        if(current_ == N)
            generateNumbers();
            
        UInt32 x = state_[current_++];
        x ^= (x >> 11);
        x ^= (x << 7) & 0x9D2C5680UL;
        x ^= (x << 15) & 0xEFC60000UL;
        return x ^ (x >> 18);
    }
    
    void generateNumbers() const;

    static UInt32 twiddle(UInt32 u, UInt32 v) 
    {
        return (((u & 0x80000000UL) | (v & 0x7FFFFFFFUL)) >> 1)
                ^ ((v & 1UL) ? 0x9908B0DFUL : 0x0UL);
    }

};

void RandomState<624>::generateNumbers() const
{
    for (int i = 0; i < (N - M); ++i)
    {
        state_[i] = state_[i + M] ^ twiddle(state_[i], state_[i + 1]);
    }
    for (int i = N - M; i < (N - 1); ++i)
    {
        state_[i] = state_[i + M - N] ^ twiddle(state_[i], state_[i + 1]);
    }
    state_[N - 1] = state_[M - 1] ^ twiddle(state_[N - 1], state_[0]);
    current_ = 0;
}

} // namespace detail

template <UInt32 StateLength = 25>
class RandomNumberGenerator
: public detail::RandomState<StateLength>
{
    mutable UInt32 normalCurrent_;
    mutable double normalState_;
    
  public:
  
    RandomNumberGenerator()
    : normalCurrent_(0),
      normalState_(0.0)
    {}

    RandomNumberGenerator(UInt32 theSeed)
    : normalCurrent_(0),
      normalState_(0.0)
    {
        seed(theSeed);
    }

        // in [0, 2^32)
    UInt32 operator()() const
    {
        return get();
    }

        // in [0,beyond)
    UInt32 uniformInt(UInt32 beyond) const
    {
        if(beyond < 2)
            return 0;

        UInt32 factor = factorForUniformInt(beyond);
        UInt32 res = get() / factor;

        // Use rejection method to avoid quantization bias.
        // On average, we will need two raw random numbers to generate one.
        while(res >= beyond)
            res = get() / factor;
        return res;
    }
    
        // in [0,1)
    double uniform53() const
    {
	    // make full use of the entire 53-bit mantissa of a double, by Isaku Wada
	    return ( (get() >> 5) * 67108864.0 + (get() >> 6)) * (1.0/9007199254740992.0); 
    }
    
        // in [0,1]
    double uniform() const
    {
        return (double)get() / 4294967295.0;
    }

        // in [lower, upper]
    double uniform(double lower, double upper) const
    {
        vigra_precondition(lower < upper,
          "RandomNumberGenerator::uniform(): lower bound must be smaller than upper bound."); 
        return uniform() * (upper-lower) + lower;
    }
    
    double normal() const
    {
        if(normalCurrent_ == 0)
        {
            normalCurrent_ = 1;

            double x1, x2, w;
            do 
            {
                 x1 = uniform(-1.0, 1.0);
                 x2 = uniform(-1.0, 1.0);
                 w = x1 * x1 + x2 * x2;
            } 
            while ( w >= 1.0 || w == 0.0);
            w = std::sqrt( -2.0 * std::log( w )  / w );
            normalState_ = x2 * w;
            return x1 * w;
            
        }
        else
        {
            normalCurrent_ = 0;
            return normalState_;
        }
    }
    
    
    double normal(double mean, double stddev) const
    {
        vigra_precondition(stddev > 0.0,
          "RandomNumberGenerator::normal(): standard deviation must be positive."); 
        return normal()*stddev + mean;
    }
    
    static RandomNumberGenerator & global()
    {
        static RandomNumberGenerator generator;
        return generator;
    }

    static UInt32 factorForUniformInt(UInt32 range)
    {
        return (range > 2147483648UL || range == 0)
                     ? 1
                     : 2*(2147483648UL / ceilPower2(range));
    }
};

typedef RandomNumberGenerator<25>  TT800; 
typedef RandomNumberGenerator<624> MT19937;

template <UInt32 StateLength>
class FunctorTraits<RandomNumberGenerator<StateLength> >
{
  public:
    typedef RandomNumberGenerator<StateLength> type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

template <class Engine =TT800>
class UniformIntRandomFunctor
{
    UInt32 lower_, difference_, factor_;
    Engine & generator_;
    bool useLowBits_;

  public:
  
    typedef UInt32 argument_type;
    typedef UInt32 result_type;

    UniformIntRandomFunctor(Engine & generator = Engine::global(),
                            bool useLowBits = false)
    : lower_(0), difference_(0xffffffff), factor_(1),
      generator_(generator),
      useLowBits_(useLowBits)
    {}
    
    UniformIntRandomFunctor(UInt32 lower, UInt32 upper, 
                            Engine & generator = Engine::global(),
                            bool useLowBits = false)
    : lower_(lower), difference_(upper-lower), 
      factor_(Engine::factorForUniformInt(difference_ + 1)),
      generator_(generator),
      useLowBits_(useLowBits)
    {
        vigra_precondition(lower < upper,
          "UniformIntRandomFunctor(): lower bound must be smaller than upper bound."); 
    }
    
    UInt32 operator()() const
    {
        if(difference_ == 0xffffffff) // lower_ is necessarily 0
            return generator_();
        else if(useLowBits_)
            return generator_() % (difference_ + 1) + lower_;
        else
        {
            UInt32 res = generator_() / factor_;

            // Use rejection method to avoid quantization bias.
            // On average, we will need two raw random numbers to generate one.
            while(res > difference_)
                res = generator_() / factor_;
            return res + lower_;
        }
    }

        /** Return a uniformly distributed integer random number such that
            <tt>0 <= i < beyond</tt>. This is a required interface for 
            <tt>std::random_shuffle</tt>. It ignores the limits specified 
            in the constructor.
        */
    UInt32 operator()(UInt32 beyond) const
    {
        if(beyond < 2)
            return 0;

        if(useLowBits_)
            return generator_() % beyond;
        else
            return generator_.uniformInt(beyond);
    }
};

template <class Engine>
class FunctorTraits<UniformIntRandomFunctor<Engine> >
{
  public:
    typedef UniformIntRandomFunctor<Engine> type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

template <class Engine =TT800>
class UniformRandomFunctor
{
    double offset_, scale_;
    Engine & generator_;

  public:
  
    typedef double result_type;

    UniformRandomFunctor(Engine & generator = Engine::global())
    : offset_(0.0),
      scale_(1.0),
      generator_(generator)
    {}

    UniformRandomFunctor(double lower, double upper, 
                         Engine & generator = Engine::global())
    : offset_(lower),
      scale_(upper - lower),
      generator_(generator)
    {
        vigra_precondition(lower < upper,
          "UniformRandomFunctor(): lower bound must be smaller than upper bound."); 
    }
    
    double operator()() const
    {
        return generator_.uniform() * scale_ + offset_;
    }
};

template <class Engine>
class FunctorTraits<UniformRandomFunctor<Engine> >
{
  public:
    typedef UniformRandomFunctor<Engine> type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

template <class Engine = TT800>
struct NormalRandomFunctor
{
    double mean_, stddev_;
    Engine & generator_;

  public:
  
    typedef double result_type;

    NormalRandomFunctor(Engine & generator = Engine::global())
    : mean_(0.0),
      stddev_(1.0),
      generator_(generator)
    {}

    NormalRandomFunctor(double mean, double stddev, 
                        Engine & generator = Engine::global())
    : mean_(mean),
      stddev_(stddev),
      generator_(generator)
    {
        vigra_precondition(stddev > 0.0,
          "NormalRandomFunctor(): standard deviation must be positive."); 
    }
    
    double operator()() const
    {
        return generator_.normal() * stddev_ + mean_;
    }

};

template <class Engine>
class FunctorTraits<NormalRandomFunctor<Engine> >
{
  public:
    typedef UniformRandomFunctor<Engine>  type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};


} // namespace vigra 

#endif // VIGRA_RANDOM_HXX
