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


/* Mersenne twister according to TT800 by M. Matsumoto */
class RandomNumberGenerator
{
    static const UInt32 stateLength = 25,
                        M           = 7, 
                        upperMask   = 1<<31, 
                        lowerMask   = ~upperMask, 
                        mask        = 0xffffffff;
    
    mutable UInt32 state_[stateLength],
                   current_, normalCurrent_;
    mutable double normalState_;
    
  public:
  
    RandomNumberGenerator()
    : current_(0),
      normalCurrent_(0),
      normalState_(0.0)
    {
        UInt32 seeds[stateLength] = { 
	        0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
	        0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
	        0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
	        0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
	        0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
        };
         
        for(UInt32 i=1; i<stateLength; ++i)
            state_[i] = seeds[i];
    }

    RandomNumberGenerator(UInt32 theSeed)
    : current_(0),
      normalCurrent_(0),
      normalState_(0.0)
    {
        seed(theSeed);
    }

    UInt32 operator()() const
    {
        if(current_ == stateLength)
        {
            generateNumbers();
            current_ = 0;
        }
        UInt32 y = state_[current_++];
        y ^= (y << 7) & 0x2b5b2500; 
        y ^= (y << 15) & 0xdb8b0000; 
        y ^= (y >> 16);
        return y;
    }
    
    double uniform() const
    {
        return (double)operator()() / double(mask);
    }

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
    
    static RandomNumberGenerator & globalGenerator()
    {
        static RandomNumberGenerator generator;
        return generator;
    }
    
    void seed(UInt32 theSeed = 5489)
    {
        state_[0] = theSeed;
        for(UInt32 i=1; i<stateLength; ++i)
        {
            state_[i] = (1812433253UL * (state_[i-1] ^ (state_[i-1] >> 30)) + i) & mask;
        }
    }

  private:

    void generateNumbers() const
    {
        UInt32 mag01[2]= { 0x0, 0x8ebfd028 };
        
        for(UInt32 i=0; i<stateLength-M; ++i)
        {
    	    state_[i] = state_[i+M] ^ (state_[i] >> 1) ^ mag01[state_[i] % 2];
        }
        for (UInt32 i=M; i<stateLength; ++i) 
        {
            state_[i] = state_[i+(M-stateLength)] ^ (state_[i] >> 1) ^ mag01[state_[i] % 2];
        }
    }
 
};

template <>
class FunctorTraits<RandomNumberGenerator>
{
  public:
    typedef RandomNumberGenerator type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

class UniformIntRandomFunctor
{
    UInt32 lower_, difference_, factor_;
    RandomNumberGenerator & generator_;
    bool useLowBits_;

  public:
  
    typedef UInt32 argument_type;
    typedef UInt32 result_type;

    UniformIntRandomFunctor(RandomNumberGenerator & generator = RandomNumberGenerator::globalGenerator(),
                            bool useLowBits = false)
    : lower_(0), difference_(0xffffffff), factor_(1),
      generator_(generator),
      useLowBits_(useLowBits)
    {}
    
    UniformIntRandomFunctor(UInt32 lower, UInt32 upper, 
                            RandomNumberGenerator & generator = RandomNumberGenerator::globalGenerator(),
                            bool useLowBits = false)
    : lower_(lower), difference_(upper-lower), factor_(computeFactor(difference_ + 1)),
      generator_(generator),
      useLowBits_(useLowBits)
    {
        vigra_precondition(lower < upper,
          "UniformIntRandomFunctor(): lower bound must be smaller than upper bound."); 
    }
    
    UInt32 operator()() const
    {
        UInt32 res = generator_();
        if(useLowBits_)
        {
            if(difference_ < 0xffffffff)
                res %= difference_ + 1;
        }
        else
        {
            res /= factor_;

            // Use rejection method to avoid quantization bias.
            // On average, we will need two raw random numbers to generate one.
            while(res > difference_)
                res = generator_() / factor_;
        }
        return res + lower_;
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

        UInt32 res = generator_();
        if(useLowBits_)
        {
            res %= beyond;
        }
        else
        {
            UInt32 factor = computeFactor(beyond);

            res /= factor;

            // Use rejection method to avoid quantization bias.
            // On average, we will need two raw random numbers to generate one.
            while(res >= beyond)
                res = generator_() / factor;
        }
        return res;
    }
    
    static UInt32 computeFactor(UInt32 range)
    {
        return (range > 2147483648UL || range == 0)
                     ? 1
                     : 2*(2147483648UL / ceilPower2(range));
    }
};

template <>
class FunctorTraits<UniformIntRandomFunctor>
{
  public:
    typedef UniformIntRandomFunctor type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

class UniformRandomFunctor
{
    double offset_, scale_;
    RandomNumberGenerator & generator_;

  public:
  
    typedef double result_type;

    UniformRandomFunctor(RandomNumberGenerator & generator = RandomNumberGenerator::globalGenerator())
    : offset_(0.0),
      scale_(1.0),
      generator_(generator)
    {}

    UniformRandomFunctor(double lower, double upper, 
                         RandomNumberGenerator & generator = RandomNumberGenerator::globalGenerator())
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

template <>
class FunctorTraits<UniformRandomFunctor>
{
  public:
    typedef UniformRandomFunctor type;
    
    typedef VigraTrueType  isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

struct NormalRandomFunctor
{
    double mean_, stddev_;
    RandomNumberGenerator & generator_;

  public:
  
    typedef double result_type;

    NormalRandomFunctor(RandomNumberGenerator & generator = RandomNumberGenerator::globalGenerator())
    : mean_(0.0),
      stddev_(1.0),
      generator_(generator)
    {}

    NormalRandomFunctor(double mean, double stddev, 
                        RandomNumberGenerator & generator = RandomNumberGenerator::globalGenerator())
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

template <>
class FunctorTraits<NormalRandomFunctor>
{
  public:
    typedef UniformRandomFunctor type;
    
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
