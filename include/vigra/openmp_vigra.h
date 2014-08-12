/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe, Christoph L. Spiel */
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

#ifndef OPENMP_VIGRA_H_INCLUDED_
#define OPENMP_VIGRA_H_INCLUDED_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vigra/diff2d.hxx>
#include <vigra/initimage.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/transformimage.hxx>
#include <vigra/combineimages.hxx>
#include <vigra/convolution.hxx>
#include <vigra/distancetransform.hxx>

#include "openmp_def.h"

#ifdef _MSC_VER
#  define VIGRA_RESTRICT __restrict
#else
#  define VIGRA_RESTRICT __restrict__
#endif


namespace vigra
{
    namespace omp
    {
#ifdef OPENMP
        template <class SrcImageIterator1, class SrcAccessor1,
                  class SrcImageIterator2, class SrcAccessor2,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        combineTwoImages(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                         SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                         DestImageIterator dest_upperleft, DestAccessor dest_acc,
                         const Functor& functor)
        {
            int num_threads = 1;
#pragma omp parallel
            {
                num_threads = omp_get_num_threads();
                const vigra::Size2D size(src1_lowerright - src1_upperleft);
                Functor f(functor);

#pragma omp for schedule(guided) nowait
                for (int y = 0; y < size.y; ++y)
                {
                    const vigra::Diff2D begin(0, y);
                    const vigra::Diff2D end(size.x, y + 1);

                    vigra::combineTwoImages(src1_upperleft + begin, src1_upperleft + end, src1_acc,
                                            src2_upperleft + begin, src2_acc,
                                            dest_upperleft + begin, dest_acc,
                                            f);
                }
            } // omp parallel
            return num_threads;
        }


        template <class SrcImageIterator1, class SrcAccessor1,
                  class SrcImageIterator2, class SrcAccessor2,
                  class MaskImageIterator, class MaskAccessor,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        combineTwoImagesIf(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                           SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                           MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
                           DestImageIterator dest_upperleft, DestAccessor dest_acc,
                           const Functor& functor)
        {
            int num_threads = 1;
#pragma omp parallel
            {
                num_threads = omp_get_num_threads();
                const vigra::Size2D size(src1_lowerright - src1_upperleft);
                Functor f(functor);

#pragma omp for schedule(guided) nowait
                for (int y = 0; y < size.y; ++y)
                {
                    const vigra::Diff2D begin(0, y);
                    const vigra::Diff2D end(size.x, y + 1);

                    vigra::combineTwoImagesIf(src1_upperleft + begin, src1_upperleft + end, src1_acc,
                                              src2_upperleft + begin, src2_acc,
                                              mask_upperleft + begin, mask_acc,
                                              dest_upperleft + begin, dest_acc,
                                              f);
                }
            } // omp parallel
            return num_threads;
        }


        template <class SrcImageIterator1, class SrcAccessor1,
                  class SrcImageIterator2, class SrcAccessor2,
                  class SrcImageIterator3, class SrcAccessor3,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        combineThreeImages(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                           SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                           SrcImageIterator3 src3_upperleft, SrcAccessor3 src3_acc,
                           DestImageIterator dest_upperleft, DestAccessor dest_acc,
                           const Functor& functor)
        {
            int num_threads = 1;
#pragma omp parallel
            {
                num_threads = omp_get_num_threads();
                const vigra::Size2D size(src1_lowerright - src1_upperleft);
                Functor f(functor);

#pragma omp for schedule(guided) nowait
                for (int y = 0; y < size.y; ++y)
                {
                	const vigra::Diff2D begin(0, y);
                    const vigra::Diff2D end(size.x, y + 1);

                    vigra::combineThreeImages(src1_upperleft + begin, src1_upperleft + end, src1_acc,
                                              src2_upperleft + begin, src2_acc,
                                              src3_upperleft + begin, src3_acc,
                                              dest_upperleft + begin, dest_acc,
                                              f);
                }
            } // omp parallel
            return num_threads;
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        inline int
        copyImage(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                  DestImageIterator dest_upperleft, DestAccessor dest_acc)
        {
            int num_threads = 1;
#pragma omp parallel
            {
                num_threads = omp_get_num_threads();
                const vigra::Size2D size(src_lowerright - src_upperleft);

#pragma omp for schedule(guided) nowait
                for (int y = 0; y < size.y; ++y)
                {
                    const vigra::Diff2D begin(0, y);
                    const vigra::Diff2D end(size.x, y + 1);

                    vigra::copyImage(src_upperleft + begin, src_upperleft + end, src_acc,
                                     dest_upperleft + begin, dest_acc);
                }
            } // omp parallel
            return num_threads;
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        transformImage(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                       DestImageIterator dest_upperleft, DestAccessor dest_acc,
                       const Functor& functor)
        {
            int num_threads = 1;
#pragma omp parallel
            {
                num_threads = omp_get_num_threads();
                const vigra::Size2D size(src_lowerright - src_upperleft);
                Functor f(functor);

#pragma omp for schedule(guided) nowait
                for (int y = 0; y < size.y; ++y)
                {
                    const vigra::Diff2D begin(0, y);
                    const vigra::Diff2D end(size.x, y + 1);

                    vigra::transformImage(src_upperleft + begin, src_upperleft + end, src_acc,
                                          dest_upperleft + begin, dest_acc,
                                          f);
                }
            } // omp parallel
            return num_threads;
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class MaskImageIterator, class MaskAccessor,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        transformImageIf(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                         MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
                         DestImageIterator dest_upperleft, DestAccessor dest_acc,
                         const Functor& functor)
        {
            int num_threads = 1;
#pragma omp parallel
            {
                num_threads = omp_get_num_threads();
                const vigra::Size2D size(src_lowerright - src_upperleft);
                Functor f(functor);

#pragma omp for schedule(guided) nowait
                for (int y = 0; y < size.y; ++y)
                {
                    const vigra::Diff2D begin(0, y);
                    const vigra::Diff2D end(size.x, y + 1);

                    vigra::transformImageIf(src_upperleft + begin, src_upperleft + end, src_acc,
                                            mask_upperleft + begin, mask_acc,
                                            dest_upperleft + begin, dest_acc,
                                            f);
                }
            } // omp parallel
            return num_threads;
        }


        namespace fh
        {
            namespace detail
            {
                template <class ValueType>
                inline static ValueType
                square(ValueType x)
                {
                    return x * x;
                }


                // Pedro F. Felzenszwalb, Daniel P. Huttenlocher
                // "Distance Transforms of Sampled Functions"


                template <class ValueType>
                struct ChessboardTransform1D
                {
                    typedef ValueType value_type;

                    int id() const {return 0;}

                    void operator()(ValueType* VIGRA_RESTRICT d, const ValueType* VIGRA_RESTRICT f, int n) const
                    {
                        vigra_fail("fh::detail::ChessboardTransform1D: not implemented");
                    }
                };


                template <class ValueType>
                struct ManhattanTransform1D
                {
                    typedef ValueType value_type;

                    int id() const {return 1;}

                    void operator()(ValueType* VIGRA_RESTRICT d, const ValueType* VIGRA_RESTRICT f, int n) const
                    {
                        const ValueType one = static_cast<ValueType>(1);

                        d[0] = f[0];
                        for (int q = 1; q < n; ++q)
                        {
                            d[q] = std::min<ValueType>(f[q], d[q - 1] + one);
                        }
                        for (int q = n - 2; q >= 0; --q)
                        {
                            d[q] = std::min<ValueType>(d[q], d[q + 1] + one);
                        }
                    }
                };


                template <class ValueType>
                struct EuclideanTransform1D
                {
                    typedef ValueType value_type;

                    int id() const {return 2;}

                    void operator()(ValueType* VIGRA_RESTRICT d, const ValueType* VIGRA_RESTRICT f, int n) const
                    {
                        typedef float math_t;

                        const math_t infinity = std::numeric_limits<math_t>::infinity();

                        int* v = static_cast<int*>(::omp::malloc(n * sizeof(int)));
                        math_t* z = static_cast<math_t*>(::omp::malloc((n + 1) * sizeof(math_t)));
                        int k = 0;

                        v[0] = 0;
                        z[0] = -infinity;
                        z[1] = infinity;

                        for (int q = 1; q < n; ++q)
                        {
                            const math_t sum_q = static_cast<math_t>(f[q]) + square(static_cast<math_t>(q));
                            math_t s = (sum_q - (f[v[k]] + square(v[k]))) / (2 * (q - v[k]));

                            while (s <= z[k])
                            {
                                --k;
                                // IMPLEMENTATION NOTE
                                //     Prefetching improves performance because we must iterate from high to
                                //     low addresses, i.e. against the cache's look-ahead algorithm.
                                HINTED_PREFETCH(z + k - 2U, PREPARE_FOR_READ, HIGH_TEMPORAL_LOCALITY);
                                s = (sum_q - (f[v[k]] + square(v[k]))) / (2 * (q - v[k]));
                            }
                            ++k;

                            v[k] = q;
                            z[k] = s;
                            z[k + 1] = infinity;
                        }

                        k = 0;
                        for (int q = 0; q < n; ++q)
                        {
                            while (z[k + 1] < static_cast<math_t>(q))
                            {
                                ++k;
                            }
                            d[q] = square(q - v[k]) + f[v[k]];
                        }

                        ::omp::free(z);
                        ::omp::free(v);
                    }
                };


                template <class SrcImageIterator, class SrcAccessor,
                          class DestImageIterator, class DestAccessor,
                          class ValueType, class Transform1dFunctor>
                int
                fhDistanceTransform(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor sa,
                                    DestImageIterator dest_upperleft, DestAccessor da,
                                    ValueType background, Transform1dFunctor transform1d)
                {
                    typedef typename Transform1dFunctor::value_type DistanceType;
                    typedef typename vigra::NumericTraits<DistanceType> DistanceTraits;
                    typedef vigra::BasicImage<DistanceType> DistanceImageType;

                    const vigra::Size2D size(src_lowerright - src_upperleft);
                    const int greatest_length = std::max(size.x, size.y);
                    DistanceImageType intermediate(size, vigra::SkipInitialization);

                    int num_threads = 1;
#pragma omp parallel
                    {
                        num_threads = omp_get_num_threads();
                        DistanceType* const f = new DistanceType[greatest_length];
                        DistanceType* const d = new DistanceType[greatest_length];

                        DistanceType* const pf_end = f + size.y;
                        const DistanceType* const pd_end = d + size.y;

                        // IMPLEMENTATION NOTE
                        //     We need "guided" schedule to reduce the waiting time at the
                        //     (implicit) barriers.  This holds true for the next OpenMP
                        //     parallelized "for" loop, too.
#pragma omp for schedule(guided)
                        for (int x = 0; x < size.x; ++x)
                        {
                            SrcImageIterator si(src_upperleft + vigra::Diff2D(x, 0));
                            for (DistanceType* pf = f; pf != pf_end; ++pf)
                            {
                                *pf = EXPECT_RESULT(sa(si) == background, false) ? DistanceTraits::max() : DistanceTraits::zero();
                                ++si.y;
                            }

                            transform1d(d, f, size.y);

                            typename DistanceImageType::column_iterator ci(intermediate.columnBegin(x));
                            for (const DistanceType* pd = d; pd != pd_end; ++pd)
                            {
                                *ci = *pd;
                                ++ci;
                                // IMPLEMENTATION NOTE
                                //     Prefetching about halves the number of stalls per instruction of this loop.
                                HINTED_PREFETCH(ci.operator->(), PREPARE_FOR_WRITE, HIGH_TEMPORAL_LOCALITY);
                            }
                        }

#pragma omp for nowait schedule(guided)
                        for (int y = 0; y < size.y; ++y)
                        {
                            transform1d(d, &intermediate(0, y), size.x);
                            DestImageIterator i(dest_upperleft + vigra::Diff2D(0, y));

                            if (transform1d.id() == 2)
                            {
                                for (DistanceType* pd = d; pd != d + size.x; ++pd, ++i.x)
                                {
                                    da.set(sqrt(*pd), i);
                                }
                            }
                            else
                            {
                                for (DistanceType* pd = d; pd != d + size.x; ++pd, ++i.x)
                                {
                                    da.set(*pd, i);
                                }
                            }
                        }

                        delete [] d;
                        delete [] f;
                    } // omp parallel
                    return num_threads;
                }
            } // namespace detail
        } // namespace fh


        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class ValueType>
        int
        distanceTransform(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor sa,
                          DestImageIterator dest_upperleft, DestAccessor da,
                          ValueType background, int norm)
        {
            switch (norm)
            {
            case 0:
                return fh::detail::fhDistanceTransform(src_upperleft, src_lowerright, sa,
                                                       dest_upperleft, da,
                                                       background,
                                                       fh::detail::ChessboardTransform1D<float>());
                break;

            case 1:
                return fh::detail::fhDistanceTransform(src_upperleft, src_lowerright, sa,
                                                       dest_upperleft, da,
                                                       background,
                                                       fh::detail::ManhattanTransform1D<float>());
                break;

            case 2: // FALLTHROUGH
            default:
                return fh::detail::fhDistanceTransform(src_upperleft, src_lowerright, sa,
                                                       dest_upperleft, da,
                                                       background,
                                                       fh::detail::EuclideanTransform1D<float>());
            }
        }


#else


        template <class SrcImageIterator1, class SrcAccessor1,
                  class SrcImageIterator2, class SrcAccessor2,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        combineTwoImages(SrcImageIterator1 src1_upperleft,
                         SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                         SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                         DestImageIterator dest_upperleft, DestAccessor dest_acc,
                         const Functor& func)
        {
            vigra::combineTwoImages(src1_upperleft, src1_lowerright, src1_acc,
                                    src2_upperleft, src2_acc,
                                    dest_upperleft, dest_acc,
                                    func);
            return 1;
        }


        template <class SrcImageIterator1, class SrcAccessor1,
                  class SrcImageIterator2, class SrcAccessor2,
                  class MaskImageIterator, class MaskAccessor,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        combineTwoImagesIf(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                           SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                           MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
                           DestImageIterator dest_upperleft, DestAccessor dest_acc,
                           const Functor& func)
        {
            vigra::combineTwoImagesIf(src1_upperleft, src1_lowerright, src1_acc,
                                      src2_upperleft, src2_acc,
                                      mask_upperleft, mask_acc,
                                      dest_upperleft, dest_acc,
                                      func);
            return 1;
        }


        template <class SrcImageIterator1, class SrcAccessor1,
                  class SrcImageIterator2, class SrcAccessor2,
                  class SrcImageIterator3, class SrcAccessor3,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        combineThreeImages(SrcImageIterator1 src1_upperleft, SrcImageIterator1 src1_lowerright, SrcAccessor1 src1_acc,
                           SrcImageIterator2 src2_upperleft, SrcAccessor2 src2_acc,
                           SrcImageIterator3 src3_upperleft, SrcAccessor3 src3_acc,
                           DestImageIterator dest_upperleft, DestAccessor dest_acc,
                           const Functor& func)
        {
            vigra::combineThreeImages(src1_upperleft, src1_lowerright, src1_acc,
                                      src2_upperleft, src2_acc,
                                      src3_upperleft, src3_acc,
                                      dest_upperleft, dest_acc,
                                      func);
            return 1;
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        inline int
        copyImage(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                  DestImageIterator dest_upperleft, DestAccessor dest_acc)
        {
            vigra::copyImage(src_upperleft, src_lowerright, src_acc, dest_upperleft, dest_acc);
            return 1;
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        transformImage(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                       DestImageIterator dest_upperleft, DestAccessor dest_acc,
                       const Functor& func)
        {
            vigra::transformImage(src_upperleft, src_lowerright, src_acc,
                                  dest_upperleft, dest_acc,
                                  func);
            return 1;
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class MaskImageIterator, class MaskAccessor,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        transformImageIf(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                         MaskImageIterator mask_upperleft, MaskAccessor mask_acc,
                         DestImageIterator dest_upperleft, DestAccessor dest_acc,
                         const Functor& func)
        {
            vigra::transformImageIf(src_upperleft, src_lowerright, src_acc,
                                    mask_upperleft, mask_acc,
                                    dest_upperleft, dest_acc,
                                    func);
            return 1;
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class ValueType>
        int
        distanceTransform(SrcImageIterator src_upperleft, SrcImageIterator src_lowerright, SrcAccessor src_acc,
                          DestImageIterator dest_upperleft, DestAccessor dest_acc,
                          ValueType background, int norm)
        {
            vigra::distanceTransform(src_upperleft, src_lowerright, src_acc,
                                     dest_upperleft, dest_acc,
                                     background, norm);
            return 1;
        }

#endif // OPENMP


        //
        // Argument Object Factory versions
        //

        template <class SrcImageIterator1, class SrcAccessor1,
                  class SrcImageIterator2, class SrcAccessor2,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        combineTwoImages(vigra::triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                         vigra::pair<SrcImageIterator2, SrcAccessor2> src2,
                         vigra::pair<DestImageIterator, DestAccessor> dest,
                         const Functor& functor)
        {
            return vigra::omp::combineTwoImages(src1.first, src1.second, src1.third,
                                                src2.first, src2.second,
                                                dest.first, dest.second,
                                                functor);
        }


        template <class SrcImageIterator1, class SrcAccessor1,
                  class SrcImageIterator2, class SrcAccessor2,
                  class MaskImageIterator, class MaskAccessor,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        combineTwoImagesIf(vigra::triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                           vigra::pair<SrcImageIterator2, SrcAccessor2> src2,
                           vigra::pair<MaskImageIterator, MaskAccessor> mask,
                           vigra::pair<DestImageIterator, DestAccessor> dest,
                           const Functor& functor)
        {
            return vigra::omp::combineTwoImagesIf(src1.first, src1.second, src1.third,
                                                  src2.first, src2.second,
                                                  mask.first, mask.second,
                                                  dest.first, dest.second,
                                                  functor);
        }                                         


        template <class SrcImageIterator1, class SrcAccessor1,
                  class SrcImageIterator2, class SrcAccessor2,
                  class SrcImageIterator3, class SrcAccessor3,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        combineThreeImages(vigra::triple<SrcImageIterator1, SrcImageIterator1, SrcAccessor1> src1,
                           vigra::pair<SrcImageIterator2, SrcAccessor2> src2,
                           vigra::pair<SrcImageIterator3, SrcAccessor3> src3,
                           vigra::pair<DestImageIterator, DestAccessor> dest,
                           const Functor& functor)
        {
            return vigra::omp::combineThreeImages(src1.first, src1.second, src1.third,
                                                  src2.first, src2.second,
                                                  src3.first, src3.second,
                                                  dest.first, dest.second,
                                                  functor);
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int
        transformImage(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                       vigra::pair<DestImageIterator, DestAccessor> dest,
                       const Functor& functor)
        {
            return vigra::omp::transformImage(src.first, src.second, src.third,
                                              dest.first, dest.second,
                                              functor);
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        inline int
        copyImage(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                  vigra::pair<DestImageIterator, DestAccessor> dest)
        {
            return vigra::omp::copyImage(src.first, src.second, src.third,
                                         dest.first, dest.second);
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class MaskImageIterator, class MaskAccessor,
                  class DestImageIterator, class DestAccessor,
                  class Functor>
        inline int 
        transformImageIf(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                         vigra::pair<MaskImageIterator, MaskAccessor> mask,
                         vigra::pair<DestImageIterator, DestAccessor> dest,
                         const Functor& functor)
        {
            return vigra::omp::transformImageIf(src.first, src.second, src.third,
                                                mask.first, mask.second,
                                                dest.first, dest.second,
                                                functor);
        }


        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor,
                  class ValueType>
        inline int
        distanceTransform(vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                          vigra::pair<DestImageIterator, DestAccessor> dest,
                          ValueType background, int norm)
        {
            return vigra::omp::distanceTransform(src.first, src.second, src.third,
                                                 dest.first, dest.second,
                                                 background, norm);
        }
    } // namespace omp
} // namespace vigra


#endif // OPENMP_VIGRA_H_INCLUDED_

// Local Variables:
// mode: c++
// End:
