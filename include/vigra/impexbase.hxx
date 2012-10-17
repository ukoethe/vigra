/************************************************************************/
/*                                                                      */
/*               Copyright 2012 Christoph Spiel                         */
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

#ifndef VIGRA_IMPEXBASE_HXX
#define VIGRA_IMPEXBASE_HXX


#include <string>
#include "inspectimage.hxx"
#include "sized_int.hxx"
#include "utilities.hxx"


// #ifdef NO_UNIQUE_PTR
#if 1
#  define UNIQUE_PTR std::auto_ptr
#else
#  include <memory>
#  define UNIQUE_PTR std::unique_ptr
#endif

namespace vigra
{
    typedef enum
    {
        UNSIGNED_INT_8,
        UNSIGNED_INT_16,
        UNSIGNED_INT_32,
        SIGNED_INT_16,
        SIGNED_INT_32,
        IEEE_FLOAT_32,
        IEEE_FLOAT_64
    } pixel_t;


    namespace detail
    {
        inline static pixel_t
        pixel_t_of_string(const std::string& pixel_type)
        {
            if (pixel_type == "UINT8")
            {
                return UNSIGNED_INT_8;
            }
            else if (pixel_type == "UINT16")
            {
                return UNSIGNED_INT_16;
            }
            else if (pixel_type == "UINT32")
            {
                return UNSIGNED_INT_32;
            }
            else if (pixel_type == "INT16")
            {
                return SIGNED_INT_16;
            }
            else if (pixel_type == "INT32")
            {
                return SIGNED_INT_32;
            }
            else if (pixel_type == "FLOAT")
            {
                return IEEE_FLOAT_32;
            }
            else if (pixel_type == "DOUBLE")
            {
                return IEEE_FLOAT_64;
            }
            else
            {
                vigra_fail("vigra_ext::detail::pixel_t_of_string: unknown pixel type");
                return UNSIGNED_INT_8; // NOT REACHED
            }
        }


        struct identity
        {
            template <typename T>
            T operator()(T x) const
            {
                return x;
            }
        };


        typedef pair<double, double> range_t;


        class linear_transform
        {
        public:
            linear_transform(const range_t& source, const range_t& destination) :
                scale_((destination.second - destination.first) / (source.second - source.first)),
                offset_(destination.first / scale_ - source.first)
            {}

            template <typename T>
            double operator()(T x) const
            {
                return scale_ * static_cast<double>(x) + offset_;
            }

        private:
            const double scale_;
            const double offset_;
        };


        template <class Iterator, class Accessor>
        inline static range_t
        find_value_range(Iterator upper_left, Iterator lower_right, Accessor accessor,
                         /* is_scalar? */ VigraTrueType)
        {
            typedef typename Accessor::value_type value_type;

            FindMinMax<value_type> extrema;

            inspectImage(upper_left, lower_right, accessor, extrema);

            return range_t(static_cast<double>(extrema.min), static_cast<double>(extrema.max));
        }


        template <class Iterator, class Accessor>
        inline static range_t
        find_value_range(Iterator upper_left, Iterator lower_right, Accessor accessor,
                         /* is_scalar? */ VigraFalseType)
        {
            typedef typename Accessor::ElementAccessor element_accessor;
            typedef typename element_accessor::value_type value_type;

            const int number_of_bands(static_cast<int>(accessor.size(upper_left)));
            FindMinMax<value_type> extrema;

            for (int i = 0; i != number_of_bands; ++i)
            {
                element_accessor band(i, accessor);

                inspectImage(upper_left, lower_right, band, extrema);
            }

            return range_t(static_cast<double>(extrema.min), static_cast<double>(extrema.max));
        }


        template <class SourceIterator, class SourceAccessor>
        inline static range_t
        find_source_value_range(const ImageExportInfo& export_info,
                                SourceIterator upper_left, SourceIterator lower_right, SourceAccessor accessor)
        {
            if (export_info.getFromMin() < export_info.getFromMax())
            {
                return range_t(export_info.getFromMin(), export_info.getFromMax());
            }
            else
            {
                typedef typename SourceAccessor::value_type SourceValueType;
                typedef typename NumericTraits<SourceValueType>::isScalar is_scalar;

                const range_t range(find_value_range(upper_left, lower_right, accessor, is_scalar()));

                if (range.first < range.second)
                {
                    return range_t(range.first, range.second);
                }
                else
                {
                    return range_t(range.first, range.first + 1.0);
                }
            }
        }


        template <typename T>
        inline static range_t
        find_destination_value_range(const ImageExportInfo& export_info)
        {
            if (export_info.getToMin() < export_info.getToMax())
            {
                return range_t(export_info.getToMin(), export_info.getToMax());
            }
            else
            {
                return range_t(static_cast<double>(NumericTraits<T>::min()),
                               static_cast<double>(NumericTraits<T>::max()));
            }
        }


        inline static range_t
        find_destination_value_range(const ImageExportInfo& export_info, pixel_t pixel_type)
        {
            switch (pixel_type)
            {
            case UNSIGNED_INT_8: return find_destination_value_range<UInt8>(export_info);
            case UNSIGNED_INT_16: return find_destination_value_range<UInt16>(export_info);
            case UNSIGNED_INT_32: return find_destination_value_range<UInt32>(export_info);
            case SIGNED_INT_16: return find_destination_value_range<Int16>(export_info);
            case SIGNED_INT_32: return find_destination_value_range<Int32>(export_info);
            case IEEE_FLOAT_32: return find_destination_value_range<float>(export_info);
            case IEEE_FLOAT_64: return find_destination_value_range<double>(export_info);
            default:
                vigra_fail("vigra_ext::detail::find_destination_value_range: not reached");
                return range_t(0.0, 0.0); // NOT REACHED
            }
        }
    } // end namespace detail
} // end namespace vigra


#endif // VIGRA_IMPEXBASE_HXX
