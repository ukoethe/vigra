/************************************************************************/
/*                                                                      */
/*                 Copyright 2014 by Ullrich Koethe                     */
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

#ifndef VIGRA_CONSTRUCT_CUSTODIAN_FOR_HPP
#define VIGRA_CONSTRUCT_CUSTODIAN_FOR_HPP

#include <boost/version.hpp>

#if BOOST_VERSION <= 105500

#include <boost/python/with_custodian_and_ward.hpp>

namespace boost { namespace python { 

template <std::size_t ward>
struct construct_custodian_for : default_call_policies
{
    BOOST_STATIC_ASSERT(ward > 0);
    
    template <class ArgumentPackage>
    static PyObject* postcall(ArgumentPackage const& args_, PyObject* result)
    {
        std::size_t arity_ = detail::arity(args_);
#if BOOST_WORKAROUND(BOOST_MSVC, < 1300)
        if (ward > arity_ )
#else
        // check if ward exceeds the arity
        // (this weird formulation avoids "always false" warnings
        // for arity_ = 0)
        if ( (std::max<std::size_t>)(0, ward) > arity_ )
#endif
        {
            PyErr_SetString(
                PyExc_IndexError
              , "boost::python::construct_custodian_for: argument index out of range"
            );
            return 0;
        }
        
        PyObject* patient = detail::get_prev<ward>::execute(args_, result);
        PyObject* nurse = detail::get_prev<1>::execute(args_.base);

        if (nurse == 0) return 0;
    
        result = default_call_policies::postcall(args_, result);
        if (result == 0)
            return 0;
            
        if (python::objects::make_nurse_and_patient(nurse, patient) == 0)
        {
            Py_XDECREF(result);
            return 0;
        }
        return result;
    }
};

}} // namespace boost::python

#  endif // BOOST_VERSION <= 105500

#endif // VIGRA_CONSTRUCT_CUSTODIAN_FOR_HPP
