#ifndef VIGRANUMPYSCALEPARAM_HXX
#define VIGRANUMPYSCALEPARAM_HXX

#include <Python.h>
#include <boost/python.hpp>
#include <vigra/tinyvector.hxx>
#include <vigra/tuple.hxx>

namespace python = boost::python;

namespace vigra
{

template <unsigned ndim>
struct pythonScaleParam1
{
    typedef TinyVector<double, ndim> p_vector;
    typedef typename p_vector::const_iterator return_type;
    p_vector vec;
    static unsigned len_check(python::object val, const char *const function_name)
    {
          unsigned count = python::len(val);
          if (count == 1)
              return 0;
          if (count == ndim)
              return 1;

          std::string msg = std::string(function_name)
	                    + "(): Parameter number be 1 or equal to the number of spatial dimensions.";
          PyErr_SetString(PyExc_ValueError, msg.c_str());
          python::throw_error_already_set();
          return 0;
    }
    pythonScaleParam1()
    {}
    pythonScaleParam1(python::object val, const char *const function_name = "pythonScaleParam1")
    {
        if (PySequence_Check(val.ptr()))
        {
            unsigned increment = len_check(val, function_name);
            unsigned i_v = 0;
            for (unsigned i = 0; i != ndim; ++i, i_v += increment)
                vec[i] = python::extract<double>(val[i_v]);
        }
	else
        {
            double x = python::extract<double>(val);
            vec = p_vector(x);
        }
    }
    return_type operator()() const
    {
        return vec.begin();
    }
};

template <unsigned ndim>
struct pythonScaleParam
{
    typedef ConvolutionOptions<ndim> return_type;
    typedef TinyVector<double, ndim> p_vector;
    pythonScaleParam1<ndim> sigma_eff;
    pythonScaleParam1<ndim> sigma_d;
    pythonScaleParam1<ndim> step_size;
    pythonScaleParam1<ndim> outer_scale;
    pythonScaleParam(python::object s_eff, python::object s_d, python::object s_size,
                     const char *const function_name = "pythonScaleParam")
        : sigma_eff(s_eff, function_name),
          sigma_d(s_d, function_name),
          step_size(s_size, function_name)
    {}
    pythonScaleParam(python::object s_eff, python::object s_d, python::object s_size,
                     python::object s_outer, const char *const function_name = "pythonScaleParam")
        : sigma_eff(s_eff, function_name),
          sigma_d(s_d, function_name),
          step_size(s_size, function_name),
	  outer_scale(s_outer, function_name)
    {}
    return_type operator()() const
    {
        return_type opt;
	return opt.stdDev(sigma_eff()).resolutionStdDev(sigma_d())
	          .stepSize(step_size()).outerScale(outer_scale());
    }
};

	
} // namespace vigra

#endif // VIGRANUMPYSCALEPARAM_HXX
