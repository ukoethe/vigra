#include <boost/python.hpp>

using namespace boost::python;

void defCanny();
void defVolume();

#include <boost/python.hpp>

BOOST_PYTHON_MODULE_INIT(hourglasscmodule)
{
    // import vigra to make sure that basic converters (e.g. for
    // Size2D) are available for default arguments:
    {
        PyObject *vigraMod = PyImport_ImportModule("vigra");
        Py_XDECREF(vigraMod);
    }
    
    defVolume();
    defCanny();
}
