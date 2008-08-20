#include <iostream>
#include <vigra/matlab.hxx>
#include <vigra/multi_resize.hxx>

using namespace vigra;
 
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{    
    matlab::InputArgumentArray inputs(nrhs, prhs);
    matlab::OutputArgumentArray outputs(nlhs, plhs);
    
    MultiArrayView<3, double> in = matlab::getMultiArray<3, double>(inputs[0]);    
    MultiArrayShape<3>::type newShape = matlab::getShape<3>(inputs[1]);
    MultiArrayView<3, double> out = matlab::createMultiArray<3, double>(newShape, outputs[0]);

    resizeMultiArraySplineInterpolation(srcMultiArrayRange(in), destMultiArrayRange(out));
#if 0 // cellarray tests    
    matlab::ConstCellArray cellin = matlab::getCellArray(inputs[2]);
    
    std::cerr << matlab::getScalar<double>(cellin[0]) << " " <<
                 matlab::getScalar<float>(cellin[1]) << " " <<
                 matlab::getScalar<Int32>(cellin[2]) << "\n";
    matlab::CellArray cells = matlab::createCellArray(3, outputs[1]);
    cells[0] = matlab::createScalar(1.0);
    MultiArrayView<1, double> cv = matlab::createMultiArray<1, double>(TinyVector<Int32, 1>(3), cells[1]);
    cv[0] = 10.0;
    cv[1] = 11.0;
    cv[2] = 12.0;
    cells[2] = matlab::createScalar(3);
#endif
}
