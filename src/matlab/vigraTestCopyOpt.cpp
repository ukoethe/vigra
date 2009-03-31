/***************************************************************************************************
**         INCLUDES AND DEFS                                                                      **
****************************************************************************************************/

#include <vigra/matlab.hxx>
#include <vigra/copyimage.hxx>
#include <vigra/multi_pointoperators.hxx>
using namespace vigra;
using namespace matlab;


template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    MultiArrayView<3, T>        Vol         =   inputs.getMultiArray<3, T>("field0", v_required());
    BasicImageView<T>           Img         =   inputs.getImage<T>("field1", v_required());
    TinyVectorView<T, 2>        TinyVec     =   inputs.getTinyVector<T, 2>("field2", v_required());
    T                           Scalar      =   inputs.getScalar<T>("field3", v_required());
    //ConstCellArray              Cell        =   inputs.getCellArray("field4", v_required());

    MultiArrayView<3,T>         outVol      =   outputs.createMultiArray<3, T>(0, v_required(), Vol.shape());
    BasicImageView<T>           outImg      =   outputs.createImage<T>(1, v_required(), Img.width(), Img.height());
    MultiArrayShape<2>::type    newShape3 (2, 1);
    MultiArrayView<2, T>        outTiny     =   outputs.createMultiArray<2, T>(2, v_required(), newShape3);
    T*                          outScal     =   outputs.createScalar<T>(3, v_required());
    //ConstCellArray              OutCell     =   outputs.createCellArray(5, v_required(), Cell.size()+1);


    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/
    vigra::copyMultiArray(srcMultiArrayRange(Vol), destMultiArray(outVol));
    outVol(0,0,0) = outVol(0,0,0)+1;
    vigra::copyImage(srcImageRange(Img), destImage(outImg));
    outImg(0,0) = outImg(0,0)+1;
    outTiny(0,0) = TinyVec[0]+1;
    outTiny(1,0) = TinyVec[1];
    outputs.createScalar<T>(4, v_required(), 42);
    *outScal = Scalar+1;

}

/***************************************************************************************************
**         VIGRA GATEWAY                                                                          **
****************************************************************************************************/
void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
   switch(inputs.typeOf("field0"))
   {
       ALLOW_FD;
       ALLOW_INT;
       ALLOW_UINT;
       default:
            mexErrMsgTxt("Invalid Type");
   }
}
/** MATLAB
NODOC
vigraBjoernsFunc1903(ProbMap, Thres);
*/
