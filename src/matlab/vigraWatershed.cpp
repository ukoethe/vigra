/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/


#include <vigra/matlab.hxx>

#include <vigra/watersheds.hxx>
#include <vigra/watersheds3d.hxx>

#include <vigra/seededregiongrowing.hxx>
#include <vigra/seededregiongrowing3d.hxx>



using namespace vigra;
using namespace matlab;


//#define RN_DEBUG
#define cP3_(a, b , c) cP3<a, b, c>::value
template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    //Change these as for your needs.
    typedef UInt32 outType;
    typedef UInt32 seedType;

    //Load Input Image
    MultiArrayView<3,T>         in3D        = inputs.getMultiArray<3,T>(0, v_required());
    BasicImageView<T>           in          = makeBasicImageView(in3D.bindOuter(0));
    Int32                       numOfDim    = inputs.getDimOfInput(0, v_required());

    //Load seed if provided
    enum Methods{UNION = 0, SEED = 1} ;
    Methods                            method       = UNION;
    MultiArrayView<3,seedType>         seed3D       = inputs.getMultiArray<3,seedType>("seeds", v_optional());
    BasicImageView<seedType>           seed         = makeBasicImageView(seed3D.bindOuter(0));
    if(seed3D.data()!= 0 && !(seed3D.shape() == in3D.shape()))
        mexErrMsgTxt("Seed Array and Input Array Dimension/ Number of Value Mismatch");
    if(seed3D.data()!= 0 )  method = SEED;


    //Get right connectivity options
    Int32                         v2Dconn[2]      = {8, 4};
    Int32                         v3Dconn[2]      = {26, 6};
    Int32                         connectivity    = (method == UNION)?
                                inputs.getScalarVals2D3D<int>("conn", v_default(8,26,numOfDim),
                                                                                    v2Dconn, v2Dconn+2,
                                                                                    v3Dconn, v3Dconn+2,
                                                                                    numOfDim)
                            :   inputs.getScalarVals2D3D<int>("conn", v_default(4,6,numOfDim),
                                                                                    v2Dconn+1, v2Dconn+2,
                                                                                    v3Dconn+1, v3Dconn+2,
                                                                                    numOfDim);
    //Get Crack options
    VIGRA_CREATE_ENUM_AND_STD_MAP2(Crack,CrackMap, complete_grow, keep_contours);
    Crack                       crack           = (Crack)inputs.getEnum("crack", v_default(UInt32(complete_grow)), CrackMap);
    SRGType                     SRGcrack        = (crack == complete_grow)? vigra::CompleteGrow : vigra::KeepContours;


    double                      CostThreshold   =  inputs.getScalar<double>("CostThreshold", v_default(-1.0));

    //Allocate space for outputArray.
    MultiArrayView<3,outType>   out3D           = outputs.createMultiArray      <3,outType>   (0, v_required(), in3D.shape());
    BasicImageView<outType>     out(out3D.data(), in3D.shape(0), in3D.shape(1));
    // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function)

    UInt32 max_region_label = 0;


    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/
    // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function)
    switch(cantorPair(method, numOfDim, connectivity)){
        //cP is the templated version o f the cantorPair function first value is Dimension of Inputimage, second the connectivity setting
        //Code is basically the code on the VIGRA-reference page
        case cP3_(UNION, IMAG, 8):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("UNION IMAG 8");
            #endif
            max_region_label = watersheds(srcImageRange(in), destImage(out));

            break;
        case cP3_(UNION, IMAG, 4):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("UNION IMAG 4");
            #endif
            max_region_label = watersheds(srcImageRange(in), destImage(out), FourNeighborCode());

            break;
        case cP3_(UNION, VOLUME, 26):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("UNION VOLUME 26");
            #endif
            max_region_label = watersheds3DTwentySix(srcMultiArrayRange(in3D), destMultiArray(out3D));

            break;
        case cP3_(UNION, VOLUME, 6):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("UNION VOLUME 6");
            #endif
            max_region_label = watersheds3DSix(srcMultiArrayRange(in3D), destMultiArray(out3D));

            break;
        case cP3_(SEED, IMAG, 4):
        {
            #ifdef RN_DEBUG
            mexWarnMsgTxt("SEED IMAG 4");
            #endif
            //find maximimum of seed Image
            FindMinMax<seedType> minmax;   // init functor
            inspectImage(srcImageRange(seed), minmax);
            max_region_label = static_cast<UInt32>(minmax.max);

            ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<seedType> >
                                            gradstat(max_region_label);

            seededRegionGrowing(    srcImageRange(in),
                                    srcImage(seed),
                                    destImage(out),
                                    gradstat,SRGcrack );
            break;
        }
        case cP3_(SEED, VOLUME, 6):
        {
            #ifdef RN_DEBUG
            mexWarnMsgTxt("SEED VOLUME 6");
            #endif
            FindMinMax<seedType> minmax;
            inspectMultiArray(srcMultiArrayRange(seed3D), minmax);
            max_region_label = static_cast<UInt32>(minmax.max);

            ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<seedType> >
                                            gradstat(max_region_label);

            seededRegionGrowing3D(  srcMultiArrayRange(in3D),
                                    srcMultiArray(seed3D),
                                    destMultiArray(out3D),
                                    gradstat,CostThreshold, SRGcrack);
            break;
        }
        default:
            mexErrMsgTxt("Something went horribily wrong - Internal error");
    }
    outputs.createScalar(1,v_optional(), max_region_label);
}


void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    /*
    FLEXIBLE_TYPE_START(0, in);
        ALLOW_D;
    FLEXIBLE_TYPE_END;
    */

    mxClassID inClass;
    FLEX_TYPE(inClass, 0, in);
    switch(inClass)
    {
        ALLOW_UINT_8;
        DEFAULT_ERROR;
    }
}

/** MATLAB
function L = vigraWatershed(inputArray)
function L = vigraWatershed(inputArray, options);

L = vigraWatershed(inputArray);
        Uses the union find algorithm to compute a label matrix identifying
        the watershed regions of the inputArray (which may be 2 or 3 dimensional).
        The elements of L are Int32 values greater than 0. Boundaries are
        represented by crack edges, i.e. there are no explicit watershed pixels.

options = struct('seeds', seedImage, 'crack' = 'keepContours');
L = vigraWatershed(inputImage, options);
        Uses seeded region growing to compute a label matrix identifying
        the watershed regions of the inputArray (which may be 2 or 3 dimensional).
        The elements of L are Int32 values greater than 0. Boundaries are
        represented by explicit watershed pixels, i.e. there is a one-pixel
        thick boundary between every pair of regions.

inputArray - a 2D or 3D array of type 'single' or 'double'

options    - is a struct with the following possible fields
    'seeds':    An Array of same size as inputArray. If supplied, seeded region growing
                shall be used. If not, the union find algorithm will be used. As of now,
                the seed array has to be of same type as the input array - this will be changed in the next update.
    'conn':     2D: 4 (default),  8 (only supported by the union find algorithm)
                3D: 6 (default), 26 (only supported by the union find algorithm)
                While using seeded region growing, only the default values can be used.
    'crack':    'completeGrow' (default), 'keepContours' (only supported by seeded region growing and 3D Volumes)
                Choose whether to keep watershed pixels or not. While using union find,
                only the default value can be used.
    'CostThreshold':  -1.0 (default) - any double value.
		If, at any point in the algorithm, the cost of the current candidate exceeds the optional
		max_cost value (which defaults to -1), region growing is aborted, and all voxels not yet
		assigned to a region remain unlabeled.
Usage:
    opt = struct('fieldname' ,'value',....);
    out = vigraWatershed(in, opt);
*/


