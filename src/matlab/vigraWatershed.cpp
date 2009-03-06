/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <string>
#include <vigra/watersheds.hxx>
#include <vigra/watersheds3d.hxx>
#include <vigra/matlab.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/seededregiongrowing3d.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just
//set outType to T



using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
template <class T>
struct data
: public base_data<T>
{
    declScalar2D3D(int, conn, 4, 6);
    declCharConstr2(crack, completeGrow, keepContours);
    declOut(Int32);
    declScalar(double, CostThreshold, -1.0);

    enum method_enum {UNION, SEED};
    method_enum method;

    SRGType SRGcrack;
    BasicImageView<T>  seed;
    MultiArrayView<3,T> seed3D;

    data(matlab::OutputArray outputs, matlab::InputArray inputs)
    :           base_data<T>(inputs),
                initOption(conn),
				initOption(crack),
				initOption(CostThreshold),
				method(UNION),
                SRGcrack(crack == completeGrow
                             ? vigra::CompleteGrow
                             : vigra::KeepContours)
    {
        if(this->options.isValid("seeds"))
        {
            method = SEED;
            if(this->numOfDim == IMAG)
            {
                if(conn != 4)
                    mexErrMsgTxt("vigraWatershed(): Connectivity for 2D seeded region growing must be 4.");

                seed = matlab::getImage<T>(this->options["seeds"]);
            }
            if(this->numOfDim == VOLUME && conn != 6)
                mexErrMsgTxt("vigraWatershed(): Connectivity for 3D seeded region growing must be 6.");

            seed3D = matlab::getMultiArray<3, T>(this->options["seeds"]);
            if(seed3D.shape() != this->in3D.shape())
                mexErrMsgTxt("igraWatershed(): Seed and input array dimension mismatch.");
        }
        else
        {
            method = UNION;
            if(this->numOfDim == IMAG && conn != 4 && conn != 8)
                mexErrMsgTxt("vigraWatershed(): Connectivity for 2D union find must be 4 or 8.");
            if(this->numOfDim == VOLUME && conn != 6 && conn != 26)
                mexErrMsgTxt("vigraWatershed(): Connectivity for 3D union find must be 6 or 26.");
        }
        if(CostThreshold != -1.0 && method == SEED)
        {
            this->numOfDim = VOLUME;
            conn = 6;
        }
        initOut_SAME(Int32);
    }
};
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cP3_(a, b , c) cP3<data<T>::a, b, c>::value
struct vigraFunctor
{
    template <class T>
    static void exec(matlab::OutputArray outputs, matlab::InputArray inputs){
        //Options
        data<T>  o(outputs, inputs);

        // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function)
        int max_region_label = -1;
        // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function)
        switch(cantorPair(o.method, o.numOfDim, o.conn)){
            //cP is the templated version o f the cantorPair function first value is Dimension of Inputimage, second the connectivity setting
            //Code is basically the code on the VIGRA-reference page
            case cP3_(UNION, IMAG, 8):
                max_region_label = watersheds(srcImageRange(o.in), destImage(o.out));

                break;
            case cP3_(UNION, IMAG, 4):
                max_region_label = watersheds(srcImageRange(o.in), destImage(o.out), FourNeighborCode());

                break;
            case cP3_(UNION, VOLUME, 26):
                max_region_label = watersheds3DTwentySix(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D));

                break;
            case cP3_(UNION, VOLUME, 6):
                max_region_label = watersheds3DSix(srcMultiArrayRange(o.in3D), destMultiArray(o.out3D));

                break;
            case cP3_(SEED, IMAG, 4):
            {
                //find maximimum of seed Image
                FindMinMax<T> minmax;   // init functor
                inspectImage(srcImageRange(o.seed), minmax);
                max_region_label = minmax.max;

                ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<T> >
                                                gradstat(max_region_label);

                seededRegionGrowing(    srcImageRange(o.in),
                                        srcImage(o.seed),
                                        destImage(o.out),
                                        gradstat,o.SRGcrack );
                break;
            }
            case cP3_(SEED, VOLUME, 6):
            {
                mexWarnMsgTxt("SEED IMAG 4");
                FindMinMax<T> minmax;
                inspectMultiArray(srcMultiArrayRange(o.seed3D), minmax);
                max_region_label = minmax.max;

                ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<T> >
                                                gradstat(max_region_label);

                seededRegionGrowing3D(  srcMultiArrayRange(o.in3D),
                                        srcMultiArray(o.seed3D),
                                        destMultiArray(o.out3D),
                                        gradstat,o.CostThreshold, o.SRGcrack);
                break;
            }
            default:
                mexErrMsgTxt("Something went wrong");
        }


    }
};


/*+++++++++++++++++++++++MexEntryFunc++++++++++++++++++++++++++++++++*/
/* Gatewayfunction - see matlab.hxx for details.
/* if a certain class is NOT supported - you will have to copy the
/* body of the callMexFunctor function and edit it here.
/* Supports (u)int[8|16|32|64], float and double.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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


