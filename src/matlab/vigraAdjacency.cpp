/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <string>
#include <set>
#include <vigra/matlab.hxx>
#include <vigra/contourcirculator.hxx>
#include <vigra/pixelneighborhood.hxx>
#include <vigra/diff2d.hxx>
#include <vigra/cellconfigurations.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/multi_pointoperators.hxx>

using namespace vigra;
using namespace matlab;


#define cP2_(a, b) cP<a, b>::value
template <class T>
static void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs)
{
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    //Change these as for your needs.
    typedef UInt32 outType;
    typedef double seedType;

    //Load Input Image
    MultiArrayView<3,T>         in3D        = inputs.getMultiArray<3,T>(0, v_required());
    BasicImageView<T>           in          = makeBasicImageView(in3D.bindOuter(0));
    int                         numOfDim    = inputs.getDimOfInput(0, v_required());


    SparseArray<Int32>            adj_matrix;
    bool                        IsSet_maxRegion;
    UInt32                          max_region_label    = inputs.getScalar<UInt32> ("max_region_label", v_optional(IsSet_maxRegion));
    if(!IsSet_maxRegion)
    {
        FindMinMax<T> minmax;
        inspectMultiArray(srcMultiArrayRange(in3D), minmax);
        max_region_label = static_cast<UInt32>(minmax.max);
        adj_matrix.assign(max_region_label, max_region_label);
    }

    bool                        hasWatershedPixel   = inputs.getBool("hasWatershedPixel", v_default(false));

    //Get Crack options
    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/
    // cantorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function)
    switch(cantorPair(numOfDim, hasWatershedPixel))
    {
        //in case <XX, 0 (crackedges)> the Pixel/Voxel underneath,
        //    on the right (and behind) are considered for neighborhood.
        //in case<IMAG, 1>  the cell configuration around current pixel is compared to
        //    the ones in cellconfigurations.h (Configuration are Freeman coded (see PhD thesis, Chap 5.7)
        // if current configuration is of Type Line then adjacent pixels are considered as neighbors.
        // in case<VOLUME, 1> the number of adjacent regions are counted. the regions are considered
        //    adjacent if only two such regions exist
        case cP2_(IMAG, 0):
        {

            ImageIterator<T> upleft = in.upperLeft();
            ImageIterator<T> downright = in.lowerRight();
            Diff2D size = downright - upleft;
            for(int y = 0; y < size.y-1; ++y){
                for(int x = 0; x< size.x-1; ++x){
                    adj_matrix(upleft(x,y)-1, upleft(x+1, y)-1)++;
                    adj_matrix(upleft(x,y)-1, upleft(x, y+1)-1)++;
                    adj_matrix(upleft(x+1, y)-1, upleft(x,y)-1)++;
                    adj_matrix(upleft(x, y+1)-1, upleft(x,y)-1)++;
                }
            }

            for(int x = 0; x < size.x-1; ++x){
                adj_matrix(upleft(x, size.y-1)-1, upleft(x+1, size.y-1)-1)++;
                adj_matrix(upleft(x+1, size.y-1)-1, upleft(x, size.y-1)-1)++;
            }
            for(int y = 0; y < size.y-1; ++y){
                adj_matrix(upleft(size.x-1, y)-1, upleft(size.x-1, y+1)-1)++;
                adj_matrix(upleft(size.x-1, y+1)-1, upleft(size.x-1, y)-1)++;
            }
            break;
        }
        case cP2_(VOLUME, 0):
        {

            MultiArrayShape<3>::type sze = in3D.shape();
            for(int ii = 0; ii < sze[0]-1; ii++){
                for(int jj = 0; jj < sze[1]-1; jj++){
                    for(int kk = 0; kk < sze[2]-1; kk++){
                        adj_matrix(in3D(ii,jj,kk)-1, in3D(ii+1,jj,kk)-1)++;
                        adj_matrix(in3D(ii,jj,kk)-1, in3D(ii,jj+1,kk)-1)++;
                        adj_matrix(in3D(ii,jj,kk)-1, in3D(ii,jj,kk+1)-1)++;

                        adj_matrix(in3D(ii+1,jj,kk)-1, in3D(ii,jj,kk)-1)++;
                        adj_matrix(in3D(ii,jj+1,kk)-1, in3D(ii,jj,kk)-1)++;
                        adj_matrix(in3D(ii,jj,kk+1)-1, in3D(ii,jj,kk)-1)++;
                    }
                }
            }
            for(int jj = 0; jj < sze[1]-1; jj++){
                for(int kk = 0; kk < sze[2]-1; kk++){
                    adj_matrix(in3D(sze[0]-1,jj,kk)-1, in3D(sze[0]-1,jj+1,kk)-1)++;
                    adj_matrix(in3D(sze[0]-1,jj,kk)-1, in3D(sze[0]-1,jj,kk+1)-1)++;

                    adj_matrix(in3D(sze[0]-1,jj+1,kk)-1, in3D(sze[0]-1,jj,kk)-1)++;
                    adj_matrix(in3D(sze[0]-1,jj,kk+1)-1, in3D(sze[0]-1,jj,kk)-1)++;
                }
            }
            for(int ii = 0; ii < sze[0]-1; ii++){
                for(int kk = 0; kk < sze[2]-1; kk++){
                    adj_matrix(in3D(ii,sze[1]-1,kk)-1, in3D(ii+1,sze[1]-1,kk)-1)++;
                    adj_matrix(in3D(ii,sze[1]-1,kk)-1, in3D(ii,sze[1]-1,kk+1)-1)++;

                    adj_matrix(in3D(ii+1,sze[1]-1,kk)-1, in3D(ii,sze[1]-1,kk)-1)++;
                    adj_matrix(in3D(ii,sze[1]-1,kk+1)-1, in3D(ii,sze[1]-1,kk)-1)++;
                }
            }
            for(int jj = 0; jj < sze[1]-1; jj++){
                for(int ii = 0; ii <sze[0] -1; ii++){
                    adj_matrix(in3D(ii,jj,sze[2]-1)-1, in3D(ii,jj+1,sze[2]-1)-1)++;
                    adj_matrix(in3D(ii,jj,sze[2]-1)-1, in3D(ii+1,jj,sze[2]-1)-1)++;

                    adj_matrix(in3D(ii,jj+1,sze[2]-1)-1, in3D(ii,jj,sze[2]-1)-1)++;
                    adj_matrix(in3D(ii+1,jj,sze[2]-1)-1, in3D(ii,jj,sze[2]-1)-1)++;
                }
            }

            for(int ii = 0; ii < sze[0]-1; ii++){
                adj_matrix(in3D(ii,sze[1]-1,sze[2]-1)-1, in3D(ii+1,sze[1]-1,sze[2]-1)-1)++;
                adj_matrix(in3D(ii+1,sze[1]-1,sze[2]-1)-1, in3D(ii,sze[1]-1,sze[2]-1)-1)++;
            }
            for(int jj = 0; jj < sze[1]-1; jj++){
                adj_matrix(in3D(sze[0]-1,jj,sze[2]-1)-1, in3D(sze[0]-1,jj+1,sze[2]-1)-1)++;
                adj_matrix(in3D(sze[0]-1,jj+1,sze[2]-1)-1, in3D(sze[0]-1,jj,sze[2]-1)-1)++;
            }
            for(int kk = 0; kk < sze[2]-1; kk++){
                adj_matrix(in3D(sze[0]-1,sze[1]-1,kk)-1, in3D(sze[0]-1,sze[1]-1,kk+1)-1)++;
                adj_matrix(in3D(sze[0]-1,sze[1]-1,kk+1)-1, in3D(sze[0]-1,sze[1]-1,kk)-1)++;
            }
            break;
        }
        case cP2_(IMAG, 1):
        {
            ImageIterator<T> upleft = in.upperLeft();
            ImageIterator<T> downright = in.lowerRight();
            Diff2D size = downright - upleft;
            for(int y = 1; y < size.y-1; ++y){
                for(int x = 1; x< size.x-1; ++x){
                    if(upleft(x,y) == 0){
                        NeighborhoodCirculator<typename BasicImageView<T>::Iterator, EightNeighborCode>
                                                            circulator(upleft+ Diff2D(x, y));
                        //circulator += 2;
                        NeighborhoodCirculator<typename BasicImageView<T>::Iterator, EightNeighborCode>
                                                            end(circulator);
                        unsigned char BitField = 0;
                        std::set<T> regions;

                        do{
                            if(*circulator == 0){
                                BitField = BitField >> 1;
                                BitField = BitField | 0x80;
                            }else{
                                regions.insert(*circulator);
                                BitField = BitField >> 1;
                                BitField = BitField | 0x00;
                            }
                        }while(++circulator != end);

                        if(cellimage::cellConfigurations[(int)BitField] == cellimage::CellTypeLine){
                            if(regions.size() == 2){
                                std::set<T>::const_iterator iter = regions.begin();
                                adj_matrix(*iter-1,*(iter++)-1)++;
                            }
                        }
                    }

                    std::ostringstream  oight;
                }
            }
        }
            break;
        case cP2_(VOLUME, 1):
        {
            char ne[26*3] =     {   1,  -1, 0,
            1,  -1, 1,
            0,  -1, 1,
            -1, -1, 1,
            -1, -1, 0,
            -1, -1, -1,
            0,  -1, -1,
            1,  -1, -1,
            1,  0,  0,
            1,  0,  1,
            0,  0,  1,
            -1, 0,  1,
            -1, 0,  0,
            -1, 0,  -1,
            0,  0,  -1,
            1,  0,  -1,
            1,  1,  0,
            1,  1,  1,
            0,  1,  1,
            -1, 1,  1,
            -1, 1,  0,
            -1, 1,  -1,
            0,  1,  -1,
            1,  1,  -1,
            0,  1,  0,
            0,  -1, 0}  ;
            std::set<T> regions;
            MultiArrayShape<3>::type sze = in3D.shape();
            for(int ii = 1; ii < sze[0]-1; ii++)
            {
                for(int jj = 1; jj < sze[1]-1; jj++)
                {
                    for(int kk = 1; kk < sze[2]-1; kk++)
                    {
                        if(in3D(ii,jj,kk) == 0)
                        {


                            regions.clear();

                            for(int ll = 0; ll < 26; ll++){
                                if(in3D(ii+ne[ll*3], jj+ne[ll*3+1], kk+ne[ll*3+2]) != 0 ){
                                    regions.insert((in3D(ii+ne[ll*3], jj+ne[ll*3+1], kk+ne[ll*3+2])));
                                }
                            }

                        if(regions.size() == 2)
                            {
                                std::set<T>::const_iterator iter = regions.begin();
                                adj_matrix(*iter-1,*(iter++)-1)++;
                            }
                        }
                    }
                }
            }

            break;
        }
    }
    for(int ii = 0; ii < max_region_label; ii++){
        for(int jj = ii; jj < max_region_label; jj++){
            if(ii == jj)adj_matrix(ii,ii) = 1;
            else if(adj_matrix.get(ii, jj) != 0)
                adj_matrix(jj,ii) = adj_matrix.get(ii,jj);
            else if(adj_matrix.get(jj, ii) != 0)
                    adj_matrix(ii, jj) = adj_matrix.get(jj, ii);
        }
    }
    adj_matrix.mapToMxArray(outputs[0]);
}


/***************************************************************************************************
**         VIGRA GATEWAY                                                                          **
****************************************************************************************************/
void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    /*
    FLEXIBLE_TYPE_START(0, in);
        ALLOW_D;
    FLEXIBLE_TYPE_END;
    */
    //Add classes as you feel
    //mxClassID inClass;
    //FLEX_TYPE(inClass, 0, in);
    switch(inputs.typeOf(0))
    {
        ALLOW_FD
	ALLOW_UINT_8_64
	ALLOW_INT_8_64
        default:
	    mexErrMsgTxt("Type of input 0 not supported");
    }
}
/** MATLAB
function D = vigraAdjacency(inputArray)
function D = vigraAdjacency(inputArray, options);

D = vigraAdjacency(inputArray) computes the Adjacencymatrix of Label images.
D = vigraAdjacency(inputImage, options)  does the same with user options.
options is a struct with possible fields: "hasWatershedPixel"

D               is a sparse matrix of size max_region_label x max_region_label.
                The entries in D correlate to the length of region borders.
                (Images with and without watershedPixels return different values)
inputArray          must be a Image or a Volume with regions labeled with positive whole numbers
                0 denotes watershed Pixels.
hasWatershedPixel:  it is advised to set this attribute. Otherwise the Function searches for 0 in the
                image.

Usage:
    opt = struct('method' ,value);
    out = vigraAdjacency(in, opt);

*/
