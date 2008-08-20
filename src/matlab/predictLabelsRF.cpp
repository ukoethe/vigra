#include <iostream>
#include "mex.h"
#include <set>
#include <vigra/random_forest.hxx>

using namespace std;
using namespace vigra;
 
 
 
 void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
                     
   /* INPUT */
   /* CHECK INPUTS */ 
   /* Check NUMBER of in- and outputs */    
    if (nrhs != 2)
        mexErrMsgTxt("Two inputs required.");
    else if (nlhs > 1)
        mexErrMsgTxt("Too many output arguments");
   /* Check type of FEATURE MATRIX (1. Input) */
    if (!(mxIsCell(prhs[0])))
        mexErrMsgTxt("First argument must be a RandomForest cell array.");
   /* Check type of LABEL MATRIX (2. Input) */
    if (!(mxIsDouble(prhs[1])))
        mexErrMsgTxt("Input feature matrix must be of type double.");    
    
    typedef MultiArrayView<2, double>::difference_type Shape;
    
    const mwSize *FeatureDimVector = mxGetDimensions(prhs[1]);
    
    UInt32 numberOfSamples = FeatureDimVector[0];
    UInt32 numberOfFeatures = FeatureDimVector[1];
    MultiArrayView<2, double> features(Shape(numberOfSamples, numberOfFeatures), mxGetPr(prhs[1]));

    // read RF parameters
    int parameterCount = 3;
    mxArray *parameterArray =  mxGetCell(prhs[0], 0);
    UInt32 * parameters = (UInt32 *)mxGetData(parameterArray);
    FeatureDimVector = mxGetDimensions(parameterArray);
    if(FeatureDimVector[0] != 3)
        mexErrMsgTxt("Parameter array must have size 3.");
    UInt32 labelCount = parameters[0];
    UInt32 Ntree = parameters[2]; 
    if(numberOfFeatures != parameters[1])
        mexErrMsgTxt("Feature array has wrong number of columns.");
     
    // read array of possible class labels
    mxArray *classLabelsArray = mxGetCell(prhs[0], 1);
    double *classLabels = mxGetPr(classLabelsArray);
    FeatureDimVector = mxGetDimensions(classLabelsArray);
    if(FeatureDimVector[0] != labelCount)
        mexErrMsgTxt("Class label array has wrong size.");
     
    ArrayVector<ArrayVector<Int32> >  trees;
    ArrayVector<ArrayVector<double> > weights;
    
    // for all decision trees
    for(UInt32 k=0; k<Ntree; ++k)
    {
        // read int tree array
        mxArray *treeArray = mxGetCell(prhs[0], 2*k+2);
        Int32 * tree = (Int32 *)mxGetData(treeArray);
        FeatureDimVector = mxGetDimensions(treeArray);
        int treeSize = FeatureDimVector[0];
        trees.push_back(ArrayVector<Int32>(tree, tree+treeSize));
        
        // read double weight/threshold array
        mxArray *weightArray = mxGetCell(prhs[0], 2*k+3);
        double *weight = mxGetPr(weightArray);
        FeatureDimVector = mxGetDimensions(weightArray);
        int weightSize = FeatureDimVector[0];
        
        weights.push_back(ArrayVector<double>(weight, weight+weightSize));
    }
    
    plhs[0] = mxCreateDoubleMatrix(numberOfSamples, 1, mxREAL);
     
    MultiArrayView<2, double> labels(Shape(numberOfSamples, 1), mxGetPr(plhs[0]));
    
    RandomForest<double> rf(classLabels, classLabels+labelCount, Ntree, numberOfFeatures, trees.begin(), weights.begin());
    
    rf.predictLabels(features, labels);
}
