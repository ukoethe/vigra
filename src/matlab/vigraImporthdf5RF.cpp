/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by Rahul Nair and Ullrich Koethe          */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#include <iostream>
#include <set>
#include <vigra/matlab.hxx>
#include "random_forest_impex.hxx"
#include <vigra/random_forest_hdf5_impex.hxx>

using namespace vigra;
using namespace matlab;



void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /* INPUT */
    if (inputs.size() != 1)
        mexErrMsgTxt("One inputs required.");

    // get RF object
   RandomForest<> rf;
   std::string filename = inputs.getString(0, v_required());
   std::string groupname = inputs.getString(1, v_default(std::string("")));
   vigra::rf_import_HDF5(rf, filename, groupname );
   matlab::exportRandomForest(rf, matlab::createCellArray(2*rf.options_.tree_count_+2, outputs[0]));
}




/***************************************************************************************************
**         VIGRA GATEWAY                                                                          **
****************************************************************************************************/
inline void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    vigraMain(outputs, inputs);
};


/**ADDITIONAL_BUILD_FLAGS
-lhdf5 -lhdf5_hl
*/

/** MATLAB
function rf = vigraImporthdf5RF( filename, groupname);

Import a previously trained Random Forest from  a hdf5 file

    rf        - MATLAB cell array representing the random forest classifier
   filename  - name of hdf5 file

   groupname    - optional: name of group which shoud be used as the base
   					path

to compile on Linux:
--------------------
  mex vigraImporthdf5RF.cpp -I../../include -lhdf5 -lhdf5_hl

to compile on Windows:
----------------------
  mex vigraImporthdf5RF.cpp -I../../include -I[HDF5PATH]/include -L[HDF5PATH]/lib -lhdf5dll -lhdf5_hldll -D_HDF5USEDLL_ -DHDF5CPP_USEDLL

hdf5 1.6.x or hdf5 1.8.x must be installed. 
*/
