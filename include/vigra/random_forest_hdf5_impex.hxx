/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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


#ifndef VIGRA_RANDOM_FOREST_IMPEX_HXX
#define VIGRA_RANDOM_FOREST_IMPEX_HXX
#include "vigra/windows.h"
#include <vigra/random_forest.hxx>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <cstdio>

namespace vigra {

template <class T>
std::auto_ptr<RandomForest<T> >
importHDF5RandomForest(hid_t parent_id)
{
    hid_t       group_id;
    hsize_t     size;
    char        name[100];

    group_id = H5Gopen (parent_id, "RandomForest", H5P_DEFAULT);
    vigra_postcondition(group_id >= 0, "importHDF5RandomForest(): Unable to open 'RandomForest' group - wrong file contents.");
    
    int labelCount;
    vigra_postcondition(H5LTread_dataset_int (group_id, "labelCount", &labelCount) >= 0,
                        "importHDF5RandomForest(): Unable to read labelCount.");
    int featureCount;
    vigra_postcondition(H5LTread_dataset_int (group_id, "featureCount", &featureCount) >= 0,
                        "importHDF5RandomForest(): Unable to read featureCount.");
    int treeCount;
    vigra_postcondition(H5LTread_dataset_int (group_id, "treeCount", &treeCount) >= 0,
                        "importHDF5RandomForest(): Unable to read treeCount.");
                        
    ArrayVector<double> labelSet(labelCount);
    vigra_postcondition(H5LTread_dataset_double (group_id, "labelSet", labelSet.begin()) >= 0,
                        "importHDF5RandomForest(): Unable to read labelSet.");

    ArrayVector<ArrayVector<Int32> >  trees;
    ArrayVector<ArrayVector<double> > weights;
    
    // for all decision trees
    for(Int32 k=1; k <= treeCount; ++k)
    {
        std::sprintf(name, "tree%04d", k);
        vigra_postcondition(H5LTget_dataset_info(group_id, name, &size, NULL, NULL) >= 0,
                "importHDF5RandomForest(): Unable to read tree array size.");
        trees.push_back(ArrayVector<Int32>((ArrayVector<Int32>::size_type)size));
        vigra_postcondition(H5LTread_dataset_int (group_id, name, trees.back().begin()) >= 0,
                "importHDF5RandomForest(): Unable to read tree array.");

        std::sprintf(name, "weights%04d", k);
        vigra_postcondition(H5LTget_dataset_info(group_id, name, &size, NULL, NULL) >= 0,
                "importHDF5RandomForest(): Unable to read weight array size.");
        weights.push_back(ArrayVector<double>((ArrayVector<double>::size_type)size));
        vigra_postcondition(H5LTread_dataset_double (group_id, name, weights.back().begin()) >= 0,
                "importHDF5RandomForest(): Unable to read weight array.");
    }
    
    vigra_postcondition(H5Gclose(group_id) >= 0,
                    "importHDF5RandomForest(): Unable to close group.");
    
    return std::auto_ptr<RandomForest<T> >(
                           new RandomForest<T>(labelSet.begin(), labelSet.end(), 
                               treeCount, featureCount, trees.begin(), weights.begin()));
}

template <class T>
std::auto_ptr<RandomForest<T> >
importHDF5RandomForest(const char* filename)
{
    hid_t       file_id;

    file_id = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    vigra_postcondition(file_id >= 0, "importHDF5RandomForest(): Unable to open file.");

    std::auto_ptr<RandomForest<T> > res(importHDF5RandomForest<T>(file_id));

    vigra_postcondition(H5Fclose(file_id) >= 0,
                    "importHDF5RandomForest(): Unable to close file.");
    
    return res;
}

template <class T>
void exportHDF5RandomForest(RandomForest<T> const & rf, hid_t parent_id)
{
    hid_t   group_id;
    hsize_t size = 1;
    char    name[100];

    group_id = H5Gcreate(parent_id, "/RandomForest", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    vigra_postcondition(group_id >= 0, "exportHDF5RandomForest(): Unable to open 'RandomForest' group.");
       
    int labelCount = rf.labelCount();
    vigra_postcondition(H5LTmake_dataset (group_id, "labelCount", 1, &size, H5T_NATIVE_UINT, &labelCount) >= 0,
                        "exportHDF5RandomForest(): Unable to write labelCount.");

    int featureCount = rf.featureCount();
    vigra_postcondition(H5LTmake_dataset (group_id, "featureCount", 1, &size, H5T_NATIVE_UINT, &featureCount) >= 0,
                        "exportHDF5RandomForest(): Unable to write featureCount.");

    int treeCount = rf.treeCount();
    vigra_postcondition(H5LTmake_dataset (group_id, "treeCount", 1, &size, H5T_NATIVE_UINT, &treeCount) >= 0,
                        "exportHDF5RandomForest(): Unable to write treeCount.");

    size = labelCount;
    vigra_postcondition(H5LTmake_dataset (group_id, "labelSet", 1, &size, H5T_NATIVE_DOUBLE, rf.classes_.begin()) >= 0,
                        "exportHDF5RandomForest(): Unable to write labelSet.");

    // for all decision trees
    for(Int32 k=0; k < treeCount; ++k)
    {
        std::sprintf(name, "tree%04d", k+1);
        size = rf.trees_[k].tree_.size();
        vigra_postcondition(H5LTmake_dataset (group_id, name, 1, &size, H5T_NATIVE_INT, rf.trees_[k].tree_.begin()) >= 0,
                        "exportHDF5RandomForest(): Unable to write tree array.");

        std::sprintf(name, "weights%04d", k+1);
        size = rf.trees_[k].terminalWeights_.size();
        vigra_postcondition(H5LTmake_dataset (group_id, name, 1, &size, H5T_NATIVE_DOUBLE, rf.trees_[k].terminalWeights_.begin()) >= 0,
                        "exportHDF5RandomForest(): Unable to write weight array.");
    }
    vigra_postcondition(H5Gclose(group_id) >= 0,
                    "exportHDF5RandomForest(): Unable to close group.");
}

template <class T>
void exportHDF5RandomForest(RandomForest<T> const & rf, const char * filename)
{
    hid_t   file_id;

    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    vigra_postcondition(file_id >= 0, "exportHDF5RandomForest(): Unable to open file.");
       
    exportHDF5RandomForest(rf, file_id);

    vigra_postcondition(H5Fclose(file_id) >= 0,
                    "exportHDF5RandomForest(): Unable to close file.");
}

} // namespace vigra


#endif // VIGRA_RANDOM_FOREST_IMPEX_HXX

