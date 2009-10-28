#ifndef RN_RF_VARIABLE_IMPORTANCE_HXX
#define RN_RF_VARIABLE_IMPORTANCE_HXX


#include <vigra/array_vector.hxx>
#include <vector>
#include <algorithm>
#include <vigra/random.hxx>
#include <vigra/matrix.hxx>

namespace rn
{
template<class T>
struct FieldProxy
{

    T* something;
    FieldProxy(T& in): something(&in){};

    FieldProxy(): something(0){};
    FieldProxy& operator= (const T rhs)
    {
        *something = rhs;
        return *this;
    }
    FieldProxy& operator= (const FieldProxy rhs)
    {
        *something = *(rhs.something);
        return *this;
    }

    operator T() const
    {
        return *something;
    }


};

template<class T>
void swap(FieldProxy<T> & rhs, FieldProxy<T> & lhs)
{
    T temp= *(lhs.something);
    *(lhs.something) = *(rhs.something);
    *(rhs.something) = temp;
}

//copied into random_forest.hxx
//
//template<class RFClassType, class T1, class C1, class T2, class C2, class T3, class C3>
//void variableImportanceMarginals(vigra::RandomForest<RFClassType> const & rf ,
//                                 vigra::MultiArrayView<2, T1, C1> & samples,
//                                 vigra::MultiArrayView<2, T2, C2>  const& labels)
//{
//    checkPreconditions(rf, samples, labels, imp_array);
//    unsigned int                        sampleCount         = rf.options_.oob_data.shape(0);
//    unsigned int                        treeCount           = rf.options_.oob_data.shape(0);
//    unsigned int                        classCount          = rf.classes_.size();
//    unsigned int                        permutationCount    = 5;
//                                                            //How often should the mariginal permutation
//                                                            //be done.
//    std::vector<double>                 backup_column;
//                                                            //Vector containing the original values
//                                                            //of a column
//    std::vector<FieldProxy<double> >    original_column;
//                                                            //Vector of "Pointers" pointing to original
//                                                            //data (see FieldProxy)
//    std::vector<size_t>                 indices_of_oobData;
//    std::vector<size_t>::iterator       iter;
//                                                            //Vector containing the indices of the oob
//                                                            //data of current tree and iterator
//    vigra::ArrayVector<double>          oob_error(classCount);
//    vigra::UniformIntRandomFunctor<vigra::RandomTT800>
//                                        rand_int;
//    vigra::MultiArrayView<2, double>    imp_array(rf.options_.variable_importance);
//
//    /**
//        Go thru each tree:
//            find oob indices and make a backup
//            Shuffle data permutationCount times and find oob error
//            //copy backup back and find the original ooberror
//            copy backup back -> oroignal oob error is calculated in the RF.learn function;
//    **/
//    for(size_t ii = 0; ii < treeCount; ++ii)
//    {
//        //Find indices of the oob Data of current tree.
//        indices_of_oobData.clear();
//        for(size_t jj = 0; jj < sampleCount; ++jj)
//            if(rf.options_.oob_data(jj, ii) != 0)
//                indices_of_oobData.push_back(jj);
//
//        //get per class and other oob_error for each tree.
//        oob_error.init(0);
//        for(iter = indices_of_oobData.begin(); iter != indices_of_oobData.end(); ++iter)
//        {
//            if(rf.predictLabel(vigra::linalg::rowVector(samples, *iter)) != labels(*iter, 0))
//            {
//                oob_error[classCount] += 1;
//                oob_error[labels(*iter,0)] += 1;
//            }
//        }
//
//        //Find the variable importance for each variable
//        for(size_t kk_feat = 0; kk_feat < rf.columnCount_; ++kk_feat)
//        {
//            //make backup of orinal column
//            backup_column.clear();
//            original_column.clear();
//            // TODO Hopefully the malloc only takes place in the first loop
//            for(iter = indices_of_oobData.begin(); iter != indices_of_oobData.end(); ++iter)
//            {
//                backup_column.push_back(samples(*iter,kk_feat));
//                original_column.push_back(samples(*iter,kk_feat));
//            }
//
//            // TODO do this multiple time for stability (is this a good idea?)
//            for(size_t gg = 0; gg < permutationCount; ++gg)
//            {
//                std::random_shuffle(original_column.begin(), original_column.end(), rand_int);
//                // Check misclassification rate.
//                for(iter = indices_of_oobData.begin(); iter != indices_of_oobData.end(); ++iter)
//                {
//                    if(rf.predictLabel(vigra::linalg::rowVector(samples, *iter)) != labels(*iter, 0))
//                    {
//                        imp_array(kk_feat, classCount) += 1;
//                        imp_array(kk_feat, labels(*iter, 0)) +=1;
//                    }
//                }
//            }
//            //load backup of the column
//            std::copy(backup_column.begin(), backup_column.end(), original_column);
//
//            //normalise;
//            imp_array(kk_feat, classCount) = imp_array(kk_feat, classCount)/permutationCount;
//            imp_array(kk_feat, classCount) = (imp_array(kk_feat, classCount) - oob_error[classCount])/oob_error[classCount];
//            for(size_t gg_class = 0; gg_class < classCount; ++gg_class)
//            {
//                imp_array(kk_feat,gg_class) = imp_array(kk_feat, gg_class)/permutationCount;
//                imp_array(kk_feat,gg_class)  = (imp_array(kk_feat, gg_class) - oob_error[gg_class])/oob_error[gg_class];
//            }
//        }
//    }
//}
}  //namespace rn;
#endif //RN_RF_VARIABLE_IMPORTANCE_HXX
