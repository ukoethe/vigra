#define HasHDF5
#include "../include/vigra/random_forest.hxx"
#include "../test/classifier/test_visitors.hxx"
using namespace vigra;
int main()
{
    vigra::TestVisitor testVisitor;
    RandomForest<double, RegressionTag> rf(RandomForestOptions().tree_count(1));
    Matrix<double> feat(100, 2);
    int ss = 0;
    Matrix<double> resp(100, 1);
    Matrix<double> resp2(100, 1);
    for(int ii = 0; ii < 10; ++ ii)
    {
        for(int jj = 0; jj < 10; ++jj)
        {
            feat(ss, 0 ) = ii;
            feat(ss, 1 ) = jj;
            resp(ss, 0) = ii +  jj;
            ++ss;
        }
    }
    rf.learn(feat, resp, create_visitor(testVisitor), RegressionSplit());
    rf.predictRaw(feat, resp2);
    for(int ii = 0; ii < 100; ++ii)
    {
        std::cerr << resp(ii, 0) << " " << resp2(ii, 0) << std::endl;
    }
    //shoudEqual(resp, resp2)
}
