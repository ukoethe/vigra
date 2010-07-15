
#include <vigra/timing.hxx>
USETICTOC;
#include <vigra/random_forest_deprec.hxx>
#include <vigra/random_forest.hxx>

using namespace vigra;



int main(int argc, char ** argv)
{
    typedef MultiArrayShape<2>::type Shp;
    MultiArray<2, double> features(Shp(1000, 50), 0.0);
    MultiArray<2, int>    labels(Shp(1000, 1), 0);
    ArrayVector<int> classes(2, 0);
    classes[1] = 1; 

    RandomMT19937 random(1); 
    RandomMT19937 random_old(1); 


    for(int ii = 0; ii < features.shape(0); ++ii)
    {
        labels(ii, 0) = random.uniform53() > 0.5; 
        for(int jj = 0; jj < features.shape(1); ++jj)
        {
            features(ii, jj) = random.uniform53();
        }
    }
    
    RandomForest<int>    rf_new(RandomForestOptions().tree_count(255));
    RandomForestDeprec<int> rf_old(classes.begin(), 
                                   classes.end());
    
    std::cerr << "Learning New Random Forest:" << std::endl;
    TIC;
    rf_new.learn(features, labels, rf_default(), rf_default(), rf_default(), random);
    std::cerr << TOC<< std::endl;
    std::cerr << "Learning Old Random Forest:" << std::endl;
    TIC;
    rf_old.learn(features, labels, random_old); 
    std::cerr << TOC << std::endl;
    return 0;

}
