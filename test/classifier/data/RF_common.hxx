#ifndef RF_TEST_COMMON
#define RF_TEST_COMMON


//#include "pina_indians_var_imp.hxx"

#include "SPECTF_features.hxx"
#include "SPECTF_labels.hxx"
#include "ecoli_features.hxx"
#include "ecoli_labels.hxx"
#include "glass_features.hxx"
#include "glass_labels.hxx"
#include "ionosphere_features.hxx"
#include "ionosphere_labels.hxx"
#include "iris_features.hxx"
#include "iris_labels.hxx"
//#include "letter_recognition_features.hxx"
//#include "letter_recognition_labels.hxx"
//#include "madelon_features.hxx"
//#include "madelon_labels.hxx"
//#include "magic04_features.hxx"
//#include "magic04_labels.hxx"
//#include "page_blocks_features.hxx"
//#include "page_blocks_labels.hxx"
#include "pina_indians_diabetes_features.hxx"
#include "pina_indians_diabetes_labels.hxx"
#include "segmentation_features.hxx"
#include "segmentation_labels.hxx"
//#include "shuttle_features.hxx"
//#include "shuttle_labels.hxx"
//#include "spambase_features.hxx"
//#include "spambase_labels.hxx"
#include "wine_features.hxx"
#include "wine_labels.hxx"
#include "wpbc_features.hxx"
#include "wpbc_labels.hxx"
//#include "yeast_features.hxx"
//#include "yeast_labels.hxx"

//#include "SPECTF_resultsRFDT.hxx"
//#include "SPECTF_resultsRFm1.hxx"
//#include "SPECTF_resultsRFdef.hxx"
//#include "SPECTF_resultsRFbagging.hxx"
//
//#include "ecoli_resultsRFDT.hxx"
//#include "ecoli_resultsRFm1.hxx"
//#include "ecoli_resultsRFdef.hxx"
//#include "ecoli_resultsRFbagging.hxx"
//
//#include "glass_resultsRFDT.hxx"
//#include "glass_resultsRFm1.hxx"
//#include "glass_resultsRFdef.hxx"
//#include "glass_resultsRFbagging.hxx"
//
//#include "ionosphere_resultsRFDT.hxx"
//#include "ionosphere_resultsRFm1.hxx"
//#include "ionosphere_resultsRFdef.hxx"
//#include "ionosphere_resultsRFbagging.hxx"
//
//#include "iris_resultsRFDT.hxx"
//#include "iris_resultsRFm1.hxx"
//#include "iris_resultsRFdef.hxx"
//#include "iris_resultsRFbagging.hxx"
//
//#include "pina_indians_diabetes_resultsRFDT.hxx"
//#include "pina_indians_diabetes_resultsRFm1.hxx"
//#include "pina_indians_diabetes_resultsRFdef.hxx"
//#include "pina_indians_diabetes_resultsRFbagging.hxx"
//
//#include "segmentation_resultsRFDT.hxx"
//#include "segmentation_resultsRFm1.hxx"
//#include "segmentation_resultsRFdef.hxx"
//#include "segmentation_resultsRFbagging.hxx"
//
//#include "wine_resultsRFDT.hxx"
//#include "wine_resultsRFm1.hxx"
//#include "wine_resultsRFdef.hxx"
//#include "wine_resultsRFbagging.hxx"

//#include "wpbc_resultsRFDT.hxx"
//#include "wpbc_resultsRFm1.hxx"
//#include "wpbc_resultsRFdef.hxx"
//#include "wpbc_resultsRFbagging.hxx"
//;
enum DataSetNames
{
    SPECTF,
    ecoli,
    glass,
    ionosphere,
    iris,
    /*letter_recognition,*/
    /*madelon,*/
    /*magic04,*/
    /*page_blocks,*/
    pina_indians_diabetes,
    segmentation,
    /*shuttle,*/
    /*spambase,*/
    wine,
    //wpbc
    /*yeast*/
};
double _statistics[] = {0.347468, 0.00964554,
                        0.233222, 0.00376475,
                        0.356637, 0.00546409,
                        0.126234, 0.00313297,
                        0.0632961, 0.00275075,
                        0.313567, 0.00252546,
                        0.176406, 0.00497856,
                        0.0964343, 0.00397317};
int dataNumFeats[] =
{
                            44,
                             7,
                             9,
                            34,
                             4,
                            /*16,*/
                           /*501,*/
                            /*10,*/
                            /*10,*/
                             8,
                            19,
                             /*9,*/
                           /*57,*/
                            13,
                            //32
                             /*8*/
};

int dataNumSamps[] =
{
                              80,
                             336,
                             214,
                             351,
                             150,
                           /*20000,*/
                            /*2000,*/
                           /*19020,*/
                            /*5473,*/
                             768,
                             210,
                           /*43500,*/
                            /*4601,*/
                             178,
                             //194
                            /*1484*/
};

struct _twotuple{
    typedef int* type;
    int* classes;
    int sze;
    int* begin(){
        return classes;
    }
    int* end(){
        return classes+sze;
    }

    int size(){
        return sze;
    }
    _twotuple(int* a, int s):
        classes(a),sze(s)
    {}
};

#endif
