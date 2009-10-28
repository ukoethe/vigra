#ifndef RF_DATA
#define RF_DATA

#include <vector>
#include <vigra/multi_array.hxx>
#include <iostream>
#include <string>
#include "RF_common.hxx"

#define LOAD_DATA(name)\
    {\
            int a = 0;\
            vigra::MultiArrayShape<2>::type LabelSize(1, dataNumSamps[name]);\
            vigra::MultiArrayShape<2>::type FeatureSize(dataNumFeats[name], dataNumSamps[name]);\
            vigra::MultiArrayShape<2>::type PredictSize(dataNumSamps[name], name##_size);\
            _labels.push_back(LabelType(LabelSize, name##_labels).transpose());\
            _features.push_back(FeatureType(FeatureSize, name##_features).transpose());\
            _ClassIter.push_back(_twotuple(name##_Classes, name##_size));\
            _names.push_back(#name);\
            _RFres.push_back(ManagedLabelType(PredictSize, a));\
    }


struct RF_Test_Training_Data
{
    typedef vigra::MultiArrayView<2,int, vigra::StridedArrayTag> LabelType;
    typedef vigra::MultiArrayView<2,double, vigra::StridedArrayTag> FeatureType;
    typedef vigra::MultiArray<2,double> ManagedLabelType;
    typedef int* ClassIterType;

    vigra::MultiArrayView<2, double, vigra::StridedArrayTag> _oobRange;
    std::vector<LabelType> _labels;
    std::vector<FeatureType> _features;
    std::vector<_twotuple> _ClassIter;
    std::vector<std::string> _names;
    std::vector<ManagedLabelType> _RFres;

    RF_Test_Training_Data()
    {
        _oobRange = vigra::MultiArrayView<2, double>( vigra::MultiArrayShape<2>::type(2, this->size()),
                                                    _statistics                                     ).transpose();
        _labels.reserve(this->size());
        _features.reserve(this->size());
        _ClassIter.reserve(this->size());
        _names.reserve(this->size());
        _RFres.reserve(this->size());
        LOAD_DATA(SPECTF)
        LOAD_DATA(ecoli)
        LOAD_DATA(glass)
        LOAD_DATA(ionosphere)
        LOAD_DATA(iris)
        //LOAD_DATA(letter_recognition)
        //LOAD_DATA(madelon)
        //LOAD_DATA(magic04)
        //LOAD_DATA(page_blocks)
        LOAD_DATA(pina_indians_diabetes)
        LOAD_DATA(segmentation)
        //LOAD_DATA(shuttle)
        //LOAD_DATA(spambase)
        LOAD_DATA(wine)
        //LOAD_DATA(wpbc)
        //LOAD_DATA(yeast)
/*
        std::cout << (int)_features[0].shape(0) << " " << _features[0].shape(1);
        std::cout << "Press Enter to continue...." << std::endl;
        std::cin.ignore(INT_MAX, '\n');
        std::cin.get();


        for(int ii = 0; ii < features(1).shape(0); ++ii)
        {
            for(int jj = 0; jj < features(1).shape(1); ++jj){
                std::cout << features(1).operator()(ii, jj) << ",\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl <<std::endl;
        for(int ii = 0; ii < labels(1).shape(0); ++ii)
        {
            for(int jj = 0; jj < labels(1).shape(1); ++jj){
                std::cout << labels(1).operator()(ii, jj) << ",\t";
            }
            std::cout << std::endl;
        }*/
    }

    LabelType& labels(int index)
    {
        vigra_precondition(index > 0 || index < this->size(),
            "index out of bounds");
        return _labels[index];
    }

    double oobError(int index)
    {
        vigra_precondition(index > 0 || index < this->size(),
            "index out of bounds");
        return _oobRange(index, 0);
    }

    double oobSTD(int index)
    {
        vigra_precondition(index > 0 || index < this->size(),
            "index out of bounds");
        return _oobRange(index, 1);
    }

    FeatureType& features(int index)
    {
        vigra_precondition(index > 0 || index < this->size(),
            "index out of bounds");
        return _features[index];
    }

    _twotuple ClassIter(int index)
    {
        vigra_precondition(index > 0 || index < this->size(),
            "index out of bounds");
        return _ClassIter[index];
    }

    const std::string names(int index)
    {
        vigra_precondition(index > 0 || index < this->size(),
            "index out of bounds");
        return _names[index];
    }

    ManagedLabelType& RFres(int index)
    {
        vigra_precondition(index > 0 || index < this->size(),
            "index out of bounds");
        return _RFres[index];
    }

    int size()
    {
        return 16-8;
    }
};

#endif //RF_DATA
