//
// C++ Interface: rf_ridge_split
//
// Description: 
//
//
// Author: Nico Splitthoff <splitthoff@zg00103>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef VIGRA_RANDOM_FOREST_RIDGE_SPLIT_H
#define VIGRA_RANDOM_FOREST_RIDGE_SPLIT_H
//#include "rf_sampling.hxx"
#include "../sampling.hxx"
#include "rf_split.hxx"
#include "rf_nodeproxy.hxx"
#include "../regression.hxx"

#define outm(v) std::cout << (#v) << ": " << (v) << std::endl;
#define outm2(v) std::cout << (#v) << ": " << (v) << ", ";

namespace vigra
{

/*template<>
class Node<i_RegrNode>
: public NodeBase
{
public:
    typedef NodeBase BT;


    Node(   BT::T_Container_type &   topology,
        BT::P_Container_type &   param,
            int nNumCols)
    :   BT(5+nNumCols,2+nNumCols,topology, param)
    {
        BT::typeID() = i_RegrNode;
    }

    Node(   BT::T_Container_type     &   topology,
        BT::P_Container_type     &   param,
        INT                   n             )
    :   BT(5,2,topology, param, n)
    {}

    Node( BT & node_)
    :   BT(5, 2, node_) 
    {}

    double& threshold()
    {
        return BT::parameters_begin()[1];
    }

    BT::INT& column()
    {
        return BT::column_data()[0];
    }

    template<class U, class C>
            BT::INT& next(MultiArrayView<2,U,C> const & feature)
            {
                return (feature(0, column()) < threshold())? child(0):child(1);
            }
};*/


template<class ColumnDecisionFunctor, class Tag = ClassificationTag>
class RidgeSplit: public SplitBase<Tag>
{
  public:


    typedef SplitBase<Tag> SB;

    ArrayVector<Int32>          splitColumns;
    ColumnDecisionFunctor       bgfunc;

    double                      region_gini_;
    ArrayVector<double>         min_gini_;
    ArrayVector<std::ptrdiff_t> min_indices_;
    ArrayVector<double>         min_thresholds_;

    int                         bestSplitIndex;
    
    //dns
    bool            m_bDoScalingInTraining;
    bool            m_bDoBestLambdaBasedOnGini;
    
    RidgeSplit()
    :m_bDoScalingInTraining(true),
    m_bDoBestLambdaBasedOnGini(true)
    {
    }

    double minGini() const
    {
        return min_gini_[bestSplitIndex];
    }
    
    int bestSplitColumn() const
    {
        return splitColumns[bestSplitIndex];
    }
    
    bool& doScalingInTraining()
    { return m_bDoScalingInTraining; }

    bool& doBestLambdaBasedOnGini()
    { return m_bDoBestLambdaBasedOnGini; }

    template<class T>
            void set_external_parameters(ProblemSpec<T> const & in)
    {
        SB::set_external_parameters(in);        
        bgfunc.set_external_parameters(in);
        int featureCount_ = in.column_count_;
        splitColumns.resize(featureCount_);
        for(int k=0; k<featureCount_; ++k)
            splitColumns[k] = k;
        min_gini_.resize(featureCount_);
        min_indices_.resize(featureCount_);
        min_thresholds_.resize(featureCount_);
    }

    
    template<class T, class C, class T2, class C2, class Region, class Random>
    int findBestSplit(MultiArrayView<2, T, C> features,
                      MultiArrayView<2, T2, C2>  multiClassLabels,
                      Region & region,
                      ArrayVector<Region>& childRegions,
                      Random & randint)
    {

    //std::cerr << "Split called" << std::endl;
    typedef typename Region::IndexIterator IndexIterator;
    typedef typename MultiArrayView <2, T, C>::difference_type fShape;
    typedef typename MultiArrayView <2, T2, C2>::difference_type lShape;
    typedef typename MultiArrayView <2, double>::difference_type dShape;
        
        // calculate things that haven't been calculated yet. 
//    std::cout << "start" << std::endl;
        // calculate things that haven't been calculated yet.
        detail::Correction<Tag>::exec(region, multiClassLabels);

        //mexPrintf("Starting Counts: %d %d\n", (int)region.classCounts()[0], (int)region.classCounts()[1] );

        // Is the region pure already?
#if 0
        region_gini_ = GiniCriterion::impurity(region.classCounts(),
                region.size());
#else
               region_gini_ = bgfunc.loss_of_region(multiClassLabels,
                                                     region.begin(),
                                                     region.end(),
                                                     region.classCounts());
#endif
        if(region_gini_ == SB::ext_param_.precision_|| region.size() < SB::ext_param_.actual_mtry_ || region.oob_size() < 2)
            return  SB::makeTerminalNode(features, multiClassLabels, region, randint);

        // select columns  to be tried.
    for(int ii = 0; ii < SB::ext_param_.actual_mtry_; ++ii)
        std::swap(splitColumns[ii], 
            splitColumns[ii+ randint(features.shape(1) - ii)]);

    //do implicit binary case
    MultiArray<2, T2> labels(lShape(multiClassLabels.shape(0),1));
      //number of classes should be >1, otherwise makeTerminalNode would have been called
      int nNumClasses=0;
      for(int n=0; n<(int)region.classCounts().size(); n++)
        nNumClasses+=((region.classCounts()[n]>0) ? 1:0);

      //convert to binary case
      if(nNumClasses>2)
      {
        int nMaxClass=0;
        int nMaxClassCounts=0;
        for(int n=0; n<(int)region.classCounts().size(); n++)
        {
          //this should occur in any case:
          //we had more than two non-zero classes in order to get here
          if(region.classCounts()[n]>nMaxClassCounts)
          {
        nMaxClassCounts=region.classCounts()[n];
        nMaxClass=n;
          }
        }
        
        //create binary labels
        for(int n=0; n<multiClassLabels.shape(0); n++)
          labels(n,0)=((multiClassLabels(n,0)==nMaxClass) ? 1:0);
      }
      else
        labels=multiClassLabels;

    //_do implicit binary case
    
    //uncomment this for some debugging
/*  int nNumCases=features.shape(0);

    typedef typename MultiArrayView <2, int>::difference_type nShape;
    MultiArray<2, int> elementCounterArray(nShape(nNumCases,1),(int)0);
    int nUniqueElements=0;
    for(int n=0; n<region.size(); n++)
        elementCounterArray[region[n]]++;
    
    for(int n=0; n<nNumCases; n++)
        nUniqueElements+=((elementCounterArray[n]>0) ? 1:0);
    
    outm(nUniqueElements);
    nUniqueElements=0;
    MultiArray<2, int> elementCounterArray_oob(nShape(nNumCases,1),(int)0);
    for(int n=0; n<region.oob_size(); n++)
        elementCounterArray_oob[region.oob_begin()[n]]++;
    for(int n=0; n<nNumCases; n++)
        nUniqueElements+=((elementCounterArray_oob[n]>0) ? 1:0);
    outm(nUniqueElements);
    
    int notUniqueElements=0;
    for(int n=0; n<nNumCases; n++)
        notUniqueElements+=(((elementCounterArray_oob[n]>0) && (elementCounterArray[n]>0)) ? 1:0);
    outm(notUniqueElements);*/
    
    //outm(SB::ext_param_.actual_mtry_);
    
    
//select submatrix of features for regression calculation
    //MultiArrayView<2, T, C> cVector( fShape(region.size()), 1 );
    MultiArray<2, T> xtrain(fShape(region.size(),SB::ext_param_.actual_mtry_));
    //we only want -1 and 1 for this
    MultiArray<2, double> regrLabels(dShape(region.size(),1));

    //copy data into a vigra data structure and centre and scale while doing so
    MultiArray<2, double> meanMatrix(dShape(SB::ext_param_.actual_mtry_,1));
    MultiArray<2, double> stdMatrix(dShape(SB::ext_param_.actual_mtry_,1));
    for(int m=0; m<SB::ext_param_.actual_mtry_; m++)
    {
        MultiArrayView<2,T,C> cVector( columnVector(features, splitColumns[m]) );
        
        //centre and scale the data
        double dCurrFeatureColumnMean=0.0;
        double dCurrFeatureColumnStd=1.0; //default value
        
        //calc mean on bootstrap data
        for(int n=0; n<region.size(); n++)
          dCurrFeatureColumnMean+=cVector[region[n]];
        dCurrFeatureColumnMean/=region.size();


        //calc scaling
        if(m_bDoScalingInTraining)
        {
          for(int n=0; n<region.size(); n++)
          {
              dCurrFeatureColumnStd+=
            (cVector[region[n]]-dCurrFeatureColumnMean)*(cVector[region[n]]-dCurrFeatureColumnMean);
          }
          //unbiased std estimator:
          dCurrFeatureColumnStd=sqrt(dCurrFeatureColumnStd/(region.size()-1));
        }

        //dCurrFeatureColumnStd is still 1.0 if we didn't want scaling
        stdMatrix(m,0)=dCurrFeatureColumnStd;
        
        meanMatrix(m,0)=dCurrFeatureColumnMean;
        
        //get feature matrix, i.e. A (note that weighting is done automatically
        //since rows can occur multiple times -> bagging)
        for(int n=0; n<region.size(); n++)
            xtrain(n,m)=(cVector[region[n]]-dCurrFeatureColumnMean)/dCurrFeatureColumnStd;
    }

    //return  SB::makeTerminalNode(features, multiClassLabels, region, randint);

//    std::cout << "middle" << std::endl;
    //get label vector (i.e. b)
    for(int n=0; n<region.size(); n++)
    {
        //we checked for/built binary case further up.
        //class labels should thus be either 0 or 1
        //-> convert to -1 and 1 for regression
        regrLabels(n,0)=((labels[region[n]]==0) ? -1:1);
    }

    MultiArray<2, double> dLambdas(dShape(11,1));
    int nCounter=0;
    for(int nLambda=-5; nLambda<=5; nLambda++)
        dLambdas[nCounter++]=pow(10.0,nLambda);
    //destination vector for regression coefficients; use same type as for xtrain
    MultiArray<2, double> regrCoef(dShape(SB::ext_param_.actual_mtry_,11));
    ridgeRegressionSeries(xtrain,regrLabels,regrCoef,dLambdas);
    
    double dMaxRidgeSum=NumericTraits<double>::min();
    double dCurrRidgeSum;
    int nMaxRidgeSumAtLambdaInd=0;

    for(int nLambdaInd=0; nLambdaInd<11; nLambdaInd++)
    {
        //just sum up the correct answers
        //(correct means >=intercept for class 1, <intercept for class 0)
        //(intercept=0 or intercept=threshold based on gini)
        dCurrRidgeSum=0.0;
        
        //assemble projection vector
        MultiArray<2, double> dDistanceFromHyperplane(dShape(features.shape(0),1));
        
        for(int n=0; n<region.oob_size(); n++)
        {
          dDistanceFromHyperplane(region.oob_begin()[n],0)=0.0;
          for (int m=0; m<SB::ext_param_.actual_mtry_; m++)
          {
            dDistanceFromHyperplane(region.oob_begin()[n],0)+=
              features(region.oob_begin()[n],splitColumns[m])*regrCoef(m,nLambdaInd);
          }
        }

        double dCurrIntercept=0.0;
        if(m_bDoBestLambdaBasedOnGini)
        {
          //calculate gini index
          bgfunc(dDistanceFromHyperplane,
              labels, 
              region.oob_begin(), region.oob_end(), 
              region.classCounts());
          dCurrIntercept=bgfunc.min_threshold_;
        }
        else
        {
          for (int m=0; m<SB::ext_param_.actual_mtry_; m++)
            dCurrIntercept+=meanMatrix(m,0)*regrCoef(m,nLambdaInd);
        }
        
        for(int n=0; n<region.oob_size(); n++)
        {
            //check what lambda performs best on oob data
            int nClassPrediction=((dDistanceFromHyperplane(region.oob_begin()[n],0) >=dCurrIntercept) ? 1:0);
            dCurrRidgeSum+=((nClassPrediction == labels(region.oob_begin()[n],0)) ? 1:0);
        }
        if(dCurrRidgeSum>dMaxRidgeSum)
        {
            dMaxRidgeSum=dCurrRidgeSum;
            nMaxRidgeSumAtLambdaInd=nLambdaInd;
        }
    }

//    std::cout << "middle2" << std::endl;
        //create a Node for output
        Node<i_HyperplaneNode>   node(SB::ext_param_.actual_mtry_, SB::t_data, SB::p_data);

    //normalise coeffs
        //data was scaled (by 1.0 or by std) -> take into account
        MultiArray<2, double> dCoeffVector(dShape(SB::ext_param_.actual_mtry_,1));
        for(int n=0; n<SB::ext_param_.actual_mtry_; n++)
          dCoeffVector(n,0)=regrCoef(n,nMaxRidgeSumAtLambdaInd)*stdMatrix(n,0);
        
        //calc norm
        //double dVnorm=columnVector(regrCoef,nMaxRidgeSumAtLambdaInd).norm();
        double dVnorm=columnVector(dCoeffVector,0).norm();

        for(int n=0; n<SB::ext_param_.actual_mtry_; n++)
            node.weights()[n]=dCoeffVector(n,0)/dVnorm;
#if 0
        mexPrintf("Norm: %f\n", (float)dVnorm);
        mexPrintf("W: ");
        for(int n=0; n<SB::ext_param_.actual_mtry_; n++)
            mexPrintf("%.2f ", (float)node.weights()[n]);
        mexPrintf("\n");
#endif

    //_normalise coeffs
    
    //save the columns
        node.column_data()[0]=SB::ext_param_.actual_mtry_;
        for(int n=0; n<SB::ext_param_.actual_mtry_; n++){
            node.column_data()[n+1]=splitColumns[n];
            //mexPrintf("Col: %d\n", splitColumns[n]);
        }

    //assemble projection vector
        //careful here: "region" is a pointer to indices...
        //all the indices in "region" need to have valid data
        //convert from "region" space to original "feature" space
        MultiArray<2, double> dDistanceFromHyperplane(dShape(features.shape(0),1));
        
        for(int n=0; n<region.size(); n++)
        {
            dDistanceFromHyperplane(region[n],0)=0.0;
            for (int m=0; m<SB::ext_param_.actual_mtry_; m++)
            {
              Int32 Cc = splitColumns[m];
              dDistanceFromHyperplane(region[n],0)+=
               features(region[n],Cc)*node.weights()[m];
            }
            //mexPrintf("Dist %d: %.2f\n", region[n], dDistanceFromHyperplane(region[n],0));
        }
        for(int n=0; n<region.oob_size(); n++)
        {
            dDistanceFromHyperplane(region.oob_begin()[n],0)=0.0;
            for (int m=0; m<SB::ext_param_.actual_mtry_; m++)
            {
              Int32 Cc = splitColumns[m];
              dDistanceFromHyperplane(region.oob_begin()[n],0)+=
            features(region.oob_begin()[n],Cc)*node.weights()[m];
            }
        }
        
    //calculate gini index
        bgfunc(dDistanceFromHyperplane,
            labels,
               region.begin(), region.end(),
               region.classCounts());

        // did not find any suitable split
    if(closeAtTolerance(bgfunc.min_gini_, region_gini_))
        return  SB::makeTerminalNode(features, multiClassLabels, region, randint);

    //take gini threshold here due to scaling, normalisation, etc. of the coefficients
    node.intercept()    = bgfunc.min_threshold_;

    //mexPrintf("Offset: %f\n\n", node.intercept() );
    SB::node_ = node;


    childRegions[0].classCounts() = bgfunc.bestCurrentCounts[0];
    childRegions[1].classCounts() = bgfunc.bestCurrentCounts[1];
    childRegions[0].classCountsIsValid = true;
    childRegions[1].classCountsIsValid = true;

#if 1

    // Save the ranges of the child stack entries.
    childRegions[0].setRange(   region.begin()  , region.begin() + bgfunc.min_index_   );
    childRegions[0].rule = region.rule;
    childRegions[0].rule.push_back(std::make_pair(1, 1.0));
    childRegions[1].setRange(   region.begin() + bgfunc.min_index_       , region.end()    );
    childRegions[1].rule = region.rule;
    childRegions[1].rule.push_back(std::make_pair(1, 1.0));

#else
    SortSamplesByDimensions<MultiArrayView<2, T, C> >
        sorter(dDistanceFromHyperplane, 0, bgfunc.min_threshold_);

    IndexIterator bestSplit =
        std::partition(region.begin(), region.end(), sorter);

    childRegions[0].setRange(   region.begin()  , bestSplit       );
    childRegions[0].rule = region.rule;
    childRegions[0].rule.push_back(std::make_pair(1, 1.0));
    childRegions[1].setRange(   bestSplit       , region.end()    );
    childRegions[1].rule = region.rule;
    childRegions[1].rule.push_back(std::make_pair(1, 1.0));
#endif

#if 0
    mexPrintf("Counts: %d %d\n", (int)childRegions[0].classCounts()[0], (int)childRegions[0].classCounts()[1] );
    mexPrintf("Counts: %d %d\n", (int)childRegions[1].classCounts()[0], (int)childRegions[1].classCounts()[1] );

    mexPrintf("Example R0, %d: %f %f\n", childRegions[0][0], (float)features(childRegions[0][0],0), (float)features(childRegions[0][0],1) );
    mexPrintf("Label R0: %f\n", (float) multiClassLabels(childRegions[0][0],0) );
#endif

    //adjust oob ranges
//    std::cout << "adjust oob" << std::endl;
    //sort the oobs
    //Region regionCopy = region;
      std::sort(region.oob_begin(), region.oob_end(),
            SortSamplesByDimensions< MultiArray<2, double> > (dDistanceFromHyperplane, 0));
            
      //find split index
      int nOOBindx;
      for(nOOBindx=0; nOOBindx<region.oob_size(); nOOBindx++)
      {
        if(dDistanceFromHyperplane(region.oob_begin()[nOOBindx],0)>=node.intercept())
          break;
      }

      childRegions[0].set_oob_range(   region.oob_begin()  , region.oob_begin() + nOOBindx   );
      childRegions[1].set_oob_range(   region.oob_begin() + nOOBindx , region.oob_end() );

      //    std::cout << "end" << std::endl;
//    outm2(region.oob_begin());outm2(nOOBindx);outm(region.oob_begin() + nOOBindx);
    //_adjust oob ranges

    return i_HyperplaneNode;
    }
};

/** Standard ridge regression split
 */
typedef RidgeSplit<BestGiniOfColumn<GiniCriterion> >  GiniRidgeSplit;


} //namespace vigra
#endif // VIGRA_RANDOM_FOREST_RIDGE_SPLIT_H
