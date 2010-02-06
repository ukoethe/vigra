/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
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


#ifndef VIGRA_RANDOM_FOREST_HXX
#define VIGRA_RANDOM_FOREST_HXX

#include <algorithm>
#include <map>
#include <numeric>
#include <iostream>
#include "vigra/mathutil.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/sized_int.hxx"
#include "vigra/matrix.hxx"
#include "vigra/random.hxx"
#include "vigra/functorexpression.hxx"

#define RandomForest RandomForestDeprec
#define DecisionTree DecisionTreeeDeprec

namespace vigra
{

namespace detail
{

template<class DataMatrix>
class RandomForestFeatureSorter
{
    DataMatrix const & data_;
    MultiArrayIndex sortColumn_;

  public:

    RandomForestFeatureSorter(DataMatrix const & data, MultiArrayIndex sortColumn)
    : data_(data),
      sortColumn_(sortColumn)
    {}

    void setColumn(MultiArrayIndex sortColumn)
    {
        sortColumn_ = sortColumn;
    }

    bool operator()(MultiArrayIndex l, MultiArrayIndex r) const
    {
        return data_(l, sortColumn_) < data_(r, sortColumn_);
    }
};

template<class LabelArray>
class RandomForestLabelSorter
{
    LabelArray const & labels_;

  public:

    RandomForestLabelSorter(LabelArray const & labels)
    : labels_(labels)
    {}

    bool operator()(MultiArrayIndex l, MultiArrayIndex r) const
    {
        return labels_[l] < labels_[r];
    }
};

template <class CountArray>
class RandomForestClassCounter
{
    ArrayVector<int> const & labels_;
    CountArray & counts_;

  public:

    RandomForestClassCounter(ArrayVector<int> const & labels, CountArray & counts)
    : labels_(labels),
      counts_(counts)
    {
        reset();
    }

    void reset()
    {
        counts_.init(0);
    }

    void operator()(MultiArrayIndex l) const
    {
        ++counts_[labels_[l]];
    }
};

struct DecisionTreeCountNonzeroFunctor
{
    double operator()(double old, double other) const
    {
        if(other != 0.0)
            ++old;
        return old;
    }
};

struct DecisionTreeNode
{
    DecisionTreeNode(int t, MultiArrayIndex bestColumn)
    : thresholdIndex(t), splitColumn(bestColumn)
    {}

    int children[2];
    int thresholdIndex;
    Int32 splitColumn;
};

template <class INT>
struct DecisionTreeNodeProxy
{
    DecisionTreeNodeProxy(ArrayVector<INT> const & tree, INT n)
    : node(const_cast<ArrayVector<INT> &>(tree).begin()+n)
    {}

    INT & child(INT l) const
    {
        return node[l];
    }

    INT & decisionWeightsIndex() const
    {
        return node[2];
    }

    typename ArrayVector<INT>::iterator decisionColumns() const
    {
        return node+3;
    }

    mutable typename ArrayVector<INT>::iterator node;
};

struct DecisionTreeAxisSplitFunctor
{
    ArrayVector<Int32> splitColumns;
    ArrayVector<double> classCounts, currentCounts[2], bestCounts[2], classWeights;
    double threshold;
    double totalCounts[2], bestTotalCounts[2];
    int mtry, classCount, bestSplitColumn;
    bool pure[2], isWeighted;

    void init(int mtry, int cols, int classCount, ArrayVector<double> const & weights)
    {
        this->mtry = mtry;
        splitColumns.resize(cols);
        for(int k=0; k<cols; ++k)
            splitColumns[k] = k;

        this->classCount = classCount;
        classCounts.resize(classCount);
        currentCounts[0].resize(classCount);
        currentCounts[1].resize(classCount);
        bestCounts[0].resize(classCount);
        bestCounts[1].resize(classCount);

        isWeighted = weights.size() > 0;
        if(isWeighted)
            classWeights = weights;
        else
            classWeights.resize(classCount, 1.0);
    }

    bool isPure(int k) const
    {
        return pure[k];
    }

    unsigned int totalCount(int k) const
    {
        return (unsigned int)bestTotalCounts[k];
    }

    int sizeofNode() const { return 4; }

    int writeSplitParameters(ArrayVector<Int32> & tree,
                                ArrayVector<double> &terminalWeights)
    {
        int currentWeightIndex = terminalWeights.size();
        terminalWeights.push_back(threshold);

        int currentNodeIndex = tree.size();
        tree.push_back(-1);  // left child
        tree.push_back(-1);  // right child
        tree.push_back(currentWeightIndex);
        tree.push_back(bestSplitColumn);

        return currentNodeIndex;
    }

    void writeWeights(int l, ArrayVector<double> &terminalWeights)
    {
        for(int k=0; k<classCount; ++k)
            terminalWeights.push_back(isWeighted
                                           ? bestCounts[l][k]
                                           : bestCounts[l][k] / totalCount(l));
    }

    template <class U, class C, class AxesIterator, class WeightIterator>
    bool decideAtNode(MultiArrayView<2, U, C> const & features,
                      AxesIterator a, WeightIterator w) const
    {
        return (features(0, *a) < *w);
    }

    template <class U, class C, class IndexIterator, class Random>
    IndexIterator findBestSplit(MultiArrayView<2, U, C> const & features,
                                ArrayVector<int> const & labels,
                                IndexIterator indices, int exampleCount,
                                Random & randint);

};


template <class U, class C, class IndexIterator, class Random>
IndexIterator
DecisionTreeAxisSplitFunctor::findBestSplit(MultiArrayView<2, U, C> const & features,
                                            ArrayVector<int> const & labels,
                                            IndexIterator indices, int exampleCount,
                                            Random & randint)
{
    // select columns to be tried for split
    for(int k=0; k<mtry; ++k)
        std::swap(splitColumns[k], splitColumns[k+randint(columnCount(features)-k)]);

    RandomForestFeatureSorter<MultiArrayView<2, U, C> > sorter(features, 0);
    RandomForestClassCounter<ArrayVector<double> > counter(labels, classCounts);
    std::for_each(indices, indices+exampleCount, counter);

	// find the best gini index
    double minGini = NumericTraits<double>::max();
    IndexIterator bestSplit;
    for(int k=0; k<mtry; ++k)
    {
        sorter.setColumn(splitColumns[k]);
        std::sort(indices, indices+exampleCount, sorter);

        currentCounts[0].init(0);
        std::transform(classCounts.begin(), classCounts.end(), classWeights.begin(),
                       currentCounts[1].begin(), std::multiplies<double>());
        totalCounts[0] = 0;
        totalCounts[1] = std::accumulate(currentCounts[1].begin(), currentCounts[1].end(), 0.0);
        for(int m = 0; m < exampleCount-1; ++m)
        {
            int label = labels[indices[m]];
            double w = classWeights[label];
            currentCounts[0][label] += w;
            totalCounts[0] += w;
            currentCounts[1][label] -= w;
			totalCounts[1] -= w;

			if (m < exampleCount-2 &&
                features(indices[m], splitColumns[k]) == features(indices[m+1], splitColumns[k]))
				continue ;

            double gini = 0.0;
            if(classCount == 2)
            {
                gini = currentCounts[0][0]*currentCounts[0][1] / totalCounts[0] +
                       currentCounts[1][0]*currentCounts[1][1] / totalCounts[1];
            }
            else
            {
                for(int l=0; l<classCount; ++l)
                    gini += currentCounts[0][l]*(1.0 - currentCounts[0][l] / totalCounts[0]) +
                            currentCounts[1][l]*(1.0 - currentCounts[1][l] / totalCounts[1]);
            }
            if(gini < minGini)
            {
                minGini = gini;
                bestSplit = indices+m;
                bestSplitColumn = splitColumns[k];
                bestCounts[0] = currentCounts[0];
                bestCounts[1] = currentCounts[1];
            }
        }



    }
        //std::cerr << minGini << " " << bestSplitColumn << std::endl;
	// split using the best feature
    sorter.setColumn(bestSplitColumn);
    std::sort(indices, indices+exampleCount, sorter);

    for(int k=0; k<2; ++k)
    {
        bestTotalCounts[k] = std::accumulate(bestCounts[k].begin(), bestCounts[k].end(), 0.0);
    }

    threshold = (features(bestSplit[0], bestSplitColumn) + features(bestSplit[1], bestSplitColumn)) / 2.0;
    ++bestSplit;

    counter.reset();
    std::for_each(indices, bestSplit, counter);
    pure[0] = 1.0 == std::accumulate(classCounts.begin(), classCounts.end(), 0.0, DecisionTreeCountNonzeroFunctor());
    counter.reset();
    std::for_each(bestSplit, indices+exampleCount, counter);
    pure[1] = 1.0 == std::accumulate(classCounts.begin(), classCounts.end(), 0.0, DecisionTreeCountNonzeroFunctor());

    return bestSplit;
}

enum  { DecisionTreeNoParent = -1 };

template <class Iterator>
struct DecisionTreeStackEntry
{
    DecisionTreeStackEntry(Iterator i, int c,
                           int lp = DecisionTreeNoParent, int rp = DecisionTreeNoParent)
    : indices(i), exampleCount(c),
      leftParent(lp), rightParent(rp)
    {}

    Iterator indices;
    int exampleCount, leftParent, rightParent;
};

class DecisionTree
{
  public:
    typedef Int32 TreeInt;
    ArrayVector<TreeInt>  tree_;
    ArrayVector<double> terminalWeights_;
    unsigned int classCount_;
    DecisionTreeAxisSplitFunctor split;

  public:


    DecisionTree(unsigned int classCount)
    : classCount_(classCount)
    {}

    void reset(unsigned int classCount = 0)
    {
        if(classCount)
            classCount_ = classCount;
        tree_.clear();
        terminalWeights_.clear();
    }

    template <class U, class C, class Iterator, class Options, class Random>
    void learn(MultiArrayView<2, U, C> const & features,
               ArrayVector<int> const & labels,
               Iterator indices, int exampleCount,
               Options const & options,
               Random & randint);

    template <class U, class C>
    ArrayVector<double>::const_iterator
    predict(MultiArrayView<2, U, C> const & features) const
    {
        int nodeindex = 0;
        for(;;)
        {
            DecisionTreeNodeProxy<TreeInt> node(tree_, nodeindex);
            nodeindex = split.decideAtNode(features, node.decisionColumns(),
                                       terminalWeights_.begin() + node.decisionWeightsIndex())
                                ? node.child(0)
                                : node.child(1);
            if(nodeindex <= 0)
                return terminalWeights_.begin() + (-nodeindex);
        }
    }

    template <class U, class C>
    int
    predictLabel(MultiArrayView<2, U, C> const & features) const
    {
        ArrayVector<double>::const_iterator weights = predict(features);
        return argMax(weights, weights+classCount_) - weights;
    }

    template <class U, class C>
    int
    leafID(MultiArrayView<2, U, C> const & features) const
    {
        int nodeindex = 0;
        for(;;)
        {
            DecisionTreeNodeProxy<TreeInt> node(tree_, nodeindex);
            nodeindex = split.decideAtNode(features, node.decisionColumns(),
                                       terminalWeights_.begin() + node.decisionWeightsIndex())
                                ? node.child(0)
                                : node.child(1);
            if(nodeindex <= 0)
                return -nodeindex;
        }
    }

    void depth(int & maxDep, int & interiorCount, int & leafCount, int k = 0, int d = 1) const
    {
        DecisionTreeNodeProxy<TreeInt> node(tree_, k);
        ++interiorCount;
        ++d;
        for(int l=0; l<2; ++l)
        {
            int child = node.child(l);
            if(child > 0)
                depth(maxDep, interiorCount, leafCount, child, d);
            else
            {
                ++leafCount;
                if(maxDep < d)
                    maxDep = d;
            }
        }
    }

    void printStatistics(std::ostream & o) const
    {
        int maxDep = 0, interiorCount = 0, leafCount = 0;
        depth(maxDep, interiorCount, leafCount);

        o << "interior nodes: " << interiorCount <<
             ", terminal nodes: " << leafCount <<
             ", depth: " << maxDep << "\n";
    }

    void print(std::ostream & o, int k = 0, std::string s = "") const
    {
        DecisionTreeNodeProxy<TreeInt> node(tree_, k);
        o << s << (*node.decisionColumns()) << " " << terminalWeights_[node.decisionWeightsIndex()] << "\n";

        for(int l=0; l<2; ++l)
        {
            int child = node.child(l);
            if(child <= 0)
                o << s << " weights " << terminalWeights_[-child] << " "
                                      << terminalWeights_[-child+1] << "\n";
            else
                print(o, child, s+" ");
        }
    }
};


template <class U, class C, class Iterator, class Options, class Random>
void DecisionTree::learn(MultiArrayView<2, U, C> const & features,
                          ArrayVector<int> const & labels,
                          Iterator indices, int exampleCount,
                          Options const & options,
                          Random & randint)
{
    ArrayVector<double> const & classLoss = options.class_weights;

    vigra_precondition(classLoss.size() == 0 || classLoss.size() == classCount_,
        "DecisionTree2::learn(): class weights array has wrong size.");

    reset();

    unsigned int mtry = options.mtry;
    MultiArrayIndex cols = columnCount(features);

    split.init(mtry, cols, classCount_, classLoss);

    typedef DecisionTreeStackEntry<Iterator> Entry;
    ArrayVector<Entry> stack;
    stack.push_back(Entry(indices, exampleCount));

    while(!stack.empty())
    {
//        std::cerr << "*";
        indices = stack.back().indices;
        exampleCount = stack.back().exampleCount;
        int leftParent  = stack.back().leftParent,
            rightParent = stack.back().rightParent;

        stack.pop_back();

        Iterator bestSplit = split.findBestSplit(features, labels, indices, exampleCount, randint);


        int currentNode = split.writeSplitParameters(tree_, terminalWeights_);

        if(leftParent != DecisionTreeNoParent)
            DecisionTreeNodeProxy<TreeInt>(tree_, leftParent).child(0) = currentNode;
        if(rightParent != DecisionTreeNoParent)
            DecisionTreeNodeProxy<TreeInt>(tree_, rightParent).child(1) = currentNode;
        leftParent = currentNode;
        rightParent = DecisionTreeNoParent;

        for(int l=0; l<2; ++l)
        {

            if(!split.isPure(l) && split.totalCount(l) >= options.min_split_node_size)
            {
                // sample is still large enough and not yet perfectly separated => split
                stack.push_back(Entry(indices, split.totalCount(l), leftParent, rightParent));
            }
            else
            {
                DecisionTreeNodeProxy<TreeInt>(tree_, currentNode).child(l) = -(TreeInt)terminalWeights_.size();

                split.writeWeights(l, terminalWeights_);
            }
            std::swap(leftParent, rightParent);
            indices = bestSplit;
        }
    }
//    std::cerr << "\n";
}

} // namespace detail

class RandomForestOptions
{
  public:
        /** Initialize all options with default values.
        */
    RandomForestOptions()
    : training_set_proportion(1.0),
      mtry(0),
      min_split_node_size(1),
      training_set_size(0),
      sample_with_replacement(true),
      sample_classes_individually(false),
      treeCount(255)
    {}

        /** Number of features considered in each node.

            If \a n is 0 (the default), the number of features tried in every node
            is determined by the square root of the total number of features.
            According to Breiman, this quantity should slways be optimized by means
            of the out-of-bag error.<br>
            Default: 0 (use <tt>sqrt(columnCount(featureMatrix))</tt>)
        */
    RandomForestOptions & featuresPerNode(unsigned int n)
    {
        mtry = n;
        return *this;
    }

        /** How to sample the subset of the training data for each tree.

            Each tree is only trained with a subset of the entire training data.
            If \a r is <tt>true</tt>, this subset is sampled from the entire training set with
            replacement.<br>
            Default: true (use sampling with replacement)</tt>)
        */
    RandomForestOptions & sampleWithReplacement(bool r)
    {
        sample_with_replacement = r;
        return *this;
    }

    RandomForestOptions & setTreeCount(unsigned int cnt)
    {
        treeCount = cnt;
        return *this;
    }
        /** Proportion of training examples used for each tree.

            If \a p is 1.0 (the default), and samples are drawn with replacement,
            the training set of each tree will contain as many examples as the entire
            training set, but some are drawn multiply and others not at all. On average,
            each tree is actually trained on about 65% of the examples in the full
            training set. Changing the proportion makes mainly sense when
            sampleWithReplacement() is set to <tt>false</tt>. trainingSetSizeProportional() gets
            overridden by trainingSetSizeAbsolute().<br>
            Default: 1.0
        */
    RandomForestOptions & trainingSetSizeProportional(double p)
    {
        vigra_precondition(p >= 0.0 && p <= 1.0,
            "RandomForestOptions::trainingSetSizeProportional(): proportion must be in [0, 1].");
        if(training_set_size == 0) // otherwise, absolute size gets priority
            training_set_proportion = p;
        return *this;
    }

        /** Size of the training set for each tree.

            If this option is set, it overrides the proportion set by
            trainingSetSizeProportional(). When classes are sampled individually,
            the number of examples is divided by the number of classes (rounded upwards)
            to determine the number of examples drawn from every class.<br>
            Default: <tt>0</tt> (determine size by proportion)
        */
    RandomForestOptions & trainingSetSizeAbsolute(unsigned int s)
    {
        training_set_size = s;
        if(s > 0)
            training_set_proportion = 0.0;
        return *this;
    }

        /** Are the classes sampled individually?

            If \a s is <tt>false</tt> (the default), the training set for each tree is sampled
            without considering class labels. Otherwise, samples are drawn from each
            class independently. The latter is especially useful in connection
            with the specification of an absolute training set size: then, the same number of
            examples is drawn from every class. This can be used as a counter-measure when the
            classes are very unbalanced in size.<br>
            Default: <tt>false</tt>
        */
    RandomForestOptions & sampleClassesIndividually(bool s)
    {
        sample_classes_individually = s;
        return *this;
    }

        /** Number of examples required for a node to be split.

            When the number of examples in a node is below this number, the node is not
            split even if class separation is not yet perfect. Instead, the node returns
            the proportion of each class (among the remaining examples) during the
            prediction phase.<br>
            Default: 1 (complete growing)
        */
    RandomForestOptions & minSplitNodeSize(unsigned int n)
    {
        if(n == 0)
            n = 1;
        min_split_node_size = n;
        return *this;
    }

        /** Use a weighted random forest.

            This is usually used to penalize the errors for the minority class.
            Weights must be convertible to <tt>double</tt>, and the array of weights
            must contain as many entries as there are classes.<br>
            Default: do not use weights
        */
    template <class WeightIterator>
    RandomForestOptions & weights(WeightIterator weights, unsigned int classCount)
    {
        class_weights.clear();
        if(weights != 0)
            class_weights.insert(weights, classCount);
        return *this;
    }

    RandomForestOptions & oobData(MultiArrayView<2, UInt8>& data)
    {
        oob_data =data;
        return *this;
    }

    MultiArrayView<2, UInt8> oob_data;
    ArrayVector<double> class_weights;
    double training_set_proportion;
    unsigned int mtry, min_split_node_size, training_set_size;
    bool sample_with_replacement, sample_classes_individually;
    unsigned int treeCount;
};

/*****************************************************************/
/*                                                               */
/*                          RandomForest                         */
/*                                                               */
/*****************************************************************/

template <class ClassLabelType>
class RandomForest
{
  public:
    ArrayVector<ClassLabelType> classes_;
    ArrayVector<detail::DecisionTree> trees_;
    MultiArrayIndex columnCount_;
    RandomForestOptions options_;

  public:

	//First two constructors are straight forward.
	//they take either the iterators to an Array of Classlabels or the values
    template<class ClassLabelIterator>
    RandomForest(ClassLabelIterator cl, ClassLabelIterator cend,
                  unsigned int treeCount = 255,
                  RandomForestOptions const & options = RandomForestOptions())
    : classes_(cl, cend),
      trees_(treeCount, detail::DecisionTree(classes_.size())),
      columnCount_(0),
      options_(options)
    {
        vigra_precondition(options.training_set_proportion == 0.0 ||
                           options.training_set_size == 0,
            "RandomForestOptions: absolute and proprtional training set sizes "
            "cannot be specified at the same time.");
        vigra_precondition(classes_.size() > 1,
            "RandomForestOptions::weights(): need at least two classes.");
        vigra_precondition(options.class_weights.size() == 0 || options.class_weights.size() == classes_.size(),
            "RandomForestOptions::weights(): wrong number of classes.");
    }

    RandomForest(ClassLabelType const & c1, ClassLabelType const & c2,
                  unsigned int treeCount = 255,
                  RandomForestOptions const & options = RandomForestOptions())
    : classes_(2),
      trees_(treeCount, detail::DecisionTree(2)),
      columnCount_(0),
      options_(options)
    {
        vigra_precondition(options.class_weights.size() == 0 || options.class_weights.size() == 2,
            "RandomForestOptions::weights(): wrong number of classes.");
        classes_[0] = c1;
        classes_[1] = c2;
    }
    //This is esp. For the CrosValidator Class
    template<class ClassLabelIterator>
    RandomForest(ClassLabelIterator cl, ClassLabelIterator cend,
                  RandomForestOptions const & options )
    : classes_(cl, cend),
      trees_(options.treeCount , detail::DecisionTree(classes_.size())),
      columnCount_(0),
      options_(options)
    {

        vigra_precondition(options.training_set_proportion == 0.0 ||
                           options.training_set_size == 0,
            "RandomForestOptions: absolute and proprtional training set sizes "
            "cannot be specified at the same time.");
        vigra_precondition(classes_.size() > 1,
            "RandomForestOptions::weights(): need at least two classes.");
        vigra_precondition(options.class_weights.size() == 0 || options.class_weights.size() == classes_.size(),
            "RandomForestOptions::weights(): wrong number of classes.");
    }

	//Not understood yet
	//Does not use the options object but the columnCount object.
    template<class ClassLabelIterator, class TreeIterator, class WeightIterator>
    RandomForest(ClassLabelIterator cl, ClassLabelIterator cend,
                  unsigned int treeCount, unsigned int columnCount,
                  TreeIterator trees, WeightIterator weights)
    : classes_(cl, cend),
      trees_(treeCount, detail::DecisionTree(classes_.size())),
      columnCount_(columnCount)
    {
        for(unsigned int k=0; k<treeCount; ++k, ++trees, ++weights)
        {
            trees_[k].tree_ = *trees;
            trees_[k].terminalWeights_ = *weights;
        }
    }

    int featureCount() const
    {
        vigra_precondition(columnCount_ > 0,
           "RandomForest::featureCount(): Random forest has not been trained yet.");
        return columnCount_;
    }

    int labelCount() const
    {
        return classes_.size();
    }

    int treeCount() const
    {
        return trees_.size();
    }

    // loss == 0.0 means unweighted random forest
    template <class U, class C, class Array, class Random>
    double learn(MultiArrayView<2, U, C> const & features, Array const & labels,
               Random const& random);

    template <class U, class C, class Array>
    double learn(MultiArrayView<2, U, C> const & features, Array const & labels)
    {
        return learn(features, labels, RandomTT800::global());
    }

    template <class U, class C>
    ClassLabelType predictLabel(MultiArrayView<2, U, C> const & features) const;

    template <class U, class C1, class T, class C2>
    void predictLabels(MultiArrayView<2, U, C1> const & features,
                       MultiArrayView<2, T, C2> & labels) const
    {
        vigra_precondition(features.shape(0) == labels.shape(0),
            "RandomForest::predictLabels(): Label array has wrong size.");
        for(int k=0; k<features.shape(0); ++k)
            labels(k,0) = predictLabel(rowVector(features, k));
    }

    template <class U, class C, class Iterator>
    ClassLabelType predictLabel(MultiArrayView<2, U, C> const & features,
                                Iterator priors) const;

    template <class U, class C1, class T, class C2>
    void predictProbabilities(MultiArrayView<2, U, C1> const & features,
                              MultiArrayView<2, T, C2> & prob) const;

    template <class U, class C1, class T, class C2>
    void predictNodes(MultiArrayView<2, U, C1> const & features,
                                                   MultiArrayView<2, T, C2> & NodeIDs) const;
};

template <class ClassLabelType>
template <class U, class C1, class Array, class Random>
double
RandomForest<ClassLabelType>::learn(MultiArrayView<2, U, C1> const & features,
                                             Array const & labels,
                                             Random const& random)
{
    unsigned int classCount = classes_.size();
    unsigned int m = rowCount(features);
    unsigned int n = columnCount(features);
    vigra_precondition((unsigned int)(m) == (unsigned int)labels.size(),
      "RandomForest::learn(): Label array has wrong size.");

    vigra_precondition(options_.training_set_size <= m || options_.sample_with_replacement,
       "RandomForest::learn(): Requested training set size exceeds total number of examples.");

    MultiArrayIndex mtry = (options_.mtry == 0)
                                ? int(std::floor(std::sqrt(double(n)) + 0.5))
                                : options_.mtry;

    vigra_precondition(mtry <= (MultiArrayIndex)n,
       "RandomForest::learn(): mtry must be less than number of features.");

    MultiArrayIndex msamples = options_.training_set_size;
    if(options_.sample_classes_individually)
        msamples = int(std::ceil(double(msamples) / classCount));

    ArrayVector<int> intLabels(m), classExampleCounts(classCount);

	// verify the input labels
    int minClassCount;
    {
        typedef std::map<ClassLabelType, int > LabelChecker;
        typedef typename LabelChecker::iterator LabelCheckerIterator;
        LabelChecker labelChecker;
        for(unsigned int k=0; k<classCount; ++k)
            labelChecker[classes_[k]] = k;

        for(unsigned int k=0; k<m; ++k)
        {
            LabelCheckerIterator found = labelChecker.find(labels[k]);
            vigra_precondition(found != labelChecker.end(),
                "RandomForest::learn(): Unknown class label encountered.");
            intLabels[k] = found->second;
            ++classExampleCounts[intLabels[k]];
        }
        minClassCount = *argMin(classExampleCounts.begin(), classExampleCounts.end());
        vigra_precondition(minClassCount > 0,
             "RandomForest::learn(): At least one class is missing in the training set.");
        if(msamples > 0 && options_.sample_classes_individually &&
                          !options_.sample_with_replacement)
        {
            vigra_precondition(msamples <= minClassCount,
                "RandomForest::learn(): Too few examples in smallest class to reach "
                "requested training set size.");
        }
    }
    columnCount_ = n;
    ArrayVector<int> indices(m);
    for(unsigned int k=0; k<m; ++k)
        indices[k] = k;

    if(options_.sample_classes_individually)
    {
        detail::RandomForestLabelSorter<ArrayVector<int> > sorter(intLabels);
        std::sort(indices.begin(), indices.end(), sorter);
    }

    ArrayVector<int> usedIndices(m), oobCount(m), oobErrorCount(m);

    UniformIntRandomFunctor<Random> randint(0, m-1, random);
    //std::cerr << "Learning a RF \n";
    for(unsigned int k=0; k<trees_.size(); ++k)
    {
       //std::cerr << "Learning tree " << k << " ...\n";

        ArrayVector<int> trainingSet;
        usedIndices.init(0);

        if(options_.sample_classes_individually)
        {
            int first = 0;
            for(unsigned int l=0; l<classCount; ++l)
            {
                int lc = classExampleCounts[l];
                int lsamples = (msamples == 0)
                                   ? int(std::ceil(options_.training_set_proportion*lc))
                                   : msamples;

                if(options_.sample_with_replacement)
                {
                    for(int ll=0; ll<lsamples; ++ll)
                    {
                        trainingSet.push_back(indices[first+randint(lc)]);
                        ++usedIndices[trainingSet.back()];
                    }
                }
                else
                {
                    for(int ll=0; ll<lsamples; ++ll)
                    {
                        std::swap(indices[first+ll], indices[first+ll+randint(lc-ll)]);
                        trainingSet.push_back(indices[first+ll]);
                        ++usedIndices[trainingSet.back()];
                    }
                    //std::sort(indices.begin(), indices.begin()+lsamples);
                }
                first += lc;
            }
        }
        else
        {
            if(msamples == 0)
                msamples = int(std::ceil(options_.training_set_proportion*m));

            if(options_.sample_with_replacement)
            {
                for(int l=0; l<msamples; ++l)
                {
                    trainingSet.push_back(indices[randint(m)]);
                    ++usedIndices[trainingSet.back()];
                }
            }
            else
            {
                for(int l=0; l<msamples; ++l)
                {
                    std::swap(indices[l], indices[l+randint(m-l)/*oikas*/]);
                    trainingSet.push_back(indices[l]);
                    ++usedIndices[trainingSet.back()];
                }


            }

        }
        trees_[k].learn(features, intLabels,
                        trainingSet.begin(), trainingSet.size(),
                        options_.featuresPerNode(mtry), randint);
//        for(unsigned int l=0; l<m; ++l)
//        {
//            if(!usedIndices[l])
//            {
//                ++oobCount[l];
//                if(trees_[k].predictLabel(rowVector(features, l)) != intLabels[l])
//                    ++oobErrorCount[l];
//            }
//        }

        for(unsigned int l=0; l<m; ++l)
        {
            if(!usedIndices[l])
            {
                ++oobCount[l];
                if(trees_[k].predictLabel(rowVector(features, l)) != intLabels[l])
				{
                    ++oobErrorCount[l];
                    if(options_.oob_data.data() != 0)
                        options_.oob_data(l, k) = 2;
				}
				else if(options_.oob_data.data() != 0)
				{
					options_.oob_data(l, k) = 1;
				}
            }
        }
        // TODO: default value for oob_data
        // TODO: implement variable importance
        //if(!options_.sample_with_replacement){
        //std::cerr << "done\n";
        //trees_[k].print(std::cerr);
        #ifdef VIGRA_RF_VERBOSE
        trees_[k].printStatistics(std::cerr);
        #endif
    }
    double oobError = 0.0;
    int totalOobCount = 0;
    for(unsigned int l=0; l<m; ++l)
        if(oobCount[l])
        {
            oobError += double(oobErrorCount[l]) / oobCount[l];
            ++totalOobCount;
        }
    return oobError / totalOobCount;
}

template <class ClassLabelType>
template <class U, class C>
ClassLabelType
RandomForest<ClassLabelType>::predictLabel(MultiArrayView<2, U, C> const & features) const
{
    vigra_precondition(columnCount(features) >= featureCount(),
        "RandomForest::predictLabel(): Too few columns in feature matrix.");
    vigra_precondition(rowCount(features) == 1,
        "RandomForest::predictLabel(): Feature matrix must have a single row.");
    Matrix<double> prob(1, classes_.size());
    predictProbabilities(features, prob);
    return classes_[argMax(prob)];
}


//Same thing as above with priors for each label !!!
template <class ClassLabelType>
template <class U, class C, class Iterator>
ClassLabelType
RandomForest<ClassLabelType>::predictLabel(MultiArrayView<2, U, C> const & features,
                                           Iterator priors) const
{
    using namespace functor;
    vigra_precondition(columnCount(features) >= featureCount(),
        "RandomForest::predictLabel(): Too few columns in feature matrix.");
    vigra_precondition(rowCount(features) == 1,
        "RandomForest::predictLabel(): Feature matrix must have a single row.");
    Matrix<double> prob(1,classes_.size());
    predictProbabilities(features, prob);
    std::transform(prob.begin(), prob.end(), priors, prob.begin(), Arg1()*Arg2());
    return classes_[argMax(prob)];
}

template <class ClassLabelType>
template <class U, class C1, class T, class C2>
void
RandomForest<ClassLabelType>::predictProbabilities(MultiArrayView<2, U, C1> const & features,
                                                   MultiArrayView<2, T, C2> & prob) const
{

	//Features are n xp
	//prob is n x NumOfLabel probability for each feature in each class

    vigra_precondition(rowCount(features) == rowCount(prob),
      "RandomForest::predictProbabilities(): Feature matrix and probability matrix size mismatch.");

	// num of features must be bigger than num of features in Random forest training
	// but why bigger?
    vigra_precondition(columnCount(features) >= featureCount(),
      "RandomForest::predictProbabilities(): Too few columns in feature matrix.");
    vigra_precondition(columnCount(prob) == (MultiArrayIndex)labelCount(),
      "RandomForest::predictProbabilities(): Probability matrix must have as many columns as there are classes.");

	//Classify for each row.
    for(int row=0; row < rowCount(features); ++row)
    {
	//contains the weights returned by a single tree???
	//thought that one tree has only one vote???
	//Pruning???
        ArrayVector<double>::const_iterator weights;

        //totalWeight == totalVoteCount!
	double totalWeight = 0.0;

	//Set each VoteCount = 0 - prob(row,l) contains vote counts until
	//further normalisation
        for(unsigned int l=0; l<classes_.size(); ++l)
            prob(row, l) = 0.0;

	//Let each tree classify...
        for(unsigned int k=0; k<trees_.size(); ++k)
        {
		//get weights predicted by single tree
            weights = trees_[k].predict(rowVector(features, row));

		//update votecount.
            for(unsigned int l=0; l<classes_.size(); ++l)
            {
                prob(row, l) += detail::RequiresExplicitCast<T>::cast(weights[l]);
                //every weight in totalWeight.
                totalWeight += weights[l];
            }
        }

	//Normalise votes in each row by total VoteCount (totalWeight
        for(unsigned int l=0; l<classes_.size(); ++l)
                prob(row, l) /= detail::RequiresExplicitCast<T>::cast(totalWeight);
    }
}


template <class ClassLabelType>
template <class U, class C1, class T, class C2>
void
RandomForest<ClassLabelType>::predictNodes(MultiArrayView<2, U, C1> const & features,
                                                   MultiArrayView<2, T, C2> & NodeIDs) const
{
    vigra_precondition(columnCount(features) >= featureCount(),
      "RandomForest::getNodesRF(): Too few columns in feature matrix.");
    vigra_precondition(rowCount(features) <= rowCount(NodeIDs),
      "RandomForest::getNodesRF(): Too few rows in NodeIds matrix");
    vigra_precondition(columnCount(NodeIDs) >= treeCount(),
      "RandomForest::getNodesRF(): Too few columns in NodeIds matrix.");
    NodeIDs.init(0);
    for(unsigned int k=0; k<trees_.size(); ++k)
    {
        for(int row=0; row < rowCount(features); ++row)
        {
            NodeIDs(row,k) = trees_[k].leafID(rowVector(features, row));
        }
    }
}

} // namespace vigra

#undef RandomForest
#undef DecisionTreee_Deprec

#endif // VIGRA_RANDOM_FOREST_HXX

