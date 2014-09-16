/************************************************************************/
/*                                                                      */
/*    Copyright 2013-2014 by Ullrich Koethe                             */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
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

#ifndef VIGRA_SKELETON_HXX
#define VIGRA_SKELETON_HXX

#include <vector>
#include <set>
#include <map>
#include "vector_distance.hxx"
#include "iteratorfacade.hxx"
#include "pixelneighborhood.hxx"

namespace vigra
{

namespace detail {

template <class Node>
struct SkeletonNode
{
    Node parent, principal_child;
    double length, salience;
    
    SkeletonNode()
    : parent(lemon::INVALID)
    , principal_child(lemon::INVALID)
    , length(0.5)
    , salience(1.0)
    {}
    
    SkeletonNode(Node const & s)
    : parent(s)
    , principal_child(lemon::INVALID)
    , length(0.5)
    , salience(1.0)
    {}
};

template <class Node>
struct SkeletonRegion
{
    typedef SkeletonNode<Node>                 SNode;
    typedef std::map<Node, SNode>              Skeleton;
    
    Node anchor, lower, upper;
    Skeleton skeleton;
    
    SkeletonRegion()
    : anchor(lemon::INVALID)
    , lower(NumericTraits<MultiArrayIndex>::max())
    , upper(NumericTraits<MultiArrayIndex>::min())
    {}
    
    void addNode(Node const & n)
    {
        skeleton[n] = SNode(n);
        anchor = n;
        lower = min(lower, n);
        upper = max(upper, n);
    }
};

template <class Graph, class Node, class NodeMap>
inline unsigned int 
neighborhoodConfiguration(Graph const & g, Node const & n, NodeMap const & labels)
{
    typedef typename Graph::OutArcIt      ArcIt;
    typedef typename NodeMap::value_type  LabelType;
    
    LabelType label = abs(labels[n]);
    unsigned int v = 0;
    for(ArcIt arc(g, n); arc != lemon::INVALID; ++arc)
    {
        v = (v << 1) | (abs(labels[g.target(*arc)]) == label ? 1 : 0);
    }
    return v;
}

template <class Node, class Cost>
struct SkeletonSimplePoint
{
    Node point;
    Cost cost;

    SkeletonSimplePoint(Node const & p, Cost c)
    : point(p), cost(c)
    {}

    bool operator<(SkeletonSimplePoint const & o) const
    {
        return cost < o.cost;
    }

    bool operator>(SkeletonSimplePoint const & o) const
    {
        return cost > o.cost;
    }
};

template <class Graph, class CostMap, class LabelMap>
void
skeletonThinning(Graph const & g, CostMap const & cost, LabelMap & labels)
{
    typedef typename Graph::Node           Node;
    typedef typename Graph::NodeIt         NodeIt;
    typedef typename Graph::OutArcIt       ArcIt;
    typedef typename LabelMap::value_type  LabelType;

    bool isSimplePoint[256] = {
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 
        0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
        1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 
        1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 
        1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 
        1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0
    };

    typedef SkeletonSimplePoint<Node, double> SP;
    // use std::greater because we need the smallest gradients at the top of the queue
    std::priority_queue<SP, std::vector<SP>, std::greater<SP> >  pqueue;
    // std::priority_queue<SP>  pqueue;

    int max_degree = g.maxDegree();
    for (NodeIt node(g); node != lemon::INVALID; ++node)
    {
        Node p = *node;
        if(g.out_degree(p) == max_degree &&
           labels[p] > 0 &&
           isSimplePoint[neighborhoodConfiguration(g, p, labels)])
        {
            pqueue.push(SP(p, cost[p]));
            labels[p] = -labels[p]; // remember that it is already in queue
        }
    }

    while(pqueue.size())
    {
        Node p = pqueue.top().point;
        pqueue.pop();

        if(!isSimplePoint[neighborhoodConfiguration(g, p, labels)])
        {
            labels[p] = abs(labels[p]);
            continue; // point may no longer be simple because its neighbors changed
        }
        
        labels[p] = 0; // delete simple point

        for (ArcIt arc(g, p); arc != lemon::INVALID; ++arc)
        {
            Node q = g.target(*arc);
            if(g.out_degree(q) == max_degree &&
               labels[q] > 0 &&
               isSimplePoint[neighborhoodConfiguration(g, q, labels)])
            {
                pqueue.push(SP(q, cost[q]));
                labels[q] = -labels[q]; // remember that it is already in queue
            }
        }
    }
}


template <class DestIterator>
int neighborhoodConfiguration(DestIterator dul)
{
    int v = 0;
    NeighborhoodCirculator<DestIterator, EightNeighborCode> c(dul, EightNeighborCode::SouthEast);
    for(int i=0; i<8; ++i, --c)
    {
        v = (v << 1) | ((*c != 0) ? 1 : 0);
    }

    return v;
}

} // namespace detail 

// FIXME: support pruneTopology()

struct SkeletonOptions
{
    enum SkeletonMode {
         DontPrune,
         PruneCenterLine,
         PruneLength,
         PruneLengthRelative,
         ReturnLength,
         PruneSalience,
         PruneSalienceRelative,
         ReturnSalience
    };
    
    SkeletonMode mode;
    double pruning_threshold;
    
    SkeletonOptions()
    : mode(PruneSalienceRelative)
    , pruning_threshold(0.2)
    {}
    
    SkeletonOptions & dontPrune()
    {
        mode = DontPrune;
        return *this;
    }
    
    SkeletonOptions & pruneCenterLine()
    {
        mode = PruneCenterLine;
        return *this;
    }
    
    SkeletonOptions & pruneLength(double threshold)
    {
        mode = PruneLength;
        pruning_threshold = threshold;
        return *this;
    }
    
    SkeletonOptions & pruneLengthRelative(double threshold)
    {
        mode = PruneLengthRelative;
        pruning_threshold = threshold;
        return *this;
    }
    
    SkeletonOptions & returnLength()
    {
        mode = ReturnLength;
        return *this;
    }
    
    SkeletonOptions & pruneSalience(double threshold)
    {
        mode = PruneSalience;
        pruning_threshold = threshold;
        return *this;
    }
    
    SkeletonOptions & pruneSalienceRelative(double threshold)
    {
        mode = PruneSalienceRelative;
        pruning_threshold = threshold;
        return *this;
    }
    
    SkeletonOptions & returnSalience()
    {
        mode = ReturnSalience;
        return *this;
    }
    
    bool isLengthMode() const
    {
        return mode == PruneLength ||
               mode == PruneLengthRelative ||
               mode == ReturnLength;

    }
    
    bool isSalienceMode() const
    {
        return mode == PruneSalience ||
               mode == PruneSalienceRelative ||
               mode == ReturnSalience;

    }
    
    bool isRelativePruning() const
    {
        return mode == PruneLengthRelative ||
               mode == PruneSalienceRelative;

    }
};

/********************************************************/
/*                                                      */
/*                       skeleton                       */
/*                                                      */
/********************************************************/

// template <unsigned int N, class T1, class S1,
                          // class T2, class S2>
// void
// skeletonize(MultiArrayView<N, T1, S1> const & labels,
            // MultiArrayView<N, T2, S2> dest,
            // SkeletonOptions const & options = SkeletonOptions())
// {

template <class T1, class S1,
          class T2, class S2>
void
skeletonize(MultiArrayView<2, T1, S1> const & labels,
            MultiArrayView<2, T2, S2> dest,
            SkeletonOptions const & options = SkeletonOptions())
{
    using namespace functor;
    static const unsigned int N = 2;
    typedef typename MultiArrayShape<N>::type  Shape;
    typedef GridGraph<N>                       Graph;
    typedef typename Graph::Node               Node;
    typedef typename Graph::NodeIt             NodeIt;
    typedef typename Graph::EdgeIt             EdgeIt;
    typedef typename Graph::OutBackArcIt       neighbor_iterator;
    typedef double                             WeightType;
    typedef detail::SkeletonNode<Node>         SNode;
    typedef std::map<Node, SNode>              Skeleton;
    
    vigra_precondition(labels.shape() == dest.shape(),
        "skeleton(): shape mismatch between input and output.");
        
    MultiArray<N, Shape> vectors(labels.shape());
    boundaryVectorDistance(labels, vectors, false, OuterBoundary);
    
    dest = 0;
    T1 maxLabel = 0;
    // find skeleton points
    {
        Graph g(labels.shape());
        for (EdgeIt edge(g); edge != lemon::INVALID; ++edge)
        {
            Node p1 = g.u(*edge),
                 p2 = g.v(*edge);
            T1 l1 = labels[p1], 
               l2 = labels[p2];
            maxLabel = max(maxLabel, max(l1, l2));
            if(l1 == l2)
            {
                if(l1 <= 0)  // only consider positive labels
                    continue;
                
                Node b1 = vectors[p1] + p1,
                     b2 = vectors[p2] + p2,
                     diff = b2 - b1,
                     dp = p2 - p1;
                if(max(abs(diff)) <= 1) // points whose support points are equal or adjacent
                    continue;           // don't belong to the skeleton
                
                // criterion 1: among p1 and p2, the point which is closer to the bisector 
                //              of the support points b1 and b2 belongs to the skeleton
                MultiArrayIndex d = squaredNorm(diff);
                MultiArrayIndex d1 = abs(2*dot(diff, p1) - d),
                                d2 = abs(2*dot(diff, p2) - d);
                if(d1 <= d2) 
                    dest[p1] = l1;
                else
                    dest[p2] = l2;
                
                // criterion 2: among p1 and p2, the point which is closer to the intersection
                //              between the bisector of the support points b1 and b2 
                //              and the grid line connecting p1 and p2 belongs to the skeleton
                // MultiArrayIndex d1 = dot(diff, dp),
                                // d2 = dot(diff, vectors[p1]+vectors[p2]);
                // if(d1*d2 <= 0)
                    // dest[p1] = l1;
                // else
                    // dest[p2] = l2;
            }
            else
            {
                if(l1 > 0 && 
                   max(abs(vectors[p1] + p1 - p2)) > 1)
                        dest[p1] = l1;
                if(l2 > 0 &&
                   max(abs(vectors[p2] + p2 - p1)) > 1)
                        dest[p2] = l2;

                // if(l1 > 0 && 
                  // squaredNorm(vectors[p1] + p1 - p2) > 1)
                        // dest[p1] = l1;
                // if(l2 > 0 &&
                   // squaredNorm(vectors[p2] + p2 - p1) > 1)
                        // dest[p2] = l2;
            }
        }
    }
    
    // {
        // using namespace multi_math;
        // MultiArray<2, int> dist = squaredNorm(vectors);
        // Graph g(labels.shape(), IndirectNeighborhood);
        // detail::skeletonThinning(g, dist, dest);
    // }
    
    if(options.mode == SkeletonOptions::DontPrune)
        return;
        
    ArrayVector<detail::SkeletonRegion<Node> > regions((size_t)maxLabel + 1);
    Graph g(labels.shape(), IndirectNeighborhood);
    typename Graph::template EdgeMap<WeightType> weights(g);
    for (NodeIt node(g); node != lemon::INVALID; ++node)
    {
        Node p1 = *node;
        T2 label = dest[p1];
        if(label > 0)
            regions[(size_t)label].addNode(p1);

        for (neighbor_iterator arc(g, p1); arc != lemon::INVALID; ++arc)
        {
            Node p2 = g.target(*arc);
            if(label > 0 && dest[p2] == label)
                weights[*arc] = norm(p1-p2);
            else
                weights[*arc] = NumericTraits<WeightType>::max();
        }
    }
    
    dest = 0;
    ShortestPathDijkstra<Graph, WeightType> pathFinder(g);
    WeightType maxWeight = g.edgeNum()*sqrt(N);
    for(unsigned int label=1; label < regions.size(); ++label)
    {
        Skeleton & skeleton = regions[label].skeleton;
        if(skeleton.size() == 0) // label doesn't exist
            continue;
        
        Node anchor = regions[label].anchor,
             lower  = regions[label].lower,
             upper  = regions[label].upper + Shape(1);
        for(int k=0; k < 2; ++k) // two iterations suffice in a tree, FIXME: check if skeleton has loops
        {
            pathFinder.run(weights, anchor, lemon::INVALID, maxWeight, lower, upper);
            anchor = pathFinder.target();
        }
        
        Polygon<Shape> center_line;
        center_line.push_back_unsafe(anchor);
        while(pathFinder.predecessors()[center_line.back()] != center_line.back())
            center_line.push_back_unsafe(pathFinder.predecessors()[center_line.back()]);
            
        if(options.mode == SkeletonOptions::PruneCenterLine)
        {
            for(unsigned int k=0; k<center_line.size(); ++k)
                dest[center_line[k]] = (T2)label;
            continue; // to next label
        }
        
        Node center = center_line[roundi(center_line.arcLengthQuantile(0.5))];
        pathFinder.run(weights, center, lemon::INVALID, maxWeight, lower, upper);
        
        bool compute_salience = options.isSalienceMode();
        ArrayVector<Node> const & raw_skeleton = pathFinder.discoveryOrder();
        // from periphery to center: create skeleton tree and compute salience
        for(int k=raw_skeleton.size()-1; k >= 0; --k)
        {
            Node p1 = raw_skeleton[k],
                 p2 = pathFinder.predecessors()[p1];
            SNode & n1 = skeleton[p1];
            SNode & n2 = skeleton[p2];
            n1.parent = p2;

            WeightType l = n1.length + norm(p1-p2);
            if(!compute_salience)
            {
                n1.salience = n1.length;
                // propagate length to parent if this is the longest subtree
                if(n2.length < l)
                {
                    n2.length = l;
                    n2.principal_child = p1;
                }
            }
            else 
            {
                const double min_length = 4.0; // salience is meaningless for shorter segments due
                                               // to quantization noise (staircasing) of the boundary
                // propagate length to parent if this is the longest subtree
                if(n2.length < l)
                {
                    n2.length = l;
                    if(n1.length < min_length)
                        n2.principal_child = p1;
                }
                // compute salience
                if(n1.length >= min_length)
                {
                    n1.salience = max(n1.salience, (n1.length + 0.5) / norm(vectors[p1]));
                
                    // propagate salience to parent if this is the most salient subtree
                    if(n2.salience < n1.salience)
                    {
                        n2.salience = n1.salience;
                        n2.principal_child = p1;
                    }
                }
            }
        }
        
        bool dont_prune = options.mode == SkeletonOptions::ReturnLength ||
                          options.mode == SkeletonOptions::ReturnSalience;
        WeightType threshold = options.isRelativePruning()
                                   ? options.pruning_threshold*skeleton[center].salience
                                   : options.pruning_threshold;
        // from center to periphery: propagate salience
        for(int k=0; k < (int)raw_skeleton.size(); ++k)
        {
            Node p1 = raw_skeleton[k];
            SNode & n1 = skeleton[p1];
            Node p2 = n1.parent;
            SNode & n2 = skeleton[p2];
            if(p1 == n2.principal_child)
                n1.salience = n2.salience;
            if(dont_prune)
                dest[p1] = n1.salience;
            else if(n1.salience >= threshold)
                dest[p1] = (T2)label;
        }
    }
}


// a truly reliable skeleton is only possible by considering the 
// interpixel boundaries between the voronoi regions of the 
// vector distance transform
template <class T1, class S1,
          class T2, class S2>
void
skeletonizeInterpixel(MultiArrayView<2, T1, S1> const & labels,
            MultiArrayView<2, T2, S2> dest,
            SkeletonOptions const & options = SkeletonOptions())
{
    using namespace functor;
    static const unsigned int N = 2;
    typedef typename MultiArrayShape<N>::type  Shape;
    typedef GridGraph<N>                       Graph;
    typedef typename Graph::Node               Node;
    typedef typename Graph::NodeIt             NodeIt;
    typedef typename Graph::EdgeIt             EdgeIt;
    typedef typename Graph::OutBackArcIt       neighbor_iterator;
    typedef double                             WeightType;
    typedef detail::SkeletonNode<Node>         SNode;
    typedef std::map<Node, SNode>              Skeleton;
    
    vigra_precondition(labels.shape() == dest.shape(),
        "skeleton(): shape mismatch between input and output.");
        
    MultiArray<N, Shape> vectors(labels.shape());
    boundaryVectorDistance(labels, vectors, false, OuterBoundary);
    
    MultiArray<N, Shape> vectors2(labels.shape());
    for(auto it = createCoupledIterator(vectors, vectors2, labels); it.isValid(); ++it)
        get<2>(*it) = get<3>(*it) > 0
                         ? get<1>(*it) + get<0>(*it)
                         : Shape(lemon::INVALID);
    MultiArray<N, UInt32> voronoi(labels.shape());
    labelMultiArrayWithBackground(vectors2, voronoi, IndirectNeighborhood, Shape(lemon::INVALID),
                                  [](Node const & p, Node const & q) {return max(abs(p-q)) <= 1; });
    exportImage(voronoi, "voronoi.tif");
    
    dest = 0;
    T1 maxLabel = 0;
    // find skeleton points
    {
        Graph g(labels.shape());
        for (EdgeIt edge(g); edge != lemon::INVALID; ++edge)
        {
            Node p1 = g.u(*edge),
                 p2 = g.v(*edge);
            T1 l1 = labels[p1], 
               l2 = labels[p2];
            maxLabel = max(maxLabel, max(l1, l2));
            if(l1 == l2)
            {
                if(l1 <= 0)  // only consider positive labels
                    continue;
                
                Node b1 = vectors[p1] + p1,
                     b2 = vectors[p2] + p2,
                     diff = b2 - b1;
                if(max(abs(diff)) <= 1) // points whose support points are equal or adjacent
                    continue;           // don't belong to the skeleton
                MultiArrayIndex d = squaredNorm(diff);
                MultiArrayIndex d1 = abs(2*dot(diff, p1) - d),
                                d2 = abs(2*dot(diff, p2) - d);
                if(d1 <= d2) // the point closer to the bisector of b1 and b2 belongs to the skeleton
                    dest[p1] = l1;
                else
                    dest[p2] = l2;
                
                // Node b1 = vectors[p1] + p1,
                     // b2 = vectors[p2] + p2,
                     // dp = p2 - p1,
                     // db = b2 - b1;
                // if(max(abs(db)) <= 1) // points whose support points are equal or adjacent
                    // continue;           // don't belong to the skeleton
                // // if(squaredNorm(db) <= 1) // points whose support points are equal or adjacent
                    // // continue;           // don't belong to the skeleton
                
                // MultiArrayIndex bp = dot(db, dp),
                                // d  = dot(db, vectors[p1]+vectors[p2]);
                // if(d*bp <= 0)
                    // dest[p1] = l1;
                // else
                    // dest[p2] = l2;
            }
            else
            {
                if(l1 > 0 && 
                   max(abs(vectors[p1] + p1 - p2)) > 1)
                        dest[p1] = l1;
                if(l2 > 0 &&
                   max(abs(vectors[p2] + p2 - p1)) > 1)
                        dest[p2] = l2;

                // if(l1 > 0 && 
                  // squaredNorm(vectors[p1] + p1 - p2) > 1)
                        // dest[p1] = l1;
                // if(l2 > 0 &&
                   // squaredNorm(vectors[p2] + p2 - p1) > 1)
                        // dest[p2] = l2;
            }
        }
    }
    
    {
        using namespace multi_math;
        MultiArray<2, int> dist = squaredNorm(vectors);
        Graph g(labels.shape(), IndirectNeighborhood);
        detail::skeletonThinning(g, dist, dest);
    }
    
    if(options.mode == SkeletonOptions::DontPrune)
        return;
        
    ArrayVector<detail::SkeletonRegion<Node> > regions((size_t)maxLabel + 1);
    Graph g(labels.shape(), IndirectNeighborhood);
    typename Graph::template EdgeMap<WeightType> weights(g);
    for (NodeIt node(g); node != lemon::INVALID; ++node)
    {
        Node p1 = *node;
        T2 label = dest[p1];
        if(label > 0)
            regions[(size_t)label].addNode(p1);

        for (neighbor_iterator arc(g, p1); arc != lemon::INVALID; ++arc)
        {
            Node p2 = g.target(*arc);
            if(label > 0 && dest[p2] == label)
                weights[*arc] = norm(p1-p2);
            else
                weights[*arc] = NumericTraits<WeightType>::max();
        }
    }
    
    dest = 0;
    ShortestPathDijkstra<Graph, WeightType> pathFinder(g);
    WeightType maxWeight = g.edgeNum()*sqrt(N);
    for(unsigned int label=1; label < regions.size(); ++label)
    {
        Skeleton & skeleton = regions[label].skeleton;
        if(skeleton.size() == 0) // label doesn't exist
            continue;
        
        Node anchor = regions[label].anchor,
             lower  = regions[label].lower,
             upper  = regions[label].upper + Shape(1);
        for(int k=0; k < 2; ++k) // two iterations suffice in a tree, FIXME: check if skeleton has loops
        {
            pathFinder.run(weights, anchor, lemon::INVALID, maxWeight, lower, upper);
            anchor = pathFinder.target();
        }
        
        Polygon<Shape> center_line;
        center_line.push_back_unsafe(anchor);
        while(pathFinder.predecessors()[center_line.back()] != center_line.back())
            center_line.push_back_unsafe(pathFinder.predecessors()[center_line.back()]);
            
        if(options.mode == SkeletonOptions::PruneCenterLine)
        {
            for(unsigned int k=0; k<center_line.size(); ++k)
                dest[center_line[k]] = (T2)label;
            continue; // to next label
        }
        
        Node center = center_line[roundi(center_line.arcLengthQuantile(0.5))];
        pathFinder.run(weights, center, lemon::INVALID, maxWeight, lower, upper);
        
        bool compute_salience = options.isSalienceMode();
        ArrayVector<Node> const & raw_skeleton = pathFinder.discoveryOrder();
        // from periphery to center: create skeleton tree and compute salience
        for(int k=raw_skeleton.size()-1; k >= 0; --k)
        {
            Node p1 = raw_skeleton[k],
                 p2 = pathFinder.predecessors()[p1];
            SNode & n1 = skeleton[p1];
            SNode & n2 = skeleton[p2];
            n1.parent = p2;

            WeightType l = n1.length + norm(p1-p2);
            if(!compute_salience)
            {
                n1.salience = n1.length;
                // propagate length to parent if this is the longest subtree
                if(n2.length < l)
                {
                    n2.length = l;
                    n2.principal_child = p1;
                }
            }
            else 
            {
                const double min_length = 4.0; // salience is meaningless for shorter segments due
                                               // to quantization noise (staircasing) of the boundary
                // propagate length to parent if this is the longest subtree
                if(n2.length < l)
                {
                    n2.length = l;
                    if(n1.length < min_length)
                        n2.principal_child = p1;
                }
                // compute salience
                if(n1.length >= min_length)
                {
                    n1.salience = max(n1.salience, (n1.length + 0.5) / norm(vectors[p1]));
                
                    // propagate salience to parent if this is the most salient subtree
                    if(n2.salience < n1.salience)
                    {
                        n2.salience = n1.salience;
                        n2.principal_child = p1;
                    }
                }
            }
        }
        
        bool dont_prune = options.mode == SkeletonOptions::ReturnLength ||
                          options.mode == SkeletonOptions::ReturnSalience;
        WeightType threshold = options.isRelativePruning()
                                   ? options.pruning_threshold*skeleton[center].salience
                                   : options.pruning_threshold;
        // from center to periphery: propagate salience
        for(int k=0; k < (int)raw_skeleton.size(); ++k)
        {
            Node p1 = raw_skeleton[k];
            SNode & n1 = skeleton[p1];
            Node p2 = n1.parent;
            SNode & n2 = skeleton[p2];
            if(p1 == n2.principal_child)
                n1.salience = n2.salience;
            if(dont_prune)
                dest[p1] = n1.salience;
            else if(n1.salience >= threshold)
                dest[p1] = (T2)label;
        }
    }
}

} //-- namespace vigra

#endif        //-- VIGRA_SKELETON_HXX
