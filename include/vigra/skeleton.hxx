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
#include "graph_algorithms.hxx"

namespace vigra
{

namespace detail {

template <class Node>
struct SkeletonNode
{
    Node parent, principal_child;
    double length, salience;
    MultiArrayIndex partial_area;
    bool is_loop;
    
    SkeletonNode()
    : parent(lemon::INVALID)
    , principal_child(lemon::INVALID)
    , length(0.0)
    , salience(1.0)
    , partial_area(0)
    , is_loop(false)
    {}
    
    SkeletonNode(Node const & s)
    : parent(s)
    , principal_child(lemon::INVALID)
    , length(0.0)
    , salience(1.0)
    , partial_area(0)
    , is_loop(false)
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
    
    LabelType label = labels[n];
    unsigned int v = 0;
    for(ArcIt arc(g, n); arc != lemon::INVALID; ++arc)
    {
        v = (v << 1) | (labels[g.target(*arc)] == label ? 1 : 0);
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

template <class CostMap, class LabelMap>
void
skeletonThinning(CostMap const & cost, LabelMap & labels,
                 bool preserve_endpoints=true)
{
    typedef GridGraph<CostMap::actual_dimension> Graph;
    typedef typename Graph::Node           Node;
    typedef typename Graph::NodeIt         NodeIt;
    typedef typename Graph::OutArcIt       ArcIt;

    Graph g(labels.shape(), IndirectNeighborhood);
    typedef SkeletonSimplePoint<Node, double> SP;
    // use std::greater because we need the smallest distances at the top of the queue
    // (std::priority_queue is a max queue by default)
    std::priority_queue<SP, std::vector<SP>, std::greater<SP> >  pqueue;
    
    bool isSimpleStrong[256] = {
        0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 
        0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 
        1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 
        1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 
        1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 
        1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 
        1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 
        1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 
    };
    
    bool isSimplePreserveEndPoints[256] = {
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 
        1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    
    bool * isSimplePoint = preserve_endpoints
                               ? isSimplePreserveEndPoints
                               : isSimpleStrong;
    
    int max_degree = g.maxDegree();
    double epsilon = 0.5/labels.size(), offset = 0;
    for (NodeIt node(g); node != lemon::INVALID; ++node)
    {
        Node p = *node;
        if(g.out_degree(p) == max_degree &&
           labels[p] > 0 &&
           isSimplePoint[neighborhoodConfiguration(g, p, labels)])
        {
            pqueue.push(SP(p, cost[p]+offset));
	    offset += epsilon;
        }
    }

    while(pqueue.size())
    {
        Node p = pqueue.top().point;
        pqueue.pop();

        if(labels[p] == 0 ||
           !isSimplePoint[neighborhoodConfiguration(g, p, labels)])
        {
            continue; // point already deleted or no longer simple
        }
        
        labels[p] = 0; // delete simple point

        for (ArcIt arc(g, p); arc != lemon::INVALID; ++arc)
        {
            Node q = g.target(*arc);
            if(g.out_degree(q) == max_degree &&
               labels[q] > 0 &&
               isSimplePoint[neighborhoodConfiguration(g, q, labels)])
            {
                pqueue.push(SP(q, cost[q]+offset));
		offset += epsilon;
            }
        }
    }
}

template <class Label, class Labels>
struct CheckForHole
{
    Label label;
    Labels const & labels;
    
    CheckForHole(Label l, Labels const & ls)
    : label(l)
    , labels(ls)
    {}
    
    template <class Node>
    bool operator()(Node const & n) const
    {
        return labels[n] == label;
    }
};

} // namespace detail 

/** \addtogroup MultiArrayDistanceTransform
*/
//@{

    /** \brief Option object for \ref skeletonizeImage()
    */
struct SkeletonOptions
{
    enum SkeletonMode {
         DontPrune = 0,
         Prune = 1,
         Relative = 2,
         PreserveTopology = 4,
         Length = 8,
         Salience = 16,
         PruneCenterLine = 32,
         PruneLength = Length + Prune,
         PruneLengthRelative = PruneLength + Relative,
         PruneSalience = Salience + Prune,
         PruneSalienceRelative = PruneSalience + Relative,
         PruneTopology = PreserveTopology + Prune
    };
    
    SkeletonMode mode;
    double pruning_threshold;
    
        /** \brief construct with default settings
        
            (default: <tt>pruneSalienceRelative(0.2, true)</tt>)
        */
    SkeletonOptions()
    : mode(SkeletonMode(PruneSalienceRelative | PreserveTopology))
    , pruning_threshold(0.2)
    {}
    
        /** \brief return the un-pruned skeletong
        */
    SkeletonOptions & dontPrune()
    {
        mode = DontPrune;
        return *this;
    }
    
        /** \brief return only the region's center line (i.e. skeleton graph diameter)
        */
    SkeletonOptions & pruneCenterLine()
    {
        mode = PruneCenterLine;
        return *this;
    }
    
        /** \brief Don't prune and return the length of each skeleton segment.
        */
    SkeletonOptions & returnLength()
    {
        mode = Length;
        return *this;
    }
    
        /** \brief prune skeleton segments whose length is below the given threshold
        
            If \a preserve_topology is <tt>true</tt> (default), skeleton loops
            (i.e. parts enclosing a hole in the region) are preserved even if their 
            length is below the threshold. Otherwise, loops are pruned as well.
        */
    SkeletonOptions & pruneLength(double threshold, bool preserve_topology=true)
    {
        mode = PruneLength;
        if(preserve_topology)
            mode = SkeletonMode(mode | PreserveTopology);
        pruning_threshold = threshold;
        return *this;
    }
    
        /** \brief prune skeleton segments whose relative length is below the given threshold
        
            This works like <tt>pruneLength()</tt>, but the threshold is specified as a 
            fraction of the maximum segment length in the skeleton.
        */
    SkeletonOptions & pruneLengthRelative(double threshold, bool preserve_topology=true)
    {
        mode = PruneLengthRelative;
        if(preserve_topology)
            mode = SkeletonMode(mode | PreserveTopology);
        pruning_threshold = threshold;
        return *this;
    }
    
        /** \brief Don't prune and return the salience of each skeleton segment.
        */
    SkeletonOptions & returnSalience()
    {
        mode = Salience;
        return *this;
    }
    
        /** \brief prune skeleton segments whose salience is below the given threshold
        
            If \a preserve_topology is <tt>true</tt> (default), skeleton loops
            (i.e. parts enclosing a hole in the region) are preserved even if their 
            salience is below the threshold. Otherwise, loops are pruned as well.
        */
    SkeletonOptions & pruneSalience(double threshold, bool preserve_topology=true)
    {
        mode = PruneSalience;
        if(preserve_topology)
            mode = SkeletonMode(mode | PreserveTopology);
        pruning_threshold = threshold;
        return *this;
    }
    
        /** \brief prune skeleton segments whose relative salience is below the given threshold
        
            This works like <tt>pruneSalience()</tt>, but the threshold is specified as a 
            fraction of the maximum segment salience in the skeleton.
        */
    SkeletonOptions & pruneSalienceRelative(double threshold, bool preserve_topology=true)
    {
        mode = PruneSalienceRelative;
        if(preserve_topology)
            mode = SkeletonMode(mode | PreserveTopology);
        pruning_threshold = threshold;
        return *this;
    }
    
        /** \brief prune such that only the topology is preserved
        
            If \a preserve_center is <tt>true</tt> (default), the eccentricity center
            of the skeleton will not be pruned, even if it is not essential for the topology.
            Otherwise, the center is only preserved if it is essential. The center is always 
            preserved (and is the only remaining point) when the region has no holes.
        */
    SkeletonOptions & pruneTopology(bool preserve_center=true)
    {
        if(preserve_center)
            mode = PruneTopology;
        else
            mode = Prune;
        return *this;
    }
};

template <class T1, class S1,
          class T2, class S2,
          class ArrayLike>
void
skeletonizeImageImpl(MultiArrayView<2, T1, S1> const & labels,
                MultiArrayView<2, T2, S2> dest,
                ArrayLike * features,
                SkeletonOptions const & options)
{
    static const unsigned int N = 2;
    typedef typename MultiArrayShape<N>::type  Shape;
    typedef GridGraph<N>                       Graph;
    typedef typename Graph::Node               Node;
    typedef typename Graph::NodeIt             NodeIt;
    typedef typename Graph::EdgeIt             EdgeIt;
    typedef typename Graph::OutArcIt           ArcIt;
    typedef typename Graph::OutBackArcIt       BackArcIt;
    typedef double                             WeightType;
    typedef detail::SkeletonNode<Node>         SNode;
    typedef std::map<Node, SNode>              Skeleton;
    
    vigra_precondition(labels.shape() == dest.shape(),
        "skeleton(): shape mismatch between input and output.");
        
    MultiArray<N, MultiArrayIndex> squared_distance;
    dest = 0;
    T1 maxLabel = 0;
    // find skeleton points
    {
        using namespace multi_math;
        
        MultiArray<N, Shape> vectors(labels.shape());
        boundaryVectorDistance(labels, vectors, false, OuterBoundary);
        squared_distance = squaredNorm(vectors);
    
        ArrayVector<Node> ends_to_be_checked;
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
                
                const Node v1 = vectors[p1],
                           v2 = vectors[p2],
                           dp = p2 - p1,
                           dv = v2 - v1 + dp;
                if(max(abs(dv)) <= 1)  // points whose support points coincide or are adjacent
                    continue;          // don't belong to the skeleton
                
                // among p1 and p2, the point which is closer to the bisector 
                // of the support points p1 + v1 and p2 + v2 belongs to the skeleton
                const MultiArrayIndex d1 = dot(dv, dp),
                                      d2 = dot(dv, v1+v2);
                if(d1*d2 <= 0)
                {
                    dest[p1] = l1;
                    if(squared_distance[p1] == 4)
                        ends_to_be_checked.push_back(p1);
                }
                else
                {
                    dest[p2] = l2;
                    if(squared_distance[p2] == 4)
                        ends_to_be_checked.push_back(p2);
                }
            }
            else
            {
                if(l1 > 0 && max(abs(vectors[p1] + p1 - p2)) > 1)
                    dest[p1] = l1;
                if(l2 > 0 && max(abs(vectors[p2] + p2 - p1)) > 1)
                    dest[p2] = l2;
            }
        }
        
        
        // add a point when a skeleton line stops short of the shape boundary
        // FIXME: can this be solved during the initial detection above?
        Graph g8(labels.shape(), IndirectNeighborhood);
        for (unsigned k=0; k<ends_to_be_checked.size(); ++k)
        {
            // The phenomenon only occurs at points whose distance from the background is 2.
            // We've put these points into ends_to_be_checked. 
            Node p1 = ends_to_be_checked[k];
            T2 label = dest[p1];
            int count = 0;
            for (ArcIt arc(g8, p1); arc != lemon::INVALID; ++arc)
            {
                Node p2 = g8.target(*arc);
                if(dest[p2] == label && squared_distance[p2] < 4)
                    ++count;
            }
            if(count == 0) // point p1 has no neighbor at the boundary => activate one
                dest[p1+vectors[p1]/2] = label;
        }

        // from here on, we only need the squared DT, not the vector DT
    }
    
    // The true skeleton is defined by the interpixel edges between the
    // Voronoi regions of the DT. Our skeleton detection algorithm affectively
    // rounds the interpixel edges to the nearest pixel such that the result
    // is mainly 8-connected and thin. However, thick skeleton pieces may still 
    // arise when two interpixel contours are only one pixel apart, i.e. a
    // Voronoi region is only one pixel wide. Since this happens rarely, we 
    // can simply remove these cases by thinning.
    detail::skeletonThinning(squared_distance, dest);
    
    if(options.mode == SkeletonOptions::PruneCenterLine)
        dest = 0;
        
    // Reduce the full grid graph to a skeleton graph by inserting infinite
    // edge weights between skeleton pixels and non-skeleton pixels.
    if(features)
        features->resize((size_t)maxLabel + 1);
    ArrayVector<detail::SkeletonRegion<Node> > regions((size_t)maxLabel + 1);
    Graph g(labels.shape(), IndirectNeighborhood);
    WeightType maxWeight = g.edgeNum()*sqrt(N),
               infiniteWeight = 0.5*NumericTraits<WeightType>::max();
    typename Graph::template EdgeMap<WeightType> weights(g);
    for (NodeIt node(g); node != lemon::INVALID; ++node)
    {
        Node p1 = *node;
        T2 label = dest[p1];
        if(label <= 0)
            continue;
            
        // FIXME: consider using an AdjacencyListGraph from here on
        regions[(size_t)label].addNode(p1);

        for (ArcIt arc(g, p1); arc != lemon::INVALID; ++arc)
        {
            Node p2 = g.target(*arc);
            if(dest[p2] == label)
                weights[*arc] = norm(p1-p2);
            else
                weights[*arc] = infiniteWeight;
        }
    }
    
    ShortestPathDijkstra<Graph, WeightType> pathFinder(g);
    // Handle the skeleton of each region individually.
    for(std::size_t label=1; label < regions.size(); ++label)
    {
        Skeleton & skeleton = regions[label].skeleton;
        if(skeleton.size() == 0) // label doesn't exist
            continue;
        
        // Find a diameter (longest path) in the skeleton.
        Node anchor = regions[label].anchor,
             lower  = regions[label].lower,
             upper  = regions[label].upper + Shape(1);
             
        pathFinder.run(lower, upper, weights, anchor, lemon::INVALID, maxWeight);
        anchor = pathFinder.target();
        pathFinder.reRun(weights, anchor, lemon::INVALID, maxWeight);
        anchor = pathFinder.target();
        
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
        
        // Perform the eccentricity transform of the skeleton
        Node center = center_line[roundi(center_line.arcLengthQuantile(0.5))];
        pathFinder.reRun(weights, center, lemon::INVALID, maxWeight);
        
        bool compute_salience = (options.mode & SkeletonOptions::Salience) != 0;
        ArrayVector<Node> raw_skeleton(pathFinder.discoveryOrder());
        // from periphery to center: create skeleton tree and compute salience
        for(int k=raw_skeleton.size()-1; k >= 0; --k)
        {
            Node p1 = raw_skeleton[k],
                 p2 = pathFinder.predecessors()[p1];
            SNode & n1 = skeleton[p1];
            SNode & n2 = skeleton[p2];
            n1.parent = p2;
            

            // remove non-skeleton edges (i.e. set weight = infiniteWeight)
            for (BackArcIt arc(g, p1); arc != lemon::INVALID; ++arc)
            {
                Node p = g.target(*arc);
                if(weights[*arc] == infiniteWeight)
                    continue; // edge never was in the graph
                if(p == p2)
                    continue; // edge belongs to the skeleton
                if(pathFinder.predecessors()[p] == p1)
                    continue; // edge belongs to the skeleton
                if(n1.principal_child == lemon::INVALID || 
                   skeleton[p].principal_child == lemon::INVALID)
                    continue; // edge may belong to a loop => test later
                weights[*arc] = infiniteWeight;
            }

            // propagate length to parent if this is the longest subtree
            WeightType l = n1.length + norm(p1-p2);
            if(n2.length < l)
            {
                n2.length = l;
                n2.principal_child = p1;
            }
            
            if(compute_salience)
            {
                const double min_length = 4.0; // salience is meaningless for shorter segments due
                                               // to quantization noise (staircasing) of the boundary
                if(n1.length >= min_length)
                {
                    n1.salience = max(n1.salience, (n1.length + 0.5) / sqrt(squared_distance[p1]));
                
                    // propagate salience to parent if this is the most salient subtree
                    if(n2.salience < n1.salience)
                        n2.salience = n1.salience;
                }
            }
            else if(options.mode == SkeletonOptions::DontPrune)
                n1.salience = dest[p1];
            else
                n1.salience = n1.length;
        }
        
        // from center to periphery: propagate salience and compute twice the partial area
        for(int k=0; k < (int)raw_skeleton.size(); ++k)
        {
            Node p1 = raw_skeleton[k];
            SNode & n1 = skeleton[p1];
            Node p2 = n1.parent;
            SNode & n2 = skeleton[p2];
            
            if(p1 == n2.principal_child)
            {
                n1.length = n2.length;
                n1.salience = n2.salience;
            }
            else
            {
                n1.length += norm(p2-p1);
            }
            n1.partial_area = n2.partial_area + (p1[0]*p2[1] - p1[1]*p2[0]);
        }

        // always treat eccentricity center as a loop, so that it cannot be pruned 
        // away unless (options.mode & PreserveTopology) is false.
        skeleton[center].is_loop = true;
        
        // from periphery to center: * find and propagate loops
        //                           * delete branches not reaching the boundary
        detail::CheckForHole<std::size_t, MultiArrayView<2, T1, S1> > hasNoHole(label, labels);
        int hole_count = 0;
        double total_length = 0.0;
        for(int k=raw_skeleton.size()-1; k >= 0; --k)
        {
            Node p1 = raw_skeleton[k];
            SNode & n1 = skeleton[p1];

            if(n1.principal_child == lemon::INVALID)
            {
                for (ArcIt arc(g, p1); arc != lemon::INVALID; ++arc)
                {
                    Node p2 = g.target(*arc);
                    SNode * n2 = &skeleton[p2];
                    
                    if(n1.parent == p2)
                        continue; // going back to the parent can't result in a loop
                    if(weights[*arc] == infiniteWeight)
                        continue; // p2 is not in the tree or the loop has already been handled
                    // compute twice the area exclosed by the potential loop
                    MultiArrayIndex area2 = abs(n1.partial_area - (p1[0]*p2[1] - p1[1]*p2[0]) - n2->partial_area);
                    if(area2 <= 3) // area is too small to enclose a hole => loop is a discretization artifact
                        continue;
                    
                    // use Dijkstra to find the loop
                    WeightType edge_length = weights[*arc];
                    weights[*arc] = infiniteWeight;
                    pathFinder.reRun(weights, p1, p2);
                    Polygon<Shape2> poly;
                    {
                        poly.push_back_unsafe(p1);
                        poly.push_back_unsafe(p2);
                        Node p = p2;
                        do 
                        {
                            p = pathFinder.predecessors()[p];
                            poly.push_back_unsafe(p);
                        }
                        while(p != pathFinder.predecessors()[p]);
                    }
                    // check if the loop contains a hole or is just a discretization artifact
                    if(!inspectPolygon(poly, hasNoHole))
                    {
                        // it's a genuine loop => mark it and propagate salience
                        ++hole_count;
                        total_length += n1.length + n2->length;
                        double max_salience = max(n1.salience, n2->salience);
                        for(int p=1; p<poly.size(); ++p)
                        {
                            SNode & n = skeleton[poly[p]];
                            n.is_loop = true;
                            n.salience = max(n.salience, max_salience);
                        }
                    }
                }
                // delete skeleton branches that are not loops and don't reach the shape border
                // (these branches are discretization artifacts)
                if(!n1.is_loop && squared_distance[p1] >= 4)
                {
                    SNode * n = &n1;
                    while(true)
                    {
                        n->salience = 0;
                        // remove all of p1's edges
                        for(ArcIt arc(g, p1); arc != lemon::INVALID; ++arc)
                        {
                            weights[*arc] = infiniteWeight;
                        }
                        if(skeleton[n->parent].principal_child != p1)
                            break;
                        p1 = n->parent;
                        n = &skeleton[p1];
                    }
                }
            }
            
            if(n1.is_loop)
                skeleton[n1.parent].is_loop = true;
        }
        
        bool dont_prune = (options.mode & SkeletonOptions::Prune) == 0;
        bool preserve_topology = (options.mode & SkeletonOptions::PreserveTopology) != 0 ||
                                 options.mode == SkeletonOptions::Prune;
        bool relative_pruning = (options.mode & SkeletonOptions::Relative) != 0;
        WeightType threshold = (options.mode == SkeletonOptions::PruneTopology ||
                                options.mode == SkeletonOptions::Prune)
                                   ? infiniteWeight
                                   : relative_pruning
                                       ? options.pruning_threshold*skeleton[center].salience
                                       : options.pruning_threshold;
        // from center to periphery: write result
        int branch_count = 0;
        double average_length = 0;
        for(int k=0; k < (int)raw_skeleton.size(); ++k)
        {
            Node p1 = raw_skeleton[k];
            SNode & n1 = skeleton[p1];
            Node p2 = n1.parent;
            SNode & n2 = skeleton[p2];
            if(n1.principal_child == lemon::INVALID && 
               n1.salience >= threshold && 
               !n1.is_loop)
            {
                ++branch_count;
                average_length += n1.length;
                total_length += n1.length;
            }
            if(dont_prune)
                dest[p1] = n1.salience;
            else if(preserve_topology)
            {
                if(!n1.is_loop && n1.salience < threshold)
                    dest[p1] = 0;
            }
            else if(p1 != center && n1.salience < threshold)
                dest[p1] = 0;
        }
        if(branch_count > 0)
            average_length /= branch_count;
        
        if(features)
        {
            (*features)[label].diameter = center_line.length();
            (*features)[label].total_length = total_length;
            (*features)[label].average_length = average_length;
            (*features)[label].branch_count = branch_count;
            (*features)[label].hole_count = hole_count;
            (*features)[label].center = center;
            (*features)[label].terminal1 = center_line.front();
            (*features)[label].terminal2 = center_line.back();
            (*features)[label].euclidean_diameter = norm(center_line.back()-center_line.front());
        }
    }
        
    if(options.mode == SkeletonOptions::Prune)
        detail::skeletonThinning(squared_distance, dest, false);
}

class SkeletonFeatures
{
  public:
    double diameter, total_length, average_length, euclidean_diameter;
    UInt32 branch_count, hole_count;
    Shape2 center, terminal1, terminal2;
    
    SkeletonFeatures()
    : diameter(0)
    , total_length(0)
    , average_length(0)
    , euclidean_diameter(0)
    , branch_count(0)
    , hole_count(0)
    {}
};

/********************************************************/
/*                                                      */
/*                     skeletonizeImage                      */
/*                                                      */
/********************************************************/

    /*
    To compute the skeleton reliably in higher dimensions, we have to work on 
    a topological grid. The tricks to work with rounded skeletons on the 
    pixel grid probably don't generalize from 2D to 3D and higher. Specifically:
    
    * Compute Voronoi regions of the vector distance transformation according to
      identical support point to make sure that disconnected Voronoi regions
      still get only a single label.
    * Merge Voronoi regions whose support points are adjacent.
    * Mark skeleton candidates on the interpixel grid after the basic merge.
    * Detect skeleton segments simply by connected components labeling in the interpixel grid.
    * Skeleton segments form hyperplanes => use this property to compute segment 
      attributes.
    * Detect holes (and therefore, skeleton segments that are critical for topology)
      by computing the depth of each region/surface in the homotopy tree.
    * Add a pruning mode where holes are only preserved if their size exceeds a threshold. 
    
    To implement this cleanly, we first need a good implementation of the topological grid graph.
    */
// template <unsigned int N, class T1, class S1,
                          // class T2, class S2>
// void
// skeletonizeImage(MultiArrayView<N, T1, S1> const & labels,
            // MultiArrayView<N, T2, S2> dest,
            // SkeletonOptions const & options = SkeletonOptions())
// {

    /** \brief Skeletonization of all regions in a labeled 2D image.

        <b> Declarations:</b>

        \code
        namespace vigra {
            template <class T1, class S1,
                      class T2, class S2>
            void
            skeletonizeImage(MultiArrayView<2, T1, S1> const & labels,
                        MultiArrayView<2, T2, S2> dest,
                        SkeletonOptions const & options = SkeletonOptions());
        }
        \endcode

        This function computes the skeleton for each region in the 2D label image \a labels 
        and paints the results into the result image \a dest. Input label 
        <tt>0</tt> is interpreted as background and always ignored. Skeletons will be 
        marked with the same label as the corresponding region (unless options 
        <tt>returnLength()</tt> or <tt>returnSalience()</tt> are selected, see below). 
        Non-skeleton pixels will receive label <tt>0</tt>.

        For each region, the algorithm proceeds in the following steps:
        <ol>
        <li>Compute the \ref boundaryVectorDistance() relative to the \ref OuterBoundary of the region.</li>
        <li>Mark the raw skeleton: find 4-adjacent pixel pairs whose nearest boundary points are neither equal 
            nor adjacent and mark one pixel of the pair as a skeleton candidate. The resulting raw skeleton
            is 8-connected and thin. Skip the remaining steps when option <tt>dontPrune()</tt> is selected.</li>
        <li>Compute the eccentricity transform of the raw skeleton and turn the skeleton into a tree 
            whose root is the eccentricity center. When option <tt>pruneCenterLine()</tt> is selected,
            delete all skeleton points that do not belong to the two longest tree branches and 
            skip the remaining steps.</li>
        <li>For each pixel on the skeleton, compute its <tt>length</tt> attribute as the depth of the
            pixel's longest subtree. Compute its <tt>salience</tt> attribute as the ratio between 
            <tt>length</tt> and <tt>distance</tt>, where <tt>distance</tt> is the pixel's distance to 
            the nearest boundary point according to the distance transform. It holds that <tt>length >= 0.5</tt>
            and <tt>salience >= 1.0</tt>.</li>
        <li>Detect skeleton branching points and define <i>skeleton segments</i> as maximal connected pieces 
            without branching points.</li>
        <li>Compute <tt>length</tt> and <tt>salience</tt> of each segment as the maximum of these
            attributes among the pixels in the segment. When options <tt>returnLength()</tt> or 
            <tt>returnSalience()</tt> are selected, skip the remaining steps and return the 
            requested segment attribute in <tt>dest</tt>. In this case, <tt>dest</tt>'s 
            <tt>value_type</tt> should be a floating point type to exactly accomodate the 
            attribute values.</li>
        <li>Detect minimal cycles in the raw skeleton that enclose holes in the region (if any) and mark
            the corresponding pixels as critical for skeleton topology.</li>
        <li>Prune skeleton segments according to the selected pruning strategy and return the result. 
            The following pruning strategies are available:
            <ul>
            <li><tt>pruneLength(threshold, preserve_topology)</tt>: Retain only segments whose length attribute
                            exceeds the given <tt>threshold</tt>. When <tt>preserve_topology</tt> is true
                            (the defult), cycles around holes are preserved regardless of their length.
                            Otherwise, they are pruned as well.</li>
            <li><tt>pruneLengthRelative(threshold, preserve_topology)</tt>: Like <tt>pruneLength()</tt>,
                            but the threshold is specified as a fraction of the maximum segment length in  
                            the present region.</li>
            <li><tt>pruneSalience(threshold, preserve_topology)</tt>: Retain only segments whose salience attribute
                            exceeds the given <tt>threshold</tt>. When <tt>preserve_topology</tt> is true
                            (the defult), cycles around holes are preserved regardless of their salience.
                            Otherwise, they are pruned as well.</li>
            <li><tt>pruneSalienceRelative(threshold, preserve_topology)</tt>: Like <tt>pruneSalience()</tt>,
                            but the threshold is specified as a fraction of the maximum segment salience in 
                            the present region.</li>
            <li><tt>pruneTopology(preserve_center)</tt>: Retain only segments that are essential for the region's
                            topology. If <tt>preserve_center</tt> is true (the default), the eccentricity
                            center is also preserved, even if it is not essential. Otherwise, it might be 
                            removed. The eccentricity center is always the only remaining point when
                            the region has no holes.</li>
            </ul></li>
        </ol>
        
        The skeleton has the following properties:
        <ul>
        <li>It is 8-connected and thin (except when two independent branches happen to run alongside 
            before they divert). Skeleton points are defined by rounding the exact Euclidean skeleton
            locations to the nearest pixel.</li>
        <li>Skeleton branches terminate either at the region boundary or at a cycle. There are no branch 
            end points in the region interior.</li>
        <li>The salience threshold acts as a scale parameter: Large thresholds only retain skeleton 
            branches characterizing the general region shape. When the threshold gets smaller, ever
            more detailed boundary bulges will be represented by a skeleton branch.</li>
        </ul>
        
        Remark: If you have an application where a skeleton graph would be more useful
        than a skeleton image, function <tt>skeletonizeImage()</tt> can be changed/extended easily.

        <b> Usage:</b>

        <b>\#include</b> \<vigra/skeleton.hxx\><br/>
        Namespace: vigra

        \code
        Shape2 shape(width, height);
        MultiArray<2, UInt32> source(shape);
        MultiArray<2, UInt32> dest(shape);
        ...

        // Skeletonize and keep only those segments that are at least 10% of the maximum
        // length (the maximum length is half the skeleton diameter).
        skeletonizeImage(source, dest,
                    SkeletonOptions().pruneLengthRelative(0.1));
        \endcode

        \see vigra::boundaryVectorDistance()
    */
doxygen_overloaded_function(template <...> void skeletonizeImage)

template <class T1, class S1,
          class T2, class S2>
void
skeletonizeImage(MultiArrayView<2, T1, S1> const & labels,
            MultiArrayView<2, T2, S2> dest,
            SkeletonOptions const & options = SkeletonOptions())
{
    skeletonizeImageImpl(labels, dest, (ArrayVector<SkeletonFeatures>*)0, options);
}

template <class T, class S>
void
extractSkeletonFeatures(MultiArrayView<2, T, S> const & labels, 
                        ArrayVector<SkeletonFeatures> & features,
                        SkeletonOptions const & options = SkeletonOptions())
{
    MultiArray<2, float> skeleton(labels.shape());
    skeletonizeImageImpl(labels, skeleton, &features, options);
}

//@}

} //-- namespace vigra

#endif        //-- VIGRA_SKELETON_HXX
