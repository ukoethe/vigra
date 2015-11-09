/************************************************************************/
/*                                                                      */
/*        Copyright 2014-2015 by Ullrich Koethe and Philip Schill       */
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
#ifndef VIGRA_BINARY_FOREST_HXX
#define VIGRA_BINARY_FOREST_HXX

#include <vector>

#include "graphs.hxx"



namespace vigra
{



/**
 * @brief BinaryForest is binary directed acyclic graph, meaning a graph that consists of multiple binary trees.
 * 
 * @note If there is an arc from node u to node v, then the arc id is 2*id(u) if v is the left child and 2*id(u)+1 if v is the right child.
 */
class BinaryForest
{
public:

    typedef Int64 index_type;
    typedef detail::NodeDescriptor<index_type> Node;
    typedef detail::ArcDescriptor<index_type> Arc;

    BinaryForest();

    Node addNode();
    Arc addArc(Node const & u, Node const & v);
    bool valid(Node const & node) const;
    bool valid(Arc const & arc) const;
    Node source(Arc const & arc) const;
    Node target(Arc const & arc) const;
    index_type id(Node const & node) const;
    index_type id(Arc const & arc) const;
    Node nodeFromId(index_type const & id) const;
    Arc arcFromId(index_type const & id) const;
    index_type maxNodeId() const;
    index_type maxArcId() const;
    size_t numNodes() const;
    size_t numArcs() const;
    size_t inDegree(Node const & node) const;
    size_t outDegree(Node const & node) const;
    size_t numParents(Node const & node) const;
    size_t numChildren(Node const & node) const;
    size_t numRoots() const;
    Node getParent(Node const & node, size_t i = 0) const;
    Node getChild(Node const & node, size_t i = 0) const;
    Node getRoot(size_t i = 0) const;
    void merge(BinaryForest const & other);

private:

    struct NodeT
    {
        NodeT()
            :
            parent(-1),
            left_child(-1),
            right_child(-1)
        {}
        index_type parent, left_child, right_child;
    };

    std::vector<NodeT> nodes_;

    // Sorted vector with the node ids of the roots.
    std::vector<index_type> root_nodes_;

    size_t num_arcs_;
};



inline BinaryForest::BinaryForest()
    : 
    nodes_(),
    root_nodes_(),
    num_arcs_(0)
{}

inline BinaryForest::Node BinaryForest::addNode()
{
    Node n = Node(nodes_.size());
    nodes_.push_back(NodeT());
    root_nodes_.push_back(n.id());
    return n;
}

inline BinaryForest::Arc BinaryForest::addArc(
    Node const & u,
    Node const & v
){
    NodeT & u_node = nodes_[u.id()];
    NodeT & v_node = nodes_[v.id()];
    index_type arc_id = 2*u.id();

    // Make sure that the arc is not inserted twice.
    if (u_node.left_child == v.id())
        return Arc(arc_id);
    if (u_node.right_child == v.id())
        return Arc(arc_id+1);

    // Add v as child of u.
    if (u_node.left_child == -1)
    {
        u_node.left_child = v.id();
    }
    else if (u_node.right_child == -1)
    {
        u_node.right_child = v.id();
        ++arc_id;
    }
    else
    {
        vigra_fail("BinaryForest::addArc(): The node u already has two children.");
    }

    // Add u as parent of v.
    v_node.parent = u.id();

    // If v was a root node, remove it from the list.
    auto it = std::lower_bound(root_nodes_.begin(), root_nodes_.end(), v.id());
    if (it != root_nodes_.end() && !(v.id() < *it))
        root_nodes_.erase(it);

    ++num_arcs_;
    return Arc(arc_id);
}

inline bool BinaryForest::valid(
    Node const & node
) const {
    return node.id() >= 0 && node.id() < nodes_.size();
}

inline bool BinaryForest::valid(
    Arc const & arc
) const {
    if (arc == lemon::INVALID)
        return false;

    index_type const uid = arc.id()/2;
    if (!valid(Node(uid)))
        return false;

    if (arc.id() % 2 == 0)
        return nodes_[uid].left_child != -1;
    else
        return nodes_[uid].right_child != -1;
}

inline BinaryForest::Node BinaryForest::source(
    Arc const & arc
) const {
    return Node(arc.id()/2);
}

inline BinaryForest::Node BinaryForest::target(
    Arc const & arc
) const {
    NodeT const & u_node = nodes_[arc.id()/2];
    if (arc.id() % 2 == 0)
        return Node(u_node.left_child);
    else
        return Node(u_node.right_child);
}

inline BinaryForest::index_type BinaryForest::id(
    Node const & node
) const {
    return node.id();
}

inline BinaryForest::index_type BinaryForest::id(
    Arc const & arc
) const {
    return arc.id();
}

inline BinaryForest::Node BinaryForest::nodeFromId(
    index_type const & id
) const {
    return Node(id);
}

inline BinaryForest::Arc BinaryForest::arcFromId(
    index_type const & id
) const {
    return Arc(id);
}

inline BinaryForest::index_type BinaryForest::maxNodeId() const
{
    return nodes_.size()-1;
}

inline BinaryForest::index_type BinaryForest::maxArcId() const
{
    return 2*maxNodeId() + 1;
}

inline size_t BinaryForest::numNodes() const
{
    return nodes_.size();
}

inline size_t BinaryForest::numArcs() const
{
    return num_arcs_;
}

inline size_t BinaryForest::inDegree(
    Node const & node
) const {
    if (nodes_[node.id()].parent == -1)
        return 0;
    else
        return 1;
}

inline size_t BinaryForest::outDegree(
    Node const & node
) const {
    NodeT const & n = nodes_[node.id()];
    if (n.left_child == -1 && n.right_child == -1)
        return 0;
    else if (n.left_child == -1 || n.right_child == -1)
        return 1;
    else
        return 2;
}

inline size_t BinaryForest::numParents(
    Node const & node
) const {
    return inDegree(node);
}

inline size_t BinaryForest::numChildren(
    Node const & node
) const {
    return outDegree(node);
}

inline size_t BinaryForest::numRoots() const
{
    return root_nodes_.size();
}

inline BinaryForest::Node BinaryForest::getParent(
    Node const & node,
    size_t i
) const {
    NodeT const & n = nodes_[node.id()];
    if (n.parent == -1 || i != 0)
        return Node(lemon::INVALID);
    else
        return Node(n.parent);
}

inline BinaryForest::Node BinaryForest::getChild(
    Node const & node,
    size_t i
) const {
    NodeT const & n = nodes_[node.id()];
    if (i == 0)
        return Node(n.left_child);
    else if (i == 1)
        return Node(n.right_child);
    else
        return Node(lemon::INVALID);
}

inline BinaryForest::Node BinaryForest::getRoot(
    size_t i
) const {
    if (i >= root_nodes_.size())
        return Node(lemon::INVALID);
    else
        return Node(root_nodes_[i]);
}

inline void BinaryForest::merge(
    BinaryForest const & other
){
    num_arcs_ += other.num_arcs_;
    size_t const offset = nodes_.size();
    nodes_.insert(nodes_.end(), other.nodes_.begin(), other.nodes_.end());
    for (size_t i = offset; i < nodes_.size(); ++i)
    {
        NodeT & n = nodes_[i];
        if (n.parent != -1)
            n.parent += offset;
        if (n.left_child != -1)
            n.left_child += offset;
        if (n.right_child != -1)
            n.right_child += offset;
    }

    size_t const root_offset = root_nodes_.size();
    root_nodes_.insert(root_nodes_.end(), other.root_nodes_.begin(), other.root_nodes_.end());
    for (size_t i = root_offset; i < root_nodes_.size(); ++i)
    {
        root_nodes_[i] += offset;
    }
}



} // namespace vigra



#endif
