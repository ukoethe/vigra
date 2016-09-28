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

/** \addtogroup GraphDataStructures
*/
//@{

/********************************************************/
/*                                                      */
/*                      BinaryForest                    */
/*                                                      */
/********************************************************/

/**
 * @brief BinaryForest stores a collection of rooted binary trees.
 *
 * Each connected component of the BinaryForest is thus a tree, and all edges are
 * directed away from the root node of the corresponding tree.
 *
 * @note If there is an arc from node <tt>u</tt> to node <tt>v</tt>, then the
 * arc ID is <tt>2*id(u)</tt> when <tt>v</tt> is the left child and <tt>2*id(u)+1</tt>
 * when <tt>v</tt> is the right child.
 */
class BinaryForest
{
public:

    typedef Int64 index_type;
    /// Node descriptor type of the present graph.
    typedef detail::NodeDescriptor<index_type> Node;
    /// Arc descriptor type of the present graph.
    typedef detail::ArcDescriptor<index_type> Arc;

    /// @brief Create an empty forest.
    BinaryForest();

    /// @brief Add a new node (its node ID will be selected automatically).
    Node addNode();
    /// @brief Add a new arc from node \a u to node \a v.
    /// The arc ID is <tt>2*id(u)</tt> if \a v is the left child of \a u, <tt>2*id(u)+1</tt> otherwise.
    Arc addArc(Node const & u, Node const & v);
    /// @brief Check if \a node exists.
    bool valid(Node const & node) const;
    /// @brief Check if \a arc exists.
    bool valid(Arc const & arc) const;
    /// @brief Find start node of \a arc.
    Node source(Arc const & arc) const;
    /// @brief Find end node of \a arc.
    Node target(Arc const & arc) const;
    /// @brief Get ID for node descriptor \a node.
    index_type id(Node const & node) const;
    /// @brief Get ID for arc descriptor \a arc.
    index_type id(Arc const & arc) const;
    /// @brief Get node descriptor for \a id.
    Node nodeFromId(index_type const & id) const;
    /// @brief Get arc descriptor for \a id.
    Arc arcFromId(index_type const & id) const;
    /// @brief Return the highest existing node ID.
    index_type maxNodeId() const;
    /// @brief Return the highest possible arc ID (equivalent to <tt>2*maxNodeId() + 1</tt>).
    index_type maxArcId() const;
    /// @brief Return the number of nodes (equivalent to <tt>maxNodeId()+1</tt>).
    size_t numNodes() const;
    /// @brief Return the number of arcs.
    ///        Always less than <tt>maxArcId()</tt> because not all arcs actually exist.
    size_t numArcs() const;
    /// @brief Return the number of incoming edges of \a node.
    ///        <tt>0</tt> for a root node, <tt>1</tt> otherwise.
    size_t inDegree(Node const & node) const;
    /// @brief Return the number of outgoing edges of \a node.
    ///        <tt>0</tt> for a leaf node, <tt>1</tt> or <tt>2</tt> otherwise.
    size_t outDegree(Node const & node) const;
    /// @brief Return the number of parents of \a node (equivalent to <tt>inDegree()</tt>).
    size_t numParents(Node const & node) const;
    /// @brief Return the number of children of \a node (equivalent to <tt>outDegree()</tt>).
    size_t numChildren(Node const & node) const;
    /// @brief Return the number of trees in the forest.
    size_t numRoots() const;
    /// @brief Create node cescriptor for ID \a i, or <tt>lemon::INVALID</tt> if
    ///        \a i is not a valid ID.
    Node getNode(size_t i) const;
    /// @brief Get the parent node descriptor of \a node, or <tt>lemon::INVALID</tt>
    ///        if \a node is a root or \a i is non-zero.
    Node getParent(Node const & node, size_t i = 0) const;
    /// @brief Get child number \a i of \a node.
    ///        Returns the left child if <tt>i=0</tt>, the right child if <tt>i=1</tt>,
    ///        and <tt>lemon::INVALID</tt> for other values of \a i or when the respective
    ///        is undefined.
    Node getChild(Node const & node, size_t i = 0) const;
    /// @brief Get the root node descriptor of tree \a i in the forest, or
    ///        <tt>lemon::INVALID</tt> if \a i is invalid.
    Node getRoot(size_t i = 0) const;
    /// @brief Merge two forests and increase the IDs of \a other to avoid ID clashes.
    ///        The function returns the offset that has been added to these IDs.
    size_t merge(BinaryForest const & other);

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
    // NOTE: The conversion to size_t is valid since we first check for >= 0.
    return node.id() >= 0 && (size_t)node.id() < nodes_.size();
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

inline BinaryForest::Node BinaryForest::getNode(
    size_t i
) const {
    if (i >= numNodes())
        return Node(lemon::INVALID);
    else
        return Node(i);
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

inline size_t BinaryForest::merge(
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
    return offset;
}

//@}

} // namespace vigra

#endif
