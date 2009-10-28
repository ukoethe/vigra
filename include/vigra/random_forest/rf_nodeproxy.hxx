/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by  Ullrich Koethe and Rahul Nair         */
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

#ifndef VIGRA_RANDOM_FOREST_NP_HXX
#define VIGRA_RANDOM_FOREST_NP_HXX

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


namespace vigra
{



enum NodeTags
{
    UnFilledNode        = 42,
    AllColumns          = 0x00000000,
    ToBePrunedTag       = 0x80000000,
    LeafNodeTag         = 0x40000000,

    i_ThresholdNode     = 0,
    i_HyperplaneNode    = 1,
    i_HypersphereNode   = 2,
    e_ConstProbNode     = 0 | LeafNodeTag,
    e_LogRegProbNode    = 1 | LeafNodeTag,
};

/** NodeBase class.

    \ingroup DecicionTree

    This class implements common features of all nodes.
    Memory Structure:
        Int32   Array:  TypeID, ParameterAddr, Child0, Child1, [ColumnData]0_
        double  Array:  NodeWeight, [Parameters]1_

*/


class NodeBase
{
  public:
    typedef Int32                               INT;
    typedef ArrayVector<INT>                    T_Container_type;
    typedef ArrayVector<double>                 P_Container_type;
    typedef T_Container_type::iterator     Topology_type;
    typedef P_Container_type::iterator     Parameter_type;


    mutable Topology_type                   topology_;
    size_t                                  topology_size_;

    mutable Parameter_type                  parameters_;
    size_t                                  parameter_size_ ;

    /** if numColumns = 0 then xrange is used as split axis
    **/
    static T_Container_type                 xrange;

        // Tree Parameters
    size_t                                  featureCount_;
    size_t                                  classCount_;

        // Node Parameters
    bool                                    hasData_;




        /**     Getters
        **/

    double &      weights()
    {
            return parameters_begin()[0];
    }

    bool          data()
    {
        return hasData_;
    }

    INT&          typeID()
    {
        return topology_[0];
    }

    INT&          parameter_addr()
    {
        return topology_[1];
    }

    /** Column Range **/
    Topology_type  column_data()
    {
        return topology_ + 4 ;
    }

    Topology_type columns_begin()
    {
        if(*column_data() == AllColumns)
            return NodeBase::xrange.begin();
        else
            return column_data()+1;
    }

    size_t      columns_size()
    {
        if(*column_data() == AllColumns)
            return featureCount_;
        else
            return *column_data();;
    }

    Topology_type  columns_end()
    {
        return columns_begin() + columns_size();
    }

    /** Topology Range **/
    Topology_type &  topology_begin()
    {
        return topology_;
    }
    Topology_type   topology_end()
    {
        return topology_begin() + topology_size();
    }
    size_t          topology_size()
    {
        return topology_size_;
    }

    /** Parameter Range **/
    Parameter_type  parameters_begin()
    {
        return parameters_;
    }
    Parameter_type  parameters_end()
    {
        return parameters_begin() + parameters_size();
    }

    size_t          parameters_size()
    {
        return parameter_size_;
    }

    INT &           child(Int32 l)
    {
        return topology_begin()[2+l];
    }

        /** Default Constructor**/
    NodeBase()
    :
                    hasData_(false)
    {}

        /** create ReadOnly Node at position n
        **/
    //TODO: switch parameter ordering - to spare a constructor.
    NodeBase(   T_Container_type    &  topology,
                P_Container_type    &  parameter,
                INT                         n)
    :
                    topology_   (topology.begin()+ n),
                    topology_size_(4),
                    parameters_  (parameter.begin() + parameter_addr()),
                    parameter_size_(1),
                    featureCount_(topology[0]),
                    classCount_(topology[1]),
                    hasData_(true)
    {
        while(xrange.size() <  featureCount_)
            xrange.push_back(xrange.size());
    }

    NodeBase(   size_t                      tLen,
                size_t                      pLen,
                T_Container_type    &  topology,
                P_Container_type    &  parameter,
                INT                         n)
    :
                    topology_   (topology.begin()+ n),
                    topology_size_(tLen),
                    parameters_  (parameter.begin() + parameter_addr()),
                    parameter_size_(pLen),
                    featureCount_(topology[0]),
                    classCount_(topology[1]),
                    hasData_(true)
    {
        while(xrange.size() <  featureCount_)
            xrange.push_back(xrange.size());
    }
        /** create new Node at end of vector+
        **/
    NodeBase(   size_t                      tLen,
                size_t                      pLen,
                T_Container_type   &        topology,
                P_Container_type   &        parameter)
    :
                    topology_size_(tLen),
                    parameter_size_(pLen),
                    featureCount_(topology[0]),
                    classCount_(topology[1]),
                    hasData_(true)
    {
        while(xrange.size() <  featureCount_)
            xrange.push_back(xrange.size());

        size_t n = topology.size();
        for(size_t ii = 0; ii < tLen; ++ii)
            topology.push_back(0);
        //topology.resize (n  + tLen);

        topology_           =   topology.begin()+ n;
        typeID()            =   UnFilledNode;

        parameter_addr()    =   parameter.size();

        //parameter.resize(parameter.size() + pLen);
        for(size_t ii = 0; ii < pLen; ++ii)
            parameter.push_back(0);

        parameters_          =   parameter.begin()+ parameter_addr();
        weights() = 1;
    }


        /** PseudoCopy Constructor  - Since each Node views on different data
            there can't be a real copy constructor (unless both objects should
            point to the same underlying data.                                  **/
    NodeBase(   NodeBase              &    toCopy,
                T_Container_type      &    topology,
                P_Container_type     &    parameter)
    :
                    topology_size_(toCopy.topology_size()),
                    parameter_size_(toCopy.parameters_size()),
                    featureCount_(topology[0]),
                    classCount_(topology[1]),
                    hasData_(true)
    {
        while(xrange.size() <  featureCount_)
            xrange.push_back(xrange.size());

        size_t n            = topology.size();
        for(size_t ii = 0; ii < toCopy.topology_size(); ++ii)
            topology.push_back(toCopy.topology_begin()[ii]);
//        topology.insert(topology.end(), toCopy.topology_begin(), toCopy.topology_end());
        topology_           =   topology.begin()+ n;
        parameter_addr()    =   parameter.size();
        for(size_t ii = 0; ii < toCopy.parameters_size(); ++ii)
            parameter.push_back(toCopy.parameters_begin()[ii]);
//        parameter.insert(parameter.end(), toCopy.parameters_begin(), toCopy.parameters_end());
        parameters_          =   parameter.begin()+ parameter_addr();
    }
};

 NodeBase::T_Container_type NodeBase::xrange;



template<NodeTags NodeType>
class Node;

template<>
class Node<i_ThresholdNode>
: public NodeBase
{


    public:
    typedef NodeBase BT;

        /**constructors **/

    Node(   BT::T_Container_type &   topology,
            BT::P_Container_type &   param)
                :   BT(5,2,topology, param)
    {
        BT::typeID() = i_ThresholdNode;
    }

    Node(   BT::T_Container_type     &   topology,
                    BT::P_Container_type     &   param,
                    INT                   n             )
                :   BT(5,2,topology, param, n)
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
};


template<>
class Node<i_HyperplaneNode>
: public NodeBase
{
    public:

    typedef NodeBase BT;

        /**constructors **/

    Node(           size_t                      nCol,
                    BT::T_Container_type    &   topology,
                    BT::P_Container_type    &   split_param)
                :   BT(nCol + 5,nCol + 2,topology, split_param)
    {
        BT::typeID() = i_HyperplaneNode;
    }

    Node(           BT::T_Container_type    &   topology,
                    BT::P_Container_type    &   split_param,
                    size_t                  n             )
                :   NodeBase(5 , 2,topology, split_param, n)
    {
        //TODO : is there a more elegant way to do this?
        BT::topology_size_ += BT::column_data()[0]== AllColumns ?
                                        0
                                    :   BT::column_data()[0];
        BT::parameter_size_ += BT::columns_size();
    }


    double& intercept()
    {
        return BT::parameters_begin()[1];
    }


    BT::Parameter_type weights()
    {
        return BT::parameters_begin()+2;
    }


    template<class U, class C>
    BT::INT next(MultiArrayView<2,U,C> const & feature)
    {
        double result = -1 * intercept();
        for(size_t ii = 0; ii < BT::columns_size(); ++ii)
        {
            result +=feature[BT::columns_begin()[ii]] * weights()[ii];
        }
        return result < 0 ? BT::child(0)
                          : BT::child(1);
    }
};



template<>
class Node<i_HypersphereNode>
: public NodeBase
{
    public:

    typedef NodeBase BT;

        /**constructors **/

    Node(           size_t                      nCol,
                    BT::T_Container_type    &   topology,
                    BT::P_Container_type    &   param)
                :   NodeBase(nCol + 5,nCol + 1,topology, param)
    {
        BT::typeID() = i_HypersphereNode;
    }

    Node(           BT::T_Container_type    &   topology,
                    BT::P_Container_type    &   param,
                    size_t                  n             )
                :   NodeBase(5, 1,topology, param, n)
    {
        BT::topology_size_ += BT::column_data()[0]== AllColumns ?
                                        0
                                    :   BT::column_data()[0];
        BT::parameter_size_ += BT::columns_size();
    }

    double& squaredRadius()
    {
        return BT::parameters_begin()[1];
    }


    BT::Parameter_type center()
    {
        return BT::parameters_begin()+2;
    }

    template<class U, class C>
    BT::INT next(MultiArrayView<2,U,C> const & feature)
    {
        double result = -1 * squaredRadius();
        for(size_t ii = 0; ii < BT::columns_size(); ++ii)
        {
            result += (feature[BT::columns_begin()[ii]] - center()[ii])*
                      (feature[BT::columns_begin()[ii]] - center()[ii]);
        }
        return result < 0 ? BT::child(0)
                          : BT::child(1);
    }
};


/** ExteriorNodeBase class.

    \ingroup DecicionTree

    This class implements common features of all interior nodes.
    All interior nodes are derived classes of ExteriorNodeBase.
*/






template<>
class Node<e_ConstProbNode>
: public NodeBase
{
    public:

    typedef     NodeBase    BT;

    Node(           BT::T_Container_type    &   topology,
                    BT::P_Container_type    &   param)
                    :
                BT(2,topology[1]+1, topology, param)

    {
        BT::typeID() = e_ConstProbNode;
    }


    Node(           BT::T_Container_type    &   topology,
                    BT::P_Container_type    &   param,
                    size_t                  n             )
                :   BT(2, topology[1]+1,topology, param, n)
    { }
    BT::Parameter_type prob_begin()
    {
        return BT::parameters_begin()+1;
    }
    BT::Parameter_type prob_end()
    {
        return prob_begin() + prob_size();
    }
    size_t prob_size()
    {
        return BT::classCount_;
    }
};

template<>
class Node<e_LogRegProbNode>;

} // namespace vigra

#endif //RF_nodeproxy
