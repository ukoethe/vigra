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

        TODO: Throw away the crappy iterators and use vigra::ArrayVectorView
             it is not like anybody else is going to use this NodeBase class
             is it?

        TODO: use the RF_Traits::ProblemSpec_t to specify the external 
             parameters instead of the options.
*/


class NodeBase
{
  public:
    typedef Int32                               INT;
    typedef ArrayVector<INT>                    T_Container_type;
    typedef ArrayVector<double>                 P_Container_type;
    typedef T_Container_type::iterator          Topology_type;
    typedef P_Container_type::iterator          Parameter_type;


    mutable Topology_type                       topology_;
    int                                         topology_size_;

    mutable Parameter_type                      parameters_;
    int                                         parameter_size_ ;

    /** if numColumns = 0 then xrange is used as split axis
    **/
    static T_Container_type                     xrange;

        // Tree Parameters
    int                                         featureCount_;
    int                                         classCount_;

        // Node Parameters
    bool                                        hasData_;




    /** get Node Weight
     */
    double &      weights()
    {
            return parameters_begin()[0];
    }

    double const &      weights() const
    {
            return parameters_begin()[0];
    }

    /** has the data been set?
     * todo: throw this out - bad design
     */
    bool          data() const
    {
        return hasData_;
    }

    /** get the node type id
     * \sa NodeTags
     */
    INT&          typeID()
    {
        return topology_[0];
    }

    INT const &          typeID() const
    {
        return topology_[0];
    }

    /** Where in the parameter_ array are the weights?
     */
    INT &          parameter_addr()
    {
        return topology_[1];
    }

    INT const &    parameter_addr() const
    {
        return topology_[1];
    }

    /** Column Range **/
    Topology_type  column_data() const
    {
        return topology_ + 4 ;
    }

    /** get the start iterator to the columns
     *  - once again - throw out - static members are crap.
     */
    Topology_type columns_begin() const
    {
            return column_data()+1;
    }

    /** how many columns?
     */
    int      columns_size() const
    {
        if(*column_data() == AllColumns)
            return featureCount_;
        else
            return *column_data();;
    }

    /** end iterator to the columns
     */
    Topology_type  columns_end() const
    {
        return columns_begin() + columns_size();
    }

    /** Topology Range - gives access to the raw Topo memory
     * the size_ member was added as a result of premature 
     * optimisation.
     */ 
    Topology_type   topology_begin() const
    {
        return topology_;
    }
    Topology_type   topology_end() const
    {
        return topology_begin() + topology_size();
    }
    int          topology_size() const
    {
        return topology_size_;
    }

    /** Parameter Range **/
    Parameter_type  parameters_begin() const
    {
        return parameters_;
    }
    Parameter_type  parameters_end() const
    {
        return parameters_begin() + parameters_size();
    }

    int          parameters_size() const
    {
        return parameter_size_;
    }


    /** where are the child nodes?
     */
    INT &           child(Int32 l)
    {
        return topology_begin()[2+l];
    }

    /** where are the child nodes?
     */
    INT const  &           child(Int32 l) const
    {
        return topology_begin()[2+l];
    }

    /** Default Constructor**/
    NodeBase()
    :
                    hasData_(false)
    {}

    /** create ReadOnly Base Node at position n (actual length is unknown)
     * only common features i.e. children etc are accessible.
     */
    NodeBase(   T_Container_type const   &  topology,
                P_Container_type const   &  parameter,
                INT                         n)
    :
                    topology_   (const_cast<Topology_type>(topology.begin()+ n)),
                    topology_size_(4),
                    parameters_  (const_cast<Parameter_type>(parameter.begin() + parameter_addr())),
                    parameter_size_(1),
                    featureCount_(topology[0]),
                    classCount_(topology[1]),
                    hasData_(true)
    {
        /*while((int)xrange.size() <  featureCount_)
            xrange.push_back(xrange.size());*/
    }

    /** create ReadOnly node with known length (the parameter range is valid)
     */
    NodeBase(   int                      tLen,
                int                      pLen,
                T_Container_type const & topology,
                P_Container_type const & parameter,
                INT                         n)
    :
                    topology_   (const_cast<Topology_type>(topology.begin()+ n)),
                    topology_size_(tLen),
                    parameters_  (const_cast<Parameter_type>(parameter.begin() + parameter_addr())),
                    parameter_size_(pLen),
                    featureCount_(topology[0]),
                    classCount_(topology[1]),
                    hasData_(true)
    {
        /*while((int)xrange.size() <  featureCount_)
            xrange.push_back(xrange.size());*/
    }
    /** create ReadOnly node with known length 
     * from existing Node
     */
    NodeBase(   int                      tLen,
                int                      pLen,
                NodeBase &               node)
    :
                    topology_   (node.topology_),
                    topology_size_(tLen),
                    parameters_  (node.parameters_),
                    parameter_size_(pLen),
                    featureCount_(node.featureCount_),
                    classCount_(node.classCount_),
                    hasData_(true)
    {
        /*while((int)xrange.size() <  featureCount_)
            xrange.push_back(xrange.size());*/
    }


   /** create new Node at end of vector
    * \param tLen number of integers needed in the topolog vector
    * \param plen number of parameters needed (this includes the node
    *           weight)*/
    NodeBase(   int                      tLen,
                int                      pLen,
                T_Container_type   &        topology,
                P_Container_type   &        parameter)
    :
                    topology_size_(tLen),
                    parameter_size_(pLen),
                    featureCount_(topology[0]),
                    classCount_(topology[1]),
                    hasData_(true)
    {
        /*while((int)xrange.size() <  featureCount_)
            xrange.push_back(xrange.size());*/

        int n = topology.size();
        for(int ii = 0; ii < tLen; ++ii)
            topology.push_back(0);
        //topology.resize (n  + tLen);

        topology_           =   topology.begin()+ n;
        typeID()            =   UnFilledNode;

        parameter_addr()    =   parameter.size();

        //parameter.resize(parameter.size() + pLen);
        for(int ii = 0; ii < pLen; ++ii)
            parameter.push_back(0);

        parameters_          =   parameter.begin()+ parameter_addr();
        weights() = 1;
    }


  /** PseudoCopy Constructor  - 
   *
   * Copy Node to the end of a container. 
   * Since each Node views on different data there can't be a real 
   * copy constructor (unless both objects should point to the 
   * same underlying data.                                  
   */
    NodeBase(   NodeBase      const  &    toCopy,
                T_Container_type      &    topology,
                P_Container_type     &    parameter)
    :
                    topology_size_(toCopy.topology_size()),
                    parameter_size_(toCopy.parameters_size()),
                    featureCount_(topology[0]),
                    classCount_(topology[1]),
                    hasData_(true)
    {
        /*while((int)xrange.size() <  featureCount_)
            xrange.push_back(xrange.size());*/

        int n            = topology.size();
        for(int ii = 0; ii < toCopy.topology_size(); ++ii)
            topology.push_back(toCopy.topology_begin()[ii]);
//        topology.insert(topology.end(), toCopy.topology_begin(), toCopy.topology_end());
        topology_           =   topology.begin()+ n;
        parameter_addr()    =   parameter.size();
        for(int ii = 0; ii < toCopy.parameters_size(); ++ii)
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

    Node(   BT::T_Container_type const     &   topology,
            BT::P_Container_type const     &   param,
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

    double const & threshold() const
    {
        return BT::parameters_begin()[1];
    }

    BT::INT& column()
    {
        return BT::column_data()[0];
    }
    BT::INT const & column() const
    {
        return BT::column_data()[0];
    }

    template<class U, class C>
    BT::INT  next(MultiArrayView<2,U,C> const & feature) const
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

    Node(           int                      nCol,
                    BT::T_Container_type    &   topology,
                    BT::P_Container_type    &   split_param)
                :   BT(nCol + 5,nCol + 2,topology, split_param)
    {
        BT::typeID() = i_HyperplaneNode;
    }

    Node(           BT::T_Container_type  const  &   topology,
                    BT::P_Container_type  const  &   split_param,
                    int                  n             )
                :   NodeBase(5 , 2,topology, split_param, n)
    {
        //TODO : is there a more elegant way to do this?
        BT::topology_size_ += BT::column_data()[0]== AllColumns ?
                                        0
                                    :   BT::column_data()[0];
        BT::parameter_size_ += BT::columns_size();
    }

    Node( BT & node_)
        :   BT(5, 2, node_) 
    {
        //TODO : is there a more elegant way to do this?
        BT::topology_size_ += BT::column_data()[0]== AllColumns ?
                                        0
                                    :   BT::column_data()[0];
        BT::parameter_size_ += BT::columns_size();
    }


    double const & intercept() const
    {
        return BT::parameters_begin()[1];
    }
    double& intercept()
    {
        return BT::parameters_begin()[1];
    }

    BT::Parameter_type weights() const
    {
        return BT::parameters_begin()+2;
    }

    BT::Parameter_type weights()
    {
        return BT::parameters_begin()+2;
    }


    template<class U, class C>
    BT::INT next(MultiArrayView<2,U,C> const & feature) const
    {
        double result = -1 * intercept();
        if(*(BT::column_data()) == AllColumns)
        {
            for(int ii = 0; ii < BT::columns_size(); ++ii)
            {
                result +=feature[ii] * weights()[ii];
            }
        }
        else
        {
            for(int ii = 0; ii < BT::columns_size(); ++ii)
            {
                result +=feature[BT::columns_begin()[ii]] * weights()[ii];
            }
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

    Node(           int                      nCol,
                    BT::T_Container_type    &   topology,
                    BT::P_Container_type    &   param)
                :   NodeBase(nCol + 5,nCol + 1,topology, param)
    {
        BT::typeID() = i_HypersphereNode;
    }

    Node(           BT::T_Container_type  const  &   topology,
                    BT::P_Container_type  const  &  param,
                    int                  n             )
                :   NodeBase(5, 1,topology, param, n)
    {
        BT::topology_size_ += BT::column_data()[0]== AllColumns ?
                                        0
                                    :   BT::column_data()[0];
        BT::parameter_size_ += BT::columns_size();
    }

    Node( BT & node_)
        :   BT(5, 1, node_) 
    {
        BT::topology_size_ += BT::column_data()[0]== AllColumns ?
                                        0
                                    :   BT::column_data()[0];
        BT::parameter_size_ += BT::columns_size();

    }

    double const & squaredRadius() const
    {
        return BT::parameters_begin()[1];
    }

    double& squaredRadius()
    {
        return BT::parameters_begin()[1];
    }

    BT::Parameter_type center() const
    {
        return BT::parameters_begin()+2;
    }

    BT::Parameter_type center()
    {
        return BT::parameters_begin()+2;
    }

    template<class U, class C>
    BT::INT next(MultiArrayView<2,U,C> const & feature) const
    {
        double result = -1 * squaredRadius();
        if(*(BT::column_data()) == AllColumns)
        {
            for(int ii = 0; ii < BT::columns_size(); ++ii)
            {
                result += (feature[ii] - center()[ii])*
                          (feature[ii] - center()[ii]);
            }
        }
        else
        {
            for(int ii = 0; ii < BT::columns_size(); ++ii)
            {
                result += (feature[BT::columns_begin()[ii]] - center()[ii])*
                          (feature[BT::columns_begin()[ii]] - center()[ii]);
            }
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


    Node(           BT::T_Container_type const &   topology,
                    BT::P_Container_type const &   param,
                    int                  n             )
                :   BT(2, topology[1]+1,topology, param, n)
    { }


    Node( BT & node_)
        :   BT(2, node_.classCount_ +1, node_) 
    {}
    BT::Parameter_type  prob_begin() const
    {
        return BT::parameters_begin()+1;
    }
    BT::Parameter_type  prob_end() const
    {
        return prob_begin() + prob_size();
    }
    int prob_size() const
    {
        return BT::classCount_;
    }
};

template<>
class Node<e_LogRegProbNode>;

} // namespace vigra

#endif //RF_nodeproxy
