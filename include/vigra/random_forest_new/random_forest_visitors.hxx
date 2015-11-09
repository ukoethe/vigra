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
#ifndef VIGRA_RANDOM_FOREST_VISITORS_HXX
#define VIGRA_RANDOM_FOREST_VISITORS_HXX



namespace vigra
{



/**
 * @brief Base class from which all random forest visitors derive.
 * 
 * @details
 * Due to the parallel training, we cannot simply use a single visitor for all trees.
 * Instead, each tree gets a copy of the original visitor. The total procedure looks as follows:
 * - Call visit_at_beginning on the original visitor.
 * - Give a copy of the visitor to each tree.
 * - In each tree, call visit_before_tree (visit_after_tree) on the copy before (after) the tree was trained.
 * - Call visit_at_end (which gets a vector with the visitor copies) on the original visitor.
 */
class RFVisitorBase
{
public:

    RFVisitorBase()
    {}

    /**
     * @brief Do something before a tree has been learned.
     */
    void visit_before_tree()
    {}

    /**
     * @brief Do something after a tree has been learned.
     */
    void visit_after_tree()
    {}

    /**
     * @brief Do something before learning starts.
     */
    void visit_at_beginning()
    {}

    /**
     * @brief Do something after all trees have been learned.
     * 
     * @param v vector of size number_of_trees with the visitor copies
     */
    template <typename VISITORS>
    void visit_at_end(VISITORS & v)
    {}

};



class RFStopVisitor : public RFVisitorBase
{

};



} // namespace vigra

#endif
