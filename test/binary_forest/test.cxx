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
#include <vigra/binary_forest.hxx>
#include <vigra/unittest.hxx>

using namespace vigra;

struct BinaryForestTests
{
    typedef BinaryForest Graph;
    typedef Graph::Node Node;
    typedef Graph::Arc Arc;

    void test_basic_attributes()
    {
        Graph gr;
        Node n0 = gr.addNode();
        Node n1 = gr.addNode();
        Node n2 = gr.addNode();
        Node n3 = gr.addNode();
        Node n4 = gr.addNode();
        Node n5 = gr.addNode();
        Arc a01 = gr.addArc(n0, n1);
        Arc a02 = gr.addArc(n0, n2);
        Arc a13 = gr.addArc(n1, n3);
        Arc a14 = gr.addArc(n1, n4);
        Arc a25 = gr.addArc(n2, n5);

        should(n0 != n1);
        should(n0 != n2);
        should(n0 != n3);
        should(n0 != n4);
        should(n0 != n5);
        should(n1 != n2);
        should(n1 != n3);
        should(n1 != n4);
        should(n1 != n5);
        should(n2 != n3);
        should(n2 != n4);
        should(n2 != n5);
        should(n3 != n4);
        should(n3 != n5);
        should(n4 != n5);

        should(a01 != a02);
        should(a01 != a13);
        should(a01 != a14);
        should(a01 != a25);
        should(a02 != a13);
        should(a02 != a14);
        should(a02 != a25);
        should(a13 != a14);
        should(a13 != a25);
        should(a14 != a25);

        should(gr.numNodes() == 6);
        should(gr.numArcs() == 5);
        should(gr.maxNodeId() == 5);
        should(gr.maxArcId() == 11);

        should(gr.source(a01) == n0);
        should(gr.source(a02) == n0);
        should(gr.source(a13) == n1);
        should(gr.source(a14) == n1);
        should(gr.source(a25) == n2);

        should(gr.target(a01) == n1);
        should(gr.target(a02) == n2);
        should(gr.target(a13) == n3);
        should(gr.target(a14) == n4);
        should(gr.target(a25) == n5);

        should(gr.valid(n0));
        should(gr.valid(n1));
        should(gr.valid(n2));
        should(gr.valid(n3));
        should(gr.valid(n4));
        should(gr.valid(n5));
        should(!gr.valid(Node(lemon::INVALID)));
        should(!gr.valid(Node(123)));

        should(gr.valid(a01));
        should(gr.valid(a02));
        should(gr.valid(a13));
        should(gr.valid(a14));
        should(gr.valid(a25));
        should(!gr.valid(Arc(lemon::INVALID)));
        should(!gr.valid(Arc(123)));

        should(gr.inDegree(n0) == 0);
        should(gr.inDegree(n1) == 1);
        should(gr.inDegree(n2) == 1);
        should(gr.inDegree(n3) == 1);
        should(gr.inDegree(n4) == 1);
        should(gr.inDegree(n5) == 1);

        should(gr.outDegree(n0) == 2);
        should(gr.outDegree(n1) == 2);
        should(gr.outDegree(n2) == 1);
        should(gr.outDegree(n3) == 0);
        should(gr.outDegree(n4) == 0);
        should(gr.outDegree(n5) == 0);

        should(gr.numRoots() == 1);
        should(gr.getRoot() == n0);
        should(!gr.valid(gr.getRoot(1)));
    }

    void test_merge()
    {
        Graph gr;
        Node n0 = gr.addNode();
        Node n1 = gr.addNode();
        Node n2 = gr.addNode();
        Node n3 = gr.addNode();
        Arc a01 = gr.addArc(n0, n1);
        Arc a02 = gr.addArc(n0, n2);
        Arc a13 = gr.addArc(n1, n3);

        Graph gr2;
        {
            Node n0 = gr2.addNode();
            Node n1 = gr2.addNode();
            Node n2 = gr2.addNode();
            Node n3 = gr2.addNode();
            gr2.addArc(n0, n1);
            gr2.addArc(n1, n2);
            gr2.addArc(n1, n3);
            gr.merge(gr2);
        }

        should(gr.numNodes() == 8);
        should(gr.numArcs() == 6);
        should(gr.numRoots() == 2);

        // Get the new nodes.
        Node n4 = gr.nodeFromId(4);
        Node n5 = gr.nodeFromId(5);
        Node n6 = gr.nodeFromId(6);
        Node n7 = gr.nodeFromId(7);

        // Check the roots.
        should(gr.getRoot(0) == n0);
        should(gr.getRoot(1) == n4);

        // Check that the old nodes and arcs are still correct.
        should(gr.valid(n0));
        should(gr.valid(n1));
        should(gr.valid(n2));
        should(gr.valid(n3));
        should(gr.valid(a01));
        should(gr.valid(a02));
        should(gr.valid(a13));
        should(gr.getChild(n0, 0) == n1);
        should(gr.getChild(n0, 1) == n2);
        should(gr.getChild(n1, 0) == n3);
        should(gr.numChildren(n1) == 1);
        should(gr.numChildren(n2) == 0);
        should(gr.numChildren(n3) == 0);
        should(gr.numParents(n0) == 0);

        // Check the new nodes and arcs are still correct.
        should(gr.getChild(n4, 0) == n5);
        should(gr.getChild(n5, 0) == n6);
        should(gr.getChild(n5, 1) == n7);
        should(gr.numChildren(n4) == 1);
        should(gr.numChildren(n5) == 2);
        should(gr.numChildren(n6) == 0);
        should(gr.numChildren(n7) == 0);
        should(gr.numParents(n4) == 0);
    }

    template <ContainerTag CTag>
    void test_property_map()
    {
        PropertyMap<Node, int, CTag> m;
        Node n0(2);
        Node n1(5);
        Node n2(10);
        Node n3(27);
        m.insert(n0, 27);
        m.insert(n1, 12);
        m.insert(n2, 73);

        should(m.size() == 3);
        should(m.at(n0) == 27);
        should(m.at(n1) == 12);
        should(m.at(n2) == 73);
        should(m[n0] == 27);
        should(m[n1] == 12);
        should(m[n2] == 73);
        
        {
            auto it = m.find(n0);
            should(it != m.end());
            should(it->first == n0);
            should(it->second == 27);
            should(m.find(n3) == m.end());
        }

        {
            PropertyMap<Node, int, CTag> m2;
            m2 = m;
            should(m2.size() == 3);
            should(m[n0] == 27);
            should(m[n1] == 12);
            should(m[n2] == 73);
        }

        {
            std::vector<Node> keys, keys_expected;// = {n0, n1, n2};
            keys_expected.push_back(n0);
            keys_expected.push_back(n1);
            keys_expected.push_back(n2);
            std::vector<int> values, values_expected;// = {27, 12, 73};
            values_expected.push_back(27);
            values_expected.push_back(12);
            values_expected.push_back(73);
            for (auto const & p : m)
            {
                keys.push_back(p.first);
                values.push_back(p.second);
            }
            shouldEqualSequence(keys.begin(), keys.end(), keys_expected.begin());
            shouldEqualSequence(values.begin(), values.end(), values_expected.begin());
        }

        m.erase(n1);
        should(m.size() == 2);
        should(m.at(n0) == 27);
        should(m.at(n2) == 73);
        should(m[n0] == 27);
        should(m[n2] == 73);

        {
            std::vector<Node> keys, keys_expected;// = {n0, n2};
            keys_expected.push_back(n0);
            keys_expected.push_back(n2);
            std::vector<int> values, values_expected;// = {27, 73};
            values_expected.push_back(27);
            values_expected.push_back(73);
            for (auto const & p : m)
            {
                keys.push_back(p.first);
                values.push_back(p.second);
            }
            shouldEqualSequence(keys.begin(), keys.end(), keys_expected.begin());
            shouldEqualSequence(values.begin(), values.end(), values_expected.begin());
        }

        m.clear();
        should(m.size() == 0);
    }
};

struct BinaryForestTestSuite : public test_suite
{
    BinaryForestTestSuite()
        :
        test_suite("BinaryForest test")
    {
        add(testCase(&BinaryForestTests::test_basic_attributes));
        add(testCase(&BinaryForestTests::test_merge));
        add(testCase(&BinaryForestTests::test_property_map<MapTag>));
        add(testCase(&BinaryForestTests::test_property_map<IndexVectorTag>));
        add(testCase(&BinaryForestTests::test_property_map<VectorTag>));
    }
};

int main(int argc, char** argv)
{
    BinaryForestTestSuite forest_test;
    int failed = forest_test.run(testsToBeExecuted(argc, argv));
    std::cout << forest_test.report() << std::endl;
    return (failed != 0);
}
