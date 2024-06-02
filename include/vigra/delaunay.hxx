/************************************************************************/
/*                                                                      */
/*               Copyright 2013-2017 by Benjamin Seppke                 */
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

/*
 * Parts of this file, namely the first section including the delaunay
 * function are taken from the LEMON distribution, lgf-gen.cc. Here, the
 * following copyright applies:
 *
 * Copyright (C) 2003-2009
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#ifndef VIGRA_DELAUNAY_HXX
#define VIGRA_DELAUNAY_HXX

#ifndef WITH_LEMON
    #error "Should only be included with flag \"WITH_LEMON\""
#endif

#include <set>
#include <lemon/list_graph.h>
#include <lemon/dim2.h>
#include <lemon/bfs.h>

namespace vigra {

    namespace detail
    {
        using namespace lemon;
        typedef dim2::Point<double> Point;

        GRAPH_TYPEDEFS(ListGraph);
        
        struct Part
        {
            int prev, curr, next;

            Part(int p, int c, int n) : prev(p), curr(c), next(n) {}
        };

        inline std::ostream& operator<<(std::ostream& os, const Part& part)
        {
            os << '(' << part.prev << ',' << part.curr << ',' << part.next << ')';
            return os;
        }

        inline double circle_point(const Point& p, const Point& q, const Point& r)
        {
            double a = p.x * (q.y - r.y) + q.x * (r.y - p.y) + r.x * (p.y - q.y);
            if (a == 0) return std::numeric_limits<double>::quiet_NaN();

            double d = (p.x * p.x + p.y * p.y) * (q.y - r.y) +
              (q.x * q.x + q.y * q.y) * (r.y - p.y) +
              (r.x * r.x + r.y * r.y) * (p.y - q.y);

            double e = (p.x * p.x + p.y * p.y) * (q.x - r.x) +
              (q.x * q.x + q.y * q.y) * (r.x - p.x) +
              (r.x * r.x + r.y * r.y) * (p.x - q.x);

            double f = (p.x * p.x + p.y * p.y) * (q.x * r.y - r.x * q.y) +
              (q.x * q.x + q.y * q.y) * (r.x * p.y - p.x * r.y) +
              (r.x * r.x + r.y * r.y) * (p.x * q.y - q.x * p.y);

            return d / (2 * a) + std::sqrt((d * d + e * e) / (4 * a * a) + f / a);
        }

        inline bool circle_form(const Point& p, const Point& q, const Point& r)
        {
            return rot90(q - p) * (r - q) < 0.0;
        }

        inline double intersection(const Point& p, const Point& q, double sx)
        {
            const double epsilon = 1e-8;

            if (p.x == q.x) return (p.y + q.y) / 2.0;

            if (sx < p.x + epsilon) return p.y;
            if (sx < q.x + epsilon) return q.y;

            double a = q.x - p.x;
            double b = (q.x - sx) * p.y - (p.x - sx) * q.y;
            double d = (q.x - sx) * (p.x - sx) * (p - q).normSquare();
            return (b - std::sqrt(d)) / a;
        }

        struct YLess
        {
            YLess(const std::vector<Point>& points, double& sweep)
              : _points(points),
                _sweep(sweep)
            {
            }

            bool operator()(const Part& l, const Part& r) const
            {
              const double epsilon = 1e-8;

              //      std::cerr << l << " vs " << r << std::endl;
              double lbx = l.prev != -1 ?
                intersection(_points[l.prev], _points[l.curr], _sweep) :
                - std::numeric_limits<double>::infinity();
              double rbx = r.prev != -1 ?
                intersection(_points[r.prev], _points[r.curr], _sweep) :
                - std::numeric_limits<double>::infinity();
              double lex = l.next != -1 ?
                intersection(_points[l.curr], _points[l.next], _sweep) :
                std::numeric_limits<double>::infinity();
              double rex = r.next != -1 ?
                intersection(_points[r.curr], _points[r.next], _sweep) :
                std::numeric_limits<double>::infinity();

              if (lbx > lex) std::swap(lbx, lex);
              if (rbx > rex) std::swap(rbx, rex);

              if (lex < epsilon + rex && lbx + epsilon < rex) return true;
              if (rex < epsilon + lex && rbx + epsilon < lex) return false;
              return lex < rex;
            }

          const std::vector<Point>& _points;
          double& _sweep;
      };

      struct BeachIt;

      typedef std::multimap<double, BeachIt> SpikeHeap;

      typedef std::multimap<Part, SpikeHeap::iterator, YLess> Beach;

        struct BeachIt
        {
            Beach::iterator it;
            BeachIt(Beach::iterator iter)
            : it(iter)
            {
            }
        };
    } //end of namespace vigra::detail

void delaunay(detail::ListGraph & g, detail::ListGraph::NodeMap<detail::Point>& coords)
{
    using namespace detail;

    typedef std::vector<std::pair<double, int> > SiteHeap;

    std::vector<Point> points;
    std::vector<Node> nodes;
    
    for (NodeIt it(g); it != INVALID; ++it)
    {
        nodes.push_back(it);
        points.push_back(coords[it]);
    }
    
    std::cerr << "Nodes: " << nodes.size() << "\n";
    
    SiteHeap siteheap(points.size());

    double sweep;

    for (int i = 0; i < int(siteheap.size()); ++i)
    {
      siteheap[i] = std::make_pair(points[i].x, i);
    }

    std::sort(siteheap.begin(), siteheap.end());
    sweep = siteheap.front().first;

    YLess yless(points, sweep);
    Beach beach(yless);

    SpikeHeap spikeheap;

    std::set<std::pair<int, int> > arcs;

    int siteindex = 0;
    {
        SiteHeap front;

        while (siteindex < int(siteheap.size()) &&
               siteheap[0].first == siteheap[siteindex].first)
        {
            front.push_back(std::make_pair(points[siteheap[siteindex].second].y,
                                           siteheap[siteindex].second));
            ++siteindex;
        }

        std::sort(front.begin(), front.end());

        for (int i = 0; i < int(front.size()); ++i)
        {
            int prev = (i == 0 ? -1 : front[i - 1].second);
            int curr = front[i].second;
            int next = (i + 1 == int(front.size()) ? -1 : front[i + 1].second);

            beach.insert(std::make_pair(Part(prev, curr, next),
                                        spikeheap.end()));
        }
    }

    while (siteindex < int(points.size()) || !spikeheap.empty())
    {

        SpikeHeap::iterator spit = spikeheap.begin();

        if ( siteindex < int(points.size()) &&
             (spit == spikeheap.end() || siteheap[siteindex].first < spit->first))
        {
            int site = siteheap[siteindex].second;
            sweep = siteheap[siteindex].first;
            
            Beach::iterator bit = beach.upper_bound(Part(site, site, site));

            if (bit->second != spikeheap.end())
            {
                spikeheap.erase(bit->second);
            }

            int prev = bit->first.prev;
            int curr = bit->first.curr;
            int next = bit->first.next;

            beach.erase(bit);

            SpikeHeap::iterator pit = spikeheap.end();
            if (prev != -1 &&
                circle_form(points[prev], points[curr], points[site]))
            {
                double x = circle_point(points[prev], points[curr], points[site]);
                pit = spikeheap.insert(std::make_pair(x, BeachIt(beach.end())));
                pit->second.it = beach.insert(std::make_pair(Part(prev, curr, site), pit));
            }
            else
            {
                beach.insert(std::make_pair(Part(prev, curr, site), pit));
            }

            beach.insert(std::make_pair(Part(curr, site, curr), spikeheap.end()));

            SpikeHeap::iterator nit = spikeheap.end();
          
            if (next != -1 &&
                circle_form(points[site], points[curr],points[next]))
            {
                double x   = circle_point(points[site], points[curr], points[next]);
                       nit = spikeheap.insert(std::make_pair(x, BeachIt(beach.end())));

                nit->second.it = beach.insert(std::make_pair(Part(site, curr, next), nit));
            }
            else
            {
                beach.insert(std::make_pair(Part(site, curr, next), nit));
            }
            ++siteindex;
        }
        else
        {
            sweep = spit->first;

            Beach::iterator bit = spit->second.it;

            int prev = bit->first.prev;
            int curr = bit->first.curr;
            int next = bit->first.next;

            
           // if(prev != -1 && next !=-1)
           //     std::cerr << "Triangle: (" << prev << ", " << curr << ", " << next << ")\n";
            {
                std::pair<int, int> arc;

                arc = prev < curr ? std::make_pair(prev, curr) : std::make_pair(curr, prev);

                if (arcs.find(arc) == arcs.end())
                {
                    arcs.insert(arc);
                    g.addEdge(nodes[prev], nodes[curr]);
                }

                arc = curr < next ? std::make_pair(curr, next) : std::make_pair(next, curr);

                if (arcs.find(arc) == arcs.end())
                {
                    arcs.insert(arc);
                    g.addEdge(nodes[curr], nodes[next]);
                }
            }

            Beach::iterator pbit = bit; --pbit;
            int ppv = pbit->first.prev;
            Beach::iterator nbit = bit; ++nbit;
            int nnt = nbit->first.next;

            if (bit->second != spikeheap.end())
            {
                spikeheap.erase(bit->second);
            }
            if (pbit->second != spikeheap.end())
            {
                spikeheap.erase(pbit->second);
            }
            if (nbit->second != spikeheap.end())
            {
                spikeheap.erase(nbit->second);
            }

            beach.erase(nbit);
            beach.erase(bit);
            beach.erase(pbit);

            SpikeHeap::iterator pit = spikeheap.end();
          
            if (ppv != -1 && ppv != next &&
              circle_form(points[ppv], points[prev], points[next]))
            {
              double x = circle_point(points[ppv], points[prev], points[next]);

              if (x < sweep)
              {
                  x = sweep;
              }
              pit = spikeheap.insert(std::make_pair(x, BeachIt(beach.end())));
              pit->second.it = beach.insert(std::make_pair(Part(ppv, prev, next), pit));
            }
            else
            {
              beach.insert(std::make_pair(Part(ppv, prev, next), pit));
            }

            SpikeHeap::iterator nit = spikeheap.end();
            if (nnt != -1 && prev != nnt &&
              circle_form(points[prev], points[next], points[nnt]))
            {
              double x = circle_point(points[prev], points[next], points[nnt]);
              if (x < sweep)
              {
                  x = sweep;
              }
              nit = spikeheap.insert(std::make_pair(x, BeachIt(beach.end())));
              nit->second.it = beach.insert(std::make_pair(Part(prev, next, nnt), nit));
            }
            else
            {
                beach.insert(std::make_pair(Part(prev, next, nnt), nit));
            }
        }
    }

    //Add missing convex hull arcs to complete triangulation
    for (Beach::iterator it = beach.begin(); it != beach.end(); ++it)
    {
        int curr = it->first.curr;
        int next = it->first.next;

        
        if (next == -1)
        {
            continue;
        }
        
        std::pair<int, int> arc;
        arc = curr < next ? std::make_pair(curr, next) : std::make_pair(next, curr);

        if (arcs.find(arc) == arcs.end())
        {
            arcs.insert(arc);
            g.addEdge(nodes[curr], nodes[next]);
        }
    }
}

std::vector<vigra::triple<int,int,int> > trianglesFromDelaunay(detail::ListGraph & g, detail::ListGraph::NodeMap<int>& nodeIDs)
{
    using namespace detail;
    
    std::vector<vigra::triple<int,int,int> > result;
    
    std::set<std::pair<int, int> > arcs;
    
    for(ArcIt aIt(g); aIt != INVALID; ++aIt)
    {
        int curr = nodeIDs[g.source(aIt)];
        int next = nodeIDs[g.target(aIt)];
        
        if (curr < next)
            arcs.insert(std::pair<int, int>(curr, next));
    }

    for(ArcIt aIt(g); aIt != INVALID; ++aIt)
    {
        int curr = nodeIDs[g.source(aIt)];
        int next = nodeIDs[g.target(aIt)];
        //Order: curr < next
        if (curr > next)
            continue;
        
        int found=0;
        
        //Find third node of the triangle
        for (NodeIt it(g); it != INVALID; ++it)
        {
            //Max. 2 triangles per edge
            if(found == 2)
                break;
            
            int cand = nodeIDs[it];
            
            //Order: curr < next < cand
            if (cand > next)
            {
                std::pair<int, int> arc1;
                std::pair<int, int> arc2;
                arc1 = curr < cand ? std::make_pair(curr, cand) : std::make_pair(cand, curr);
                arc2 = next < cand ? std::make_pair(next, cand) : std::make_pair(cand, next);
               
                if (arcs.find(arc1) != arcs.end() && arcs.find(arc2) != arcs.end())
                {
                    //std::cerr << "Triangle: (" << curr << ", " << next << ", " << cand <<")\n";
                    result.push_back(vigra::triple<int,int,int>(curr, next, cand));
                    found++;
                }
            }
        }
    }
    return result;
}

} //end of namespace vigra

#endif //VIGRA_DELAUNAY_HXX
