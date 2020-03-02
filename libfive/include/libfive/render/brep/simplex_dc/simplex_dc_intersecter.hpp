/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <atomic>
#include <array>
#include <vector>

#include "libfive/render/axes.hpp"
#include "libfive/eval/tape.hpp"
#include "libfive/tree/tree.hpp"

namespace libfive {

// Forward declarations
template <unsigned N, class Leaf> class SimplexTree; 
template <unsigned N> struct SimplexLeafSubspace;
template <unsigned N> struct SimplexDCLeaf;
template <unsigned N> struct SimplexDCMinEdge;
template <unsigned N> struct SimplexDCIntersection;
template <typename... T> class ObjectPool;
template <unsigned N> class PerThreadBRep;
template <unsigned N> class BRep;
class Evaluator;

/*  This dual-walker class finds the intersections between the surface and
 *  each edge of the simplices in the input tree, storing them in the
 *  appropriate simplices.  The returned output contains vertices only (no
 *  branes), and only those for the intersections (which serve to assist
 *  triangulation of the polygons that will eventually be built around them).*/

template <unsigned N>
class SimplexDCIntersecter
{
public:
    using PerThreadOutput = PerThreadBRep<N>;
    using Output = BRep<N>;
    using Input = SimplexTree<N, SimplexDCLeaf<N>>;

    using Perp = Eigen::Array<double, 3 - N, 1>;

    /*
     *  Constructs an intersecter that owns an evaluator,
     *  which is built from the given tree.
     */
    SimplexDCIntersecter(PerThreadOutput& m, Tree t,
                         ObjectPool<SimplexDCIntersection<N>>& pool,
                         Perp perp);

    /*
     *  Constructs an intersecter that has borrowed an evaluator,
     *  which is useful in cases where constructing evaluators
     *  is expensive and they should be re-used.
     */
    SimplexDCIntersecter(PerThreadOutput& m, Evaluator* es,
                         ObjectPool<SimplexDCIntersection<N>>& pool,
                         Perp perp);

    SimplexDCIntersecter(SimplexDCIntersecter&&) = default;

    ~SimplexDCIntersecter();

    /* Empty cell loader, called by Dual::walk */
    void load(Input* input) {}

    /*  Called by Dual::walk to intersect the surface with simplex edges
     *  whose lowest-dimensional vertex is a face vertex. Should not be
     *  called in 2d.*/
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 2)>& ts);

    /*  Called by Dual::walk to intersect the surface with simplex edges
     *  whose lowest-dimensional vertex is an edge vertex. */
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 1)>& ts);


    /*  Called by Dual::walk to intersect the surface with simplex edges
     *  whose lowest-dimensional vertex is a corner vertex. */
    void load(const std::array<Input*, 1 << N>& ts);


    /*
     *  Simplex DC meshing needs to walk the top edges of the tree,
     *  because those include tets that may need to be handled.
     */
    static bool needsTopEdges() { return true; }

protected:
    /*
     *  Performs a binary search along a particular edge using the
     *  provided tape, and returns a SimplexDCIntersection (allocated
     *  by pool) with refcount 0.
     */
    SimplexDCIntersection<N>* searchEdge(SimplexLeafSubspace<N>* inside,
                                         SimplexLeafSubspace<N>* outside,
                                         const std::shared_ptr<Tape>& tape);

    ObjectPool<SimplexDCIntersection<N>>& parent_pool;
    ObjectPool<SimplexDCIntersection<N>> my_pool;
    Evaluator* eval;
    bool owned;
    Perp perp;

    PerThreadBRep<N>& m;

    struct CmpPts { // For debugging, remove
        bool operator()(Eigen::Matrix<float, 3, 1> a, Eigen::Matrix<float, 3, 1> b) const {
            return std::array{ a.x(), a.y(), a.z() } < std::array{ b.x(), b.y(), b.z() };
        }
        bool operator()(Eigen::Matrix<float, 2, 1> a, Eigen::Matrix<float, 2, 1> b) const {
            return std::array{ a.x(), a.y() } < std::array{ b.x(), b.y() };
        }
    };
    std::map<Eigen::Matrix<float, N, 1>, std::pair<Eigen::Matrix<double, N, 1>, Eigen::Matrix<double, N, 1>>, CmpPts> sources;

};

////////////////////////////////////////////////////////////////////////////////

}   // namespace libfive
