/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <array>

#include <Eigen/Eigen>

#include "libfive/render/axes.hpp"
#include "libfive/eval/tape.hpp"
#include "libfive/tree/tree.hpp"

namespace libfive {

// Forward declarations
template <unsigned N, class Leaf> class SimplexTree;
template <unsigned N> struct SimplexLeaf;
template <unsigned N> class PerThreadBRep;
class Evaluator;
class Mesh;

class SimplexMesher
{
public:
    using Output = Mesh;
    using PerThreadOutput = PerThreadBRep<3>;
    using Input = const SimplexTree<3, SimplexLeaf<3>>;

    /*
     *  Constructs a mesher that owns an evaluator,
     *  which is built from the given tree.
     */
    SimplexMesher(PerThreadBRep<3>& m, Tree t);

    /*
     *  Constructs a mesher that has borrowed an evaluator,
     *  which is useful in cases where constructing evaluators
     *  is expensive and they should be re-used.
     */
    SimplexMesher(PerThreadBRep<3>& m, Evaluator* es);

    ~SimplexMesher();

    /* Empty cell loader, called by Dual::walk */
    void load(Input* input) {}

    /*
     *  Called by Dual::walk to construct the triangle mesh
     */
    template <Axis::Axis A>
    void load(const std::array<const Input*, 4>& ts);

    /*
     *  Empty face loader, called by Dual::walk
     */
    template <Axis::Axis A>
    void load(const std::array<const Input*, 2> & ts) {}

    /* Empty corner loader */
    void load(const std::array<const Input*, 8> & ts) {}

    /*
     *  Simplex meshing needs to walk the top edges of the tree,
     *  because those include tets that we have to run MT on.
     */
    static bool needsTopEdges() { return true; }

protected:
    /*
     *  Performs a binary search along a particular edge using the
     *  provided tape.  Stores the resulting vertex into the Mesh m,
     *  and returns its index.
     */
    uint64_t searchEdge(Eigen::Vector3d inside, Eigen::Vector3d outside,
                        const std::shared_ptr<Tape>& tape);

    PerThreadBRep<3>& m;
    Evaluator* eval;
    bool owned;
};

////////////////////////////////////////////////////////////////////////////////

}   // namespace libfive
