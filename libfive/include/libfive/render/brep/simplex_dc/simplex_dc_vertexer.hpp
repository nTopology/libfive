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

namespace libfive {

// Forward declarations
template <unsigned N, class Leaf> class SimplexTree;
template <unsigned N> struct SimplexDCLeaf;
template <unsigned N> class PerThreadBRep;
template <unsigned N> class BRep;
template <unsigned N> struct SimplexLeafSubspace;
template <unsigned N> struct DCSimplex;

template <unsigned N>
class SimplexDCVertexer
{
public:
    using Input = SimplexTree<N, SimplexDCLeaf<N>>;
    // The BRep output by the vertexer has vertices only, no branes.
    using Output = BRep<N>;
    using PerThreadOutput = PerThreadBRep<N>;
    /*
     *  Constructor.
     */
    SimplexDCVertexer(PerThreadBRep<N>& m, size_t vertsOffset) 
        : m(m), offset(vertsOffset) {}

    /* 
     *  Empty cell loader, called by Dual::walk 
     */
    void load(Input* input) {}

    /*  
     *  Empty face loader.  Should not be called in 2d.
     */
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 2)> & ts) { assert(N == 3); };

    /*  
     *  Edge loader; calculates and loads the vertices for all N!*2^N simplices
     *  adjoining this edge. 
     */
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 1)> & ts);


    /*  
     *  Empty corner walker.
     */
    void load(const std::array<Input*, 1 << N> & ts) {}


    /*
     *  Simplex DC meshing needs to walk the top edges of the tree,
     *  because those include tets that may need to be handled.
     */
    static bool needsTopEdges() { return true; }

protected:
    /*
     * Represents the vertices from subspaces forming a particular simplex.
     */

    using SubspaceVertArray = std::array<const SimplexLeafSubspace<N>*, N + 1>;
    /*
     *  Calculates and stores the vertex for a given simplex.
     */
    void calcAndStoreVert(
        DCSimplex<N>& simplex, SubspaceVertArray vertsFromSubspaces);

    PerThreadBRep<N>& m;

    size_t offset; // Offset when storing vertices, since early vertex numbers
                   // are taken by intersection vertices.
};

////////////////////////////////////////////////////////////////////////////////

}   // namespace libfive
