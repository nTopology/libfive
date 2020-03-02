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
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"

namespace libfive {

// Forward declarations
template <unsigned N, class Leaf> class SimplexTree;
template <unsigned N> struct SimplexDCLeaf;
template <unsigned N> class PerThreadBRep;
template <unsigned N> struct DCSimplex;
template <unsigned N, bool constTree>
class SimplexDCCirculator;
class Mesh;

/*  This dual-walker class takes a simplex DC tree that is complete (had its
 *  vertices set and assigned to indices and creates the appropriate branes for
 *  the resulting mesh.*/

class SimplexDCMesher
{
public:
    using Input = const SimplexTree<3, SimplexDCLeaf<3>>;
    // The BRep output by the vertexer has branes only, no vertices.
    using Output = Mesh;
    using PerThreadOutput = PerThreadBRep<3>;

    /*
     *  Constructor.
     */
    SimplexDCMesher(PerThreadOutput& m) : m(m) {}

    /* 
     *  Empty cell loader, called by Dual::walk 
     */
    void load(Input* input) {}

    /*  Called by Dual::walk to make branes around simplex edges
     *  whose lowest-dimensional vertex is a face vertex.
     */
    template <Axis::Axis A>
    void load(const std::array<Input*, 2>& ts);

    /*  Called by Dual::walk to make branes around simplex edges
     *  whose lowest-dimensional vertex is an edge vertex. 
     */
    template <Axis::Axis A>
    void load(const std::array<Input*, 4>& ts);


    /*  Called by Dual::walk to make branes around simplex edges
     *  whose lowest-dimensional vertex is a corner vertex. 
     */
    void load(const std::array<Input*, 8>& ts);


    /*
     *  Simplex DC meshing needs to walk the top edges of the tree,
     *  because those include tets that may need to be handled.
     */
    static bool needsTopEdges() { return true; }

protected:
    /*
     *  Adds branes composing a polygon.
     */

    void addPolygon(const SimplexDCCirculator<3, true> circulator);

    PerThreadOutput& m;
};

////////////////////////////////////////////////////////////////////////////////

}   // namespace libfive
