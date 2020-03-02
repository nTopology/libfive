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
template <unsigned N> struct SimplexDCLeaf;
template <unsigned N> struct SimplexLeafSubspace;
class Evaluator;

/*  This dual-walker class determines when subspace vertices (calculated when
 *  originally building the tree) need to be "collapsed" to reference 
 *  lower-dimensional subspaces, and does so.  Because the collapsing of
 *  higher-dimensional spaces can depend on which lower-dimensional spaces
 *  have already been collapsed (and how), it is done in two stages; first
 *  with dimensions 0-2 (which can be done together, as the rules for collapsing
 *  spaces of dimensions 0-1 are simple enough to do on the fly), and then 
 *  (for N == 3) for dimension 3.  Thus, there is a second template parameter K,
 *  which should fulfill 2 <= K <= N.*/

template <unsigned N, unsigned K>
class SimplexDCCollapser
{
public:
    struct PerThreadOutput {
        /*  Constructor needed for dual-walking, but we don't actually
         *  use the atomic int, since we're storing edge-specific
         *  information directly in the input tree.*/
        PerThreadOutput(std::atomic<uint32_t>& c) {}
    };
    struct Output {
        /*  Another method needed only to fulfill the dual-walking algorithm */
        void collect(std::vector<PerThreadOutput>) {}
    };    using Input = SimplexTree<N, SimplexDCLeaf<N>>;

    /*
     *  Constructs an intersecter that owns an evaluator,
     *  which is built from the given tree.
     */
    SimplexDCCollapser(PerThreadOutput& m, Tree t, double boundCutoff);

    /*
     *  Constructs an intersecter that has borrowed an evaluator,
     *  which is useful in cases where constructing evaluators
     *  is expensive and they should be re-used.
     */
    SimplexDCCollapser(PerThreadOutput& m, Evaluator* es, double boundCutoff);

    ~SimplexDCCollapser();

    /*  Cell loader, called by Dual::walk.  Collapses cell vertices.
     *  Trivial for K < N, and redirects to the face loader for
     *  N == 2.*/
    void load(Input* t);

    /*  Face loader, collapses face vertices.  Called by Dual::Walk for
     *  N == 3, and by the cell loader (with A == Axis::Z) for N == 2.*/
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 2)>& ts);

    /*  Edge loader, collapses edge vertices. */
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 1)>& ts);

    /*  Empty corner loader. */
    void load(const std::array<Input*, 1 << N>& ts) {}

    /*
     *  Simplex DC meshing needs to walk the top edges of the tree,
     *  because those include subspaces that may need to be handled.  However,
     *  if K == 3, we are only handling cells, and there are no cells in the
     *  top of the tree
     */
    static bool needsTopEdges() { return K == 2; }

protected:
    /*  Compares two subspaces that a higher-dimensional subspace may collapse
     *  to, in order to determine which is the better candidate.  This will be
     *  called only in the case where a high-dimension subspace has a vertex
     *  that should collapse to a lower-dimensional subspace, but there are
     *  multiple intermediate subspaces that do not collapse at all (so it 
     *  should be fairly rare).  Returns true if candidate A is better,
     *  and false if candidate B is better; these should be the neighbor
     *  indices (not relative to toCollapse) of the candidate subspaces to
     *  collapse to.  Result is undefined if the two candidates do not
     *  share at least a corner, or if they are the same.*/
    bool compareCollapseCandidates(
        const Input* cell, const SimplexLeafSubspace<N>& toCollapse,
        unsigned candidateA, unsigned candidateB);

    Evaluator* eval;
    bool owned;
    double bound_cutoff;
};

////////////////////////////////////////////////////////////////////////////////

// Some methods are trivial, and can be included here to allow them to be
// inlined and skip a function call.

template <>
inline void SimplexDCCollapser<2, 2>::load(Input* input) 
{ load<Axis::Z>(std::array{input}); }

template <>
inline void SimplexDCCollapser<3, 2>::load(Input*) {}

template <>
template <Axis::Axis A>
void SimplexDCCollapser<3, 3>::load(const std::array<Input*, 2>&) {}

template <>
template <Axis::Axis A>
void SimplexDCCollapser<3, 3>::load(const std::array<Input*, 4>&) {}

}   // namespace libfive
