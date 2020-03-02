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

namespace libfive {

// Forward declarations
template <unsigned N, class Leaf> class SimplexTree;
template <unsigned N> struct SimplexDCLeaf;
template <unsigned N> struct SimplexDCMinEdge;
template <unsigned N> struct SimplexDCIntersection;
template <typename... T> class ObjectPool;

/*  It is possible that, while building the tree, a given subspace will have
 *  different copies assigned to the different cells adjacent to it.  This
 *  can make it harder to determine if two subspaces are the same, and 
 *  additionally there seems to be a bug somewhere that can sometimes (at
 *  least when the leaf cells are high-level) result in them actually having
 *  different vertex positions.  Thus, this walker class consolidates all of
 *  the copies of each subspace to equal that from the lowest-indexed cell 
 *  associated with it.*/

template <unsigned N>
class SimplexDCSubConsolidator
{
public:
    struct PerThreadOutput {
        /*  Constructor needed for dual-walking, but we don't actually
         *  use the atomic int, since we're storing information directly 
         *  in the input tree.*/
        PerThreadOutput(std::atomic<uint32_t>& c) {}
    };
    struct Output {
        /*  Another method needed only to fulfill the dual-walking algorithm */
        void collect(std::vector<PerThreadOutput>) {}
    };
    using Input = SimplexTree<N, SimplexDCLeaf<N>>;

    /*
     *  Constructor
     */
    SimplexDCSubConsolidator(PerThreadOutput& m) {}

    /* Empty cell loader, called by Dual::walk */
    void load(Input* input) {}

    /*  Called by Dual::walk to consolidate edge subspaces. */
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 1)> & ts);

    /* Likewise for face subspaces.*/
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 2)>& ts);

    /* And corners. */
    void load(const std::array<Input*, 1 << N>& ts);

    /*
     *  An edge or corner on the top may still have multiple adjacent
     *  cells and need to be consolidated.
     */
    static bool needsTopEdges() { return true; }
};

////////////////////////////////////////////////////////////////////////////////

}   // namespace libfive
