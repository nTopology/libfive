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
#include <assert.h>

#include "libfive/render/axes.hpp"

namespace libfive {

// Forward declarations
template <unsigned N, class Leaf> class SimplexTree;
template <unsigned N> struct SimplexDCLeaf;
template <unsigned N> struct SimplexDCEdge;
template <unsigned N> struct SimplexDCIntersection;
template <typename... T> class ObjectPool;

/*  This dual-walker class simply allocates and assigns instances of the 
 *  SimplexDCEdge class to appropriate cells in the subtree.  We do this
 *  as a separate stage rather than as part of leaf construction in order
 *  to ensure that all cells bordering an edge will have the same 
 *  SimplexDCEdge struct (allowing us to later use this fact to ensure
 *  we are making each simplex vertex exactly once and it will be accessible
 *  to all its edges).  (If we used the neighbors system, it would be likely
 *  but not certain that any given edge would have the same object for all its
 *  neighbors.)  Cells whose edge contains smaller edges use a lock-free stack
 *  of non-owning edge pointers instead, and this likewise adds the edge to
 *  such stacks.*/

template <unsigned N>
class SimplexDCEdges
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
    };
    using Input = SimplexTree<N, SimplexDCLeaf<N>>;

    using Pool = ObjectPool<SimplexDCEdge<N>, SimplexDCIntersection<N>>;


    /*
     *  Constructor
     */
    SimplexDCEdges(PerThreadOutput& m, Pool& pool);

    SimplexDCEdges(SimplexDCEdges&&) = default;

    /*
     *  Destructor used to ensure the parent pool claims the pool of this. 
     */
    ~SimplexDCEdges();

    /* Empty cell loader, called by Dual::walk */
    void load(Input* input) {}

    /*  Called by Dual::walk to assign the appropriate edges. */
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 1)> & ts);

    /* Empty face loader, called by Dual::walk, used in 3d only.*/
    template <Axis::Axis A>
    void load(const std::array<Input*, 1 << (N - 2)> & ts) { assert(N == 3); }

    /* Empty corner loader */
    void load(const std::array<Input*, 1 << N> & ts) {}

    /*
     *  Simplex DC meshing needs to walk the top edges of the tree,
     *  because those include tets that may need to be handled.
     */
    static bool needsTopEdges() { return true; }

protected:
    Pool& parent_pool;
    Pool my_pool;
};

////////////////////////////////////////////////////////////////////////////////

}   // namespace libfive
