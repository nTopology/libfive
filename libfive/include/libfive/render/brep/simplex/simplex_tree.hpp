/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <array>
#include <atomic>
#include <iostream>
#include <stack>

#include <cstdint>

#include <Eigen/Eigen>
#include <Eigen/StdVector>

#include "libfive/eval/evaluator.hpp"

#include "libfive/render/brep/util.hpp"
#include "libfive/render/brep/xtree.hpp"
#include "libfive/render/brep/object_pool.hpp"
#include "libfive/render/brep/default_new_delete.hpp"
#include "libfive/render/brep/simplex/qef.hpp"
#include "libfive/render/brep/simplex/surface_edge_map.hpp"

namespace libfive {

/* Forward declarations */
template <unsigned N, class Leaf> class SimplexNeighbors;
template<unsigned N, class Leaf> class SimplexTree;
template <unsigned N> class Region;
struct BRepSettings;

template <unsigned N>
struct SimplexLeafSubspace {
    SimplexLeafSubspace();
    void reset();

    /*  Subspace vertex position */
    Eigen::Matrix<double, 1, N> vert;

    /*  Subspace vertex state */
    bool inside;

    /*   Global indices for subspace vertices  */
    std::atomic<uint64_t> index;

    /*  Per-subspace QEF */
    QEF<N> qef;

    /*  SimplexLeafSubspace objects are allocated and released to an
     *  object pool, but can be stored by more than one SubspaceLeaf
     *  at a time (since they represent shared spaces).  We use a
     *  homebrew reference counting system to avoid releasing them to
     *  the pool while they're still in use.  */
    std::atomic<uint32_t> refcount;

    /*  In Simplex DC meshing, it is possible for one subspace vertex to be
     *  collapsed into another one.  This tracks that collapsing; it uses
     *  the same ternary method as subs, but includes only those axes that
     *  are floating in this sub; in the case of a face in 3 dimensions, the
     *  ordering of the axes for this purpose is Q(A), R(A), where A is the 
     *  axis of the face.*/
    unsigned collapseRef;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <unsigned N>
struct SimplexLeaf
{
    SimplexLeaf();
    void reset();

    using Pool = ObjectPool<SimplexLeaf, SimplexLeafSubspace<N>>;

    using ParentPool = ObjectPool<SimplexTree<N, SimplexLeaf>, 
                                  SimplexLeaf, 
                                  SimplexLeafSubspace<N>>;

    void releaseTo(Pool& object_pool);

    /*  One QEF structure per subspace in the leaf, shared between neighbors.
     *  These pointers are owned by an object pool, for fast allocation
     *  and re-use. */
    std::array<std::atomic<SimplexLeafSubspace<N>*>, ipow(3, N)> sub;

    /*  Tape used for evaluation within this leaf */
    std::shared_ptr<Tape> tape;

    /*  Indices of surface vertices, populated when meshing.
     *
     *  The index is a pair of subspace vertex indices.
     *
     *  We can't simply store a fixed number of edges because of
     *  how neighboring cells of varying sizes are meshed.  Instead,
     *  we use a pair of small_vectors to act as a stack-allocated
     *  ordered map. */
    SurfaceEdgeMap<32> surface;

    DEFAULT_OPERATORS_NEW_AND_DELETE

    /*  Represents how far from minimum-size leafs we are */
    unsigned level;

    constexpr static double qefShrink = 1. - 1e-9;
};

template <unsigned N, class Leaf = SimplexLeaf<N>>
class SimplexTree : public XTree<N, SimplexTree<N, Leaf>, Leaf>
{
public:

  using Pool = typename Leaf::ParentPool;

    /*
     *  Simple constructor
     *
     *  Pointers are initialized to nullptr, but other members
     *  are invalid until reset() is called.
     */
    explicit SimplexTree();

    /*
     *  Complete constructor
     */
    explicit SimplexTree(SimplexTree<N, Leaf>* parent, unsigned index,
                         const Region<N>&);

    /*
     *  Constructs an empty SimplexTree
     *
     *  This is called during dual walking, to represent the trees outside
     *  the original bounding box.  The returned tree has an invalid parent 
     *  pointer and Interval::UNKNOWN as its type, and therefore must be
     *  handled accordingly by the walker.
     */
    static std::unique_ptr<SimplexTree> empty();

    /*
     *  Populates type, setting corners, manifold, and done if this region is
     *  fully inside or outside the mode.
     *
     *  Returns a shorter version of the tape that ignores unambiguous clauses.
     */
    std::shared_ptr<Tape> evalInterval(Evaluator* eval,
                                       const std::shared_ptr<Tape>& tape,
                                       Pool& object_pool,
                                       const BRepSettings& settings);

    /*
     *  Evaluates and stores a result at every corner of the cell.
     *  Sets type to FILLED / EMPTY / AMBIGUOUS based on the corner values.
     *  Then, solves for vertex position, populating AtA / AtB / BtB.
     */
    void evalLeaf(Evaluator* eval,
                  const std::shared_ptr<Tape>& tape,
                  Pool& object_pool,
                  const SimplexNeighbors<N, Leaf>& neighbors,
                  const BRepSettings& settings);

    /*
     *  If all children are present, then collapse based on the error
     *  metrics from the combined QEF (or interval filled / empty state).
     *
     *  Returns false if any children are yet to come, true otherwise.
     */
    bool collectChildren(Evaluator* eval,
                         const std::shared_ptr<Tape>& tape,
                         Pool& object_pool,
                         const BRepSettings& settings);

    /*  Looks up the cell's level for purposes of vertex placement,
     *  returning 0 or more for LEAF / EMPTY / FILLED cells (depending
     *  on how many other leafs were merged into them; 0 is the smallest
     *  leaf).
     *
     *  Returns UINT32_MAX for UNKNOWN cells, which should only be created
     *  with SimplexTree::empty() and are used around the borders of the
     *  model to include those edges.
     *
     *  Triggers an assertion failure if called on a BRANCH cell.
     */
    uint32_t leafLevel() const;

    /*
     *  Assigns leaf->sub[*]->index to a array of unique integers for every leaf
     *  in the tree, starting at 1.  This provides a globally unique
     *  identifier for every subspace vertex, which is used when making edges.
     *
     *  Settings are used for cancellation and worker count.
     */
    void assignIndices(const BRepSettings& settings) const;

    /*
     *  Releases this tree and any leaf objects to the given object pool
     */
    void releaseTo(Pool& object_pool);

    /*  Boilerplate for an object that contains an Eigen struct  */
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /*  Helper typedef for N-dimensional column vector */
    typedef Eigen::Matrix<double, N, 1> Vec;

    static bool hasSingletons() { return false; }
    static SimplexTree<N, Leaf>* singletonEmpty() { return nullptr; }
    static SimplexTree<N, Leaf>* singletonFilled() { return nullptr; }
    static bool isSingleton(const SimplexTree<N, Leaf>*) { return false; }

    static constexpr unsigned leafSize = ipow(3, N);

    /*
     *  Calculate and store whether each vertex is inside or outside
     *  This populates leaf->sub[i]->inside, for i in 0..ipow(3, N)
     */
    void saveVertexSigns(Evaluator* eval,
        const Tape::Handle& tape,
        const std::array<bool, leafSize>& already_solved);
protected:

    /*  We make a copy of the children when collecting them, in order to avoid
     *  atomics.  This typedef makes it easier to pass that on to other methods.
     */
    using RawChildArray = std::array<SimplexTree<N, Leaf>*, 1 << N>;

    /*
     *  Sets this->type to EMPTY / FILLED / AMBIGUOUS depending on
     *  the vertex signs, which must be populated.
     */
    void checkVertexSigns();

    /*
     *  Populates this->leaf->sub[i]->qef for every corner subspace,
     *  then solves for vertex position and signs.
     *
     *  Only corners are evaluated + populated; faces / edges / volumes
     *  are initialized to zero, because they'll be constructed by accumulation
     *  as we walk up the tree.
     */
    void findLeafVertices(Evaluator* eval,
                          const Tape::Handle& tape,
                          Pool& object_pool,
                          const SimplexNeighbors<N, Leaf>& neighbors,
                          const BRepSettings& settings);

    /*  
     *  This helper struct allows us to easily reference a subspace of a child.
     *  It is actually independent of the SimplexTree template arguments, so at
     *  some point we may wish to move it to its own file.
     */
    template <unsigned K>
    struct ChildSubArray
    {
        ChildSubArray(bool b) : components{ b, b, b, b, b } {}
        ChildSubArray() = default;
        // Helper function to perform template differentiation.
        static auto getComponentTypeRep() {
            if constexpr (K == 1) {
                return false; // Child is bool
            }
            else {
                // Child is next dimension down.
                return ChildSubArray<K - 1>(false); 
            }
        }

        using ComponentType = decltype(ChildSubArray<K>::getComponentTypeRep());

        // We have 5 elements in each direction: Fixed lower of the lower
        // child, floating of the lower child, fixed lower of the upper child/
        // upper of the lower child, floating of the upper child, and fixed
        // upper of the upper child.
        std::array<ComponentType, 5> components;

        bool& operator[](std::array<int, K> ar) {
            if constexpr (K == 1) {
                return components[ar[0]];
            }
            else {
                std::array<int, K - 1> newArray;
                std::copy(ar.begin() + 1, ar.end(), newArray.begin());
                return components[ar[0]][newArray];
            }
        }
        bool operator[](std::array<int, K> ar) const {
            if constexpr (K == 1) {
                return components[ar[0]];
            }
            else {
                std::array<int, K - 1> newArray;
                std::copy(ar.begin() + 1, ar.end(), newArray.begin());
                return components[ar[0]][newArray];
            }
        }
        void flip() {
            for (auto& component : components) {
                if constexpr (K == 1) {
                    component = !component;
                }
                else {
                    component.flip();
                }
            }
        }
    };

    /*
     *  Unwraps the atomic pointers, returning an array of plain pointers
     *  Only use this if multiple threads won't be messing with the pointers
     *  simultaneously!
     */
    std::array<SimplexLeafSubspace<N>*, ipow(3, N)> getLeafSubs() const;

    /*
     *  Indicates whether this node can be collapsed without changing the
     *  topology of the piece.  "children" should refer to the children of
     *  *this, and must all be leaves.  *this must also have its leaf set,
     *  and its subspaces populated.
     */
     bool isCollapsible(const RawChildArray& children);

    /*
     *  Indicates whether the designated subspace can be collapsed without
     *  changing topology, provided that its sub-subspaces can be collapsed.
     *  (Otherwise, it is still safe to call.)  Conditions on *this are the 
     *  same as for "isCollapsible".
     */
     bool subIsCollapsible(unsigned subIndex, 
                           const ChildSubArray<N>& childSubsInside);
};

}   // namespace libfive
