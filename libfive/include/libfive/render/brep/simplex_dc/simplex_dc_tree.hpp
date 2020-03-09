/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once
#include "libfive/render/brep/edge_tables.hpp"
#include "libfive/render/brep/simplex/simplex_tree.hpp"
#include "libfive/render/brep/dc/intersection.hpp"
#include "libfive/render/axes.hpp"

namespace libfive {

/* Forward declaration */
template <unsigned N> struct Intersection;

constexpr int _factorial(unsigned N)
{
  if (N == 0) { return 1; }
  else { return N * _factorial(N - 1); }
}

/*  Returns the number of simplices adjacent to an edge in n dimensions */
constexpr int _edgeSimplices(unsigned N)
{
  return ipow(2, N) * _factorial(N - 1);
}

/*  Returns the number of edges in an N-dimensional simplex */
constexpr int _simplexEdges(unsigned N)
{
  return (N * (N + 1) / 2);
}


template <unsigned N>
struct OrientationChecker; // A helper class to allow the SimplexDCIntersection
                           // to check if two points are oriented as expected.

template <>
struct OrientationChecker<2>
{
    Eigen::Matrix2d endpoints;
    void reset() { endpoints.setZero(); }
    void setEndpts(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
        { endpoints << a, b; }

    // Checks if the line between a and b goes between the two endpoints 
    // (assuming they are on opposite sides of the infinite line generated 
    // by those endpoints; otherwise, result is undefined).   Returns 0 if it 
    // does, 1 if it is past endpt1, and 2 if it is past the other endpoint.  
    // The order in which the endpoints were set is not relevent here, and 
    // neither is the order of a and b.  If the line between a and b hits an 
    // endpoint, it is considered to go past it.
    int check(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
        const Eigen::Vector2d& endpt1) const;
};

template <>
struct OrientationChecker<3>
{
    Eigen::Vector3d direction; // Outward direction.
    void reset() { direction.setZero(); }
    void setEndpts(
        const Eigen::Vector3d& outside, const Eigen::Vector3d& inside)
    { direction = outside - inside; }

    // Checks whether the points passed are indeed clockwise and 
    // counterclockwise when looking from the outside point toward the
    // inside one.  Returns true if they are, false if they aren't (false
    // if all four points are coplanar).  Assumes that "center" is on the
    // line between the inside and outside points, in order to avoid having
    // to explicitly store the two points (as center can be taken from the
    // SimplexDCIntersection that holds this).
    bool check(const Eigen::Vector3d& clockwisePt, 
               const Eigen::Vector3d& counterclockwisePt,
               const Eigen::Vector3d& center) const;
};

template <unsigned N>
struct SimplexDCIntersection : Intersection<N>
{
    /*  SimplexDCIntersection objects are allocated and released to an
     *  object pool, but can be stored by more than one SubspaceLeaf
     *  at a time (since they represent intersections on shared edges).  
     *  We use  a homebrew reference counting system to avoid releasing
     *  them to the pool while they're still in use.  */
    SimplexDCIntersection();
    void reset();

    std::atomic<uint32_t> refcount = 0;

    /*  Global indices for edge intersection vertices */
    std::atomic<uint64_t> index = 0;

    OrientationChecker<N> orientationChecker;
};

template <unsigned N>
struct DCSimplex
{
    std::array<std::atomic<SimplexDCIntersection<N>*>, 
               _simplexEdges(N)> intersections;

    /* The vertices of the simplex are indexed by the dimensionalities of the 
     * subspaces they came from; this gives the intersection for a given 
     * pair of vertices.*/
    const SimplexDCIntersection<N>* intersection(unsigned a, unsigned b) const;

    /* Assigns the given intersection, if it is currently empty.  Returns the 
     * intersection currently held for the appropriate edge, and a boolean
     * indicating whether the insertion was successful.  This handles both the
     * ref counting for the inserted intersection and thread-safety.*/
    std::pair<const SimplexDCIntersection<N>*, bool> insertIntersection(
        unsigned a, unsigned b, SimplexDCIntersection<N>* intersection);

    /* Number of non-null intersections, for assertions/debugging.*/
    int intersectionCount() const {
        return std::count_if(intersections.begin(), intersections.end(),
            [](const auto& atomic_ptr) {return atomic_ptr.load() != nullptr; });
    }

    /*  Simplex vertex position */
    Eigen::Matrix<double, N, 1> vert;

    /*  Global indices for simplex vertices  */
    std::atomic<uint64_t> index = 0;

    DEFAULT_OPERATORS_NEW_AND_DELETE
};

template <unsigned N>
struct SimplexDCMinEdge
{
    /* The index for a simplex is assigned with the lowest N bits indicating
     * the cell center and corner used; in 2d this is the lowest bit for the
     * corner and the next-lowest for the cell, while in 3d the second- and
     * third-lowest are for the Q and R values (respectively, with respect to
     * the edge's axis) of the cell.  The rest of the index is a number from
     * 0 to _factorial(N - 1); in 2d this is trivial, while in 3d it refers
     * to whether the face's axis is Q of the edge's axis (1) or R or the
     * edge's axis (0).  If this would result in a simplex using a face that
     * is internal to a merged cell, the value of the resulting simplex is
     * undefined.*/

    std::array<DCSimplex<N>, _edgeSimplices(N)> simplices;

    std::atomic<uint32_t> refcount;

    /* The lowest value along the edge of the coordinate corresponding to the
     * edge's axis.  Used to sort. */
    double lowPt;

    SimplexDCMinEdge();
    void reset();

    using Pool = ObjectPool<SimplexDCMinEdge<N>, SimplexDCIntersection<N>>;

    void releaseTo(Pool& object_pool);

    /*  Gets the simplex corresponding to a particular corner, particular cell,
     *  and particular face axis.  Behavior is undefined if edgeAxis does not
     *  correspond to the actual axis of this edge (which will be the axis 
     *  passed to SimplexDCLeaf::edge to get this edge), if faceAxis is the 
     *  same as edgeAxis, if either is not a valid axis, or if cell is not 
     *  a valid N-dimensional bitfield.  The bit of cell corresponding to
     *  edgeAxis is ignored; in 2 dimensions, faceAxis is ignored as well.*/
    constexpr DCSimplex<N>& simplex(Axis::Axis edgeAxis, Axis::Axis faceAxis, 
                                    unsigned cell, bool isUpperCorner);

    /*  Gets the simplex assuming the cell index has been "reduced" to N-1
     *  dimensions by ignoring the axis corresponding to the edge.  The
     *  remaining bits are rotated so that the lowest bit is the one after
     *  the edge axis (rotating if necessary).*/
    constexpr DCSimplex<N>& simplexWithEdgeReducedCell(
        Axis::Axis edgeAxis, Axis::Axis faceAxis,
        unsigned cell, bool isUpperCorner);

    /*  This variant assumes that the cell and corner are relative to the
     *  face rather than the edge; the cell is thus a boolean, and the corner
     *  has two bits.  This method has undefined behavior for N == 2.
     */
    constexpr DCSimplex<N>& simplexWithFaceReducedCell(
        Axis::Axis edgeAxis, Axis::Axis faceAxis,
        bool isUpperCell, unsigned reducedCorner);

    DEFAULT_OPERATORS_NEW_AND_DELETE

    /*  Non-owning vector of edges.  It is thread-safe with respect to multiple
     *  pushes/sorts, but not with respect to reading and pushing at the same 
     *  time.  Reading and sorting at the same time is thread-safe only if
     *  the sort() function has already returned (after being called on this
     *  object) at least once.  It should not be pushed to after being sorted;
     *  if it is, the validity of the sort is no longer guaranteed, even if 
     *  sort() is called again.*/
    class EdgeVec {
    public:
        void push_back(SimplexDCMinEdge* ptr);
        auto begin() const { return vec.begin(); }
        auto end() const { return vec.end(); }
        auto rbegin() const { return vec.rbegin(); }
        auto rend() const { return vec.rend(); }
        auto empty() const { return vec.empty(); }
        // Sorts according to the lowPt member of its members.
        void sort();
    private:
        bool isSorted;
        std::mutex mut;
        std::vector<SimplexDCMinEdge*> vec;
    };
};

/*  This is intentionally set to use the stack version by default.  That way,
 *  the potentially-multithreaded use of the stack will not need to use any
 *  thread-unsafe variant operations, while switching to the single-pointer
 *  option will only happen when the corresponding edge is minimal and 
 *  therefore will only be accessed once.  Pointers in the EdgeVec are
 *  non-owning, while those held directly are owning.*/
template <unsigned N>
using SimplexDCEdge = std::variant<typename SimplexDCMinEdge<N>::EdgeVec,
                                   SimplexDCMinEdge<N>*>;

template <unsigned N>
struct SimplexDCLeaf
{
    SimplexDCLeaf();
    void reset();

    using Pool = ObjectPool<SimplexDCLeaf, 
                            SimplexLeafSubspace<N>,
                            SimplexDCMinEdge<N>,
                            SimplexDCIntersection<N>>;

    using ParentPool = ObjectPool<SimplexTree<N, SimplexDCLeaf>, 
                                  SimplexDCLeaf, 
                                  SimplexLeafSubspace<N>,
                                  SimplexDCMinEdge<N>,
                                  SimplexDCIntersection<N>>;

    void releaseTo(Pool& object_pool);

    /*  One QEF structure per subspace in the leaf, shared between neighbors.
     *  These pointers are owned by an object pool, for fast allocation
     *  and re-use. */
    std::array<std::atomic<SimplexLeafSubspace<N>*>, ipow(3, N)> sub;

    /*  Tape used for evaluation within this leaf */
    std::shared_ptr<Tape> tape;

    /*  In addition to general subspace information, we also need edges to hold
     *  the actual simplices; by holding them by edge rather than cell, we
     *  ensure that there is a fixed number per container regardless of varying
     *  cell size.  The index is expressed as toIndex(A) * 2^(N-1) + B, where A 
     *  is the index of the axis of the edge, and B is N bits corresponding in 
     *  2d to whether the perpendicular value is at the lower (0) or upper (1) 
     *  end of the cell, and in 3d is 2R+Q, where Q and R represent the value 
     *  of Q(A) and R(A).*/
    std::array<SimplexDCEdge<N>, _edges(N)> edges;

    DEFAULT_OPERATORS_NEW_AND_DELETE

    /*  Represents how far from minimum-size leafs we are */
    unsigned level;

    /*  Gets the edge along a given axis that connects to the specified corner.
     *  edge(axis) and edge(axis^corner) thus give the same result.  If axis
     *  is not a valid axis, or corner is not a valid n-dimensional corner,
     *  behavior is undefined.*/
    constexpr SimplexDCEdge<N>& edge(Axis::Axis axis, unsigned corner);

    /*  This variant takes a reduced N-1-dimensional corner, referring to the
     *  non-axis dimensions beginning with the one after the axis dimension
     *  (and wrapping around if needed).*/
    constexpr SimplexDCEdge<N>& edgeFromReduced(Axis::Axis axis, 
                                                unsigned reducedCorner) 
    {
        return edges[(Axis::toIndex(axis) << (N - 1)) | reducedCorner];
    }

    /*  This variant gets the edge based on adjacency to a face (determined
     *  by an axis and whether it takes the higher or lower value on that axis),
     *  an axis for the edge, and whether the edge takes the higher or lower
     *  value on the unused axis.  This has undefined behavior when N == 3, and
     *  if edgeAxis == faceAxis or either is not a valid axis.*/
    constexpr SimplexDCEdge<N>& edgeFromFaceAndIndex(
        Axis::Axis edgeAxis, Axis::Axis faceAxis,
        bool upperEdge, bool upperFace);
};

template <unsigned N>
using SimplexDCTree = SimplexTree<N, SimplexDCLeaf<N>>;

/////////////// Constexpr method definitions ///////////////////

namespace {
constexpr unsigned lowestBits(unsigned N) {
    auto out = 0;
    for (auto i = 0u; i < N; ++i) {
        out |= (1 << i);
    }
    return out;
}
}

template<unsigned N>
constexpr SimplexDCEdge<N>& SimplexDCLeaf<N>::edge(Axis::Axis axis, unsigned corner)
{
    corner |= (corner << N);
    corner >>= (Axis::toIndex(axis) + 1);
    corner &= lowestBits(N - 1);
    return edgeFromReduced(axis, corner);
}

template<unsigned N>
inline constexpr SimplexDCEdge<N>& SimplexDCLeaf<N>::edgeFromFaceAndIndex(
    Axis::Axis edgeAxis, Axis::Axis faceAxis, bool upperEdge, bool upperFace)
{
    auto corner = upperEdge ? lowestBits(N) : 0;
    if (upperFace) {
        corner |= faceAxis;
    }
    else {
        corner &= ~faceAxis;
    }
    return edge(edgeAxis, corner);
}

template<unsigned N>
inline constexpr DCSimplex<N>& SimplexDCMinEdge<N>::simplex(
    Axis::Axis edgeAxis, Axis::Axis faceAxis, unsigned cell, bool isUpperCorner)
{
    cell |= (cell << N);
    cell >>= (Axis::toIndex(edgeAxis) + 1);
    cell &= lowestBits(N - 1);
    return simplexWithEdgeReducedCell(edgeAxis, faceAxis, cell, isUpperCorner);
}

template<unsigned N>
inline constexpr DCSimplex<N>& SimplexDCMinEdge<N>::simplexWithEdgeReducedCell(
    Axis::Axis edgeAxis, Axis::Axis faceAxis, unsigned cell, bool isUpperCorner)
{
    auto index = (cell << 1) | int(isUpperCorner);
    if constexpr (N == 3) {
        if (faceAxis == Axis::Q(edgeAxis)) {
            index |= (1 << N);
        }
    }
    return simplices[index];
}

template<unsigned N>
inline constexpr DCSimplex<N>& SimplexDCMinEdge<N>::simplexWithFaceReducedCell(
    Axis::Axis edgeAxis, Axis::Axis faceAxis, 
    bool isUpperCell, unsigned reducedCorner)
{
    auto shiftValue = Axis::toIndex(faceAxis) + 1;
    auto fullCorner = reducedCorner << shiftValue;
    fullCorner |= (reducedCorner >> (N - shiftValue));
    if (!isUpperCell) {
        fullCorner |= faceAxis;
    }
    auto cell = (~fullCorner) & lowestBits(N);
    return simplex(edgeAxis, faceAxis, cell, bool(fullCorner & edgeAxis));
}

}   // namespace libfive
