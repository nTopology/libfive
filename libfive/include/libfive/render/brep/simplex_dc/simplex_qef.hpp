/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <Eigen/Eigen>
#include <optional>

#include "libfive/render/brep/util.hpp"
#include "libfive/render/brep/region.hpp"
#include "libfive/render/brep/indexes.hpp"

#ifdef LIBFIVE_VERBOSE_QEF_DEBUG
#include <iostream>
#endif

namespace libfive {

template <unsigned N> struct Intersection;

/*
 *  This functions fulfills a similar purpose as the QEF class, but over an 
 *  n-dimensional simplex rather than a cube, and only solves for
 *  the DC form (an N-dimensional, rather than N+1-dimensional, QEF),
 *  albeit in a constrained manner.  It also does not bother with eigenvalue
 *  cutoffs; they can cause issues when normalizing gradients, and are not
 *  particularly necessary when constraining the solution (though they do still
 *  help keep vertices more evenly spaced, and so may be worth re-adding later).
 *  Some of the QEF functionality for summing various spaces is also not needed,
 *  since once we have simplices we are no longer combining regions.
 */
template <unsigned N>
class SimplexQEF
{
public:
    /*  Proportional distance that the constraints are moved toward the 
     *  barycenter of the simplex, in order to avoid the solutions on disjoint
     *  simplexes being too close to each other.  Larger values produce less
     *  accurate meshes, but have fewer small or thin triangles. 
     */    
    using Point = Eigen::Matrix<double, N, 1>;

    using PointOpt = std::optional<Point>;

    using Matrix = Eigen::Matrix<double, N, N>;
    using RowVector = Eigen::Matrix<double, 1, N>;
    using Vector = Eigen::Matrix<double, N, 1>;

    /*
     * Each constraint is expressed by the requirement that a particular
     * affine expression evaluate to a negative value; this then takes the
     * form A(x,1) < 0, where A is an N + 1-dimensional row vector, and
     * x,1 is the N + 1 dimensional vector consisting of the point in question
     * followed by 1; we can thus express it using only A (though we add helper
     * methods as well).
     */
    struct alignas(16) Constraint {
        int strongestDimension() const {
            int out;
            matrix.array().template head<N>().abs().maxCoeff(&out);
            return out;
        }
        double val(Point pt) const {
            return matrix.template head<N>().dot(pt) + matrix[N];
        }
        Eigen::Matrix<double, N + 1, 1> matrix;
    };

    /*
     *  The input is N+1 N-dimensional vertices.  For best results, the first 
     *  K + 1 columns of the matrix should be identical in N - K rows, for 
     *  1 <= K <= N.  The simplex defined by the vertices should not be 
     *  degenerate (equivalently, the matrix should have full rank); if this is
     *  not fulfilled, results of solve() are undefined.
     */
    SimplexQEF(Eigen::Matrix<double, N, N + 1>&& vertices, 
               double padding, double cutoffRate);

    /*
     *  Accumulate QEFs by summing
     */
     SimplexQEF& operator+=(const Intersection<N>& intersection);

    /*
     *  Solves the QEF, trying to find the point where all of the planes 
     *  intersect at a distance-field value of 0
     */
    Point solve() const;

protected:

    /*
     *  Templated method to emulate a constexpr for loop for unrolling.
     */
    template <unsigned current_mask> 
    PointOpt solveFromMask(
        unsigned currentMasks, 
        unsigned& nextMasksAnyFailed, 
        unsigned& nextMasksAllFailed,
        unsigned cutoffLimit) const;

    /*  Gets which constraints are failed by a given point.*/
    unsigned constraintFailures(const Point& pt) const;

    /*  Checks if each vertex passes the constraints of all the others,
     *  asserting false if not.  Is a no-op when asserts are ignored.*/
    void checkVertConstraints();

    /*
     *  Returns a new QEF corresponding to the subspace defined by the 
     *  constraints set in a bitfield mask field.
     *
     *  For example, this lets you go from a 2D QEF to a 1D QEF by dropping
     *  one constraint, and thus one axis from the matrices.
     */
    template <unsigned mask>
    const SimplexQEF<N - bitcount(mask)>& sub() const;

    /*
     *  Takes a result returned by a subspace QEF, and puts it in this QEF's 
     *  space. 
     */
    template <unsigned mask>
    Point convertFromSub(Eigen::Matrix<double, N - bitcount(mask), 1> pt) const;

    /*
     *  Core QEF solver, which knows nothing of our fancy constraint
     *  systems or distance-field details.  It just takes the matrix and vector
     *  (and a point to minimize towards) and does the math.  It doesn't return
     *  an error, since we don't need it; we're not reducing anything further,
     *  and if we go from higher-dimensional spaces to lower ones, ignoring
     *  those that are not contained in a visited space of the next-highest
     *  dimension or that restrict a visited space by a constraint that its
     *  best point did not fail, the first one we hit with a minimum that
     *  fits all remaining constraints will be our solution, because QEFs 
     *  on convex polytopes are nice like that.
     */
    Point solveRaw(unsigned cutoffLimit) const;

    Point massPoint() const 
    { return mass_point.template head<N>() / mass_point[N]; }

    // Here are our AtA etc.
    Matrix AtA;
    Vector AtB;
    float BtB;

    // And here are our (padded) constraints.
    std::array<Constraint, N + 1> constraints;

    // Mass point, letting us use distance to the mass point as a tiebreaker
    // when one is needed.  (Will be more significant if eigenvalue pruning is
    // restored.)
    Eigen::Matrix<double, N + 1, 1> mass_point;

    // The vertices of the simplex, pulled inward by padding.  While in 
    // principle these are implied by the constraints, storing them explicitly
    // can help avoid numeric error in some cases.
    Eigen::Matrix<double, N, N + 1> paddedVerts;

    int8_t highest_rank = -1;

    double cutoff_rate;

    // Cached subspace QEFs.  Somewhat wasteful of space due to nesting 
    // ignoring commutivity, but we're not keeping them around long-term 
    // and N is a maximum of 3, so it shouldn't be an issue.  We do need
    // special treatment for N == 0 to avoid endless recursion, though.

    // This gets tricky, since we need the special treatment inside the
    // class to avoid a compiler error from recursing below 0, but we
    // can't specialize within the class.  So instead we use constexpr if,
    // auto return type, and decltype.

    struct SubArrayAtPositive {
        std::array<std::optional<SimplexQEF<N - 1>>, N + 1> subs;
        constexpr auto& operator[](size_t n) { return subs[n]; }
    };

    struct SubArrayAtZero{
      // We should not be indexing into it in this case.
    };

    static auto getSubArrayType() {
        assert(false); // This should never be called.
        if constexpr (N == 0) {
            return SubArrayAtZero{};
        }
        else {
            return SubArrayAtPositive{};
        }
    }

    using SubArray = decltype(getSubArrayType());

    mutable SubArray subs;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    friend class SimplexQEF<0>;
    friend class SimplexQEF<1>;
    friend class SimplexQEF<2>;
    friend class SimplexQEF<3>;
};

}   // namespace libfive
