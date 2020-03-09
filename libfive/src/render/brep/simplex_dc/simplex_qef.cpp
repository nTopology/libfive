/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "libfive/render/brep/simplex_dc/simplex_qef.hpp"
#include "libfive/render/brep/dc/intersection.hpp"

namespace libfive {

template <unsigned N>
SimplexQEF<N>::SimplexQEF(Eigen::Matrix<double, N, N + 1>&& vertices, 
                          double padding) :
    AtA(Matrix::Zero()), AtB(Vector::Zero()), BtB(0.), 
    mass_point(Eigen::Matrix<double, N + 1, 1>::Zero())
{
    auto center = vertices.rowwise().mean();
    auto diff = (vertices.colwise() - center);
    paddedVerts = vertices -= padding * diff;

    // Our formula to produce a constraint only gives us where it is 0; to
    // determine which side is negative, we need to test one of them.  We
    // choose to test constraint N + 1 (the one associated with vertices
    // from 0 to N), since that will ideally correspond to an axis-aligned
    // side of the simplex and therefore minimize the effect of numeric
    // error.
    auto negateAtN = false;
    for (int i = N; i >= 0; --i) {
        auto shiftedVert = [&](int idx)
        {return paddedVerts.col((i + idx) % (N + 1)); };
        static_assert(N <= 3);
        Vector perp;
        double val;
        if constexpr (N >= 1) {
            if constexpr (N >= 2) {
                auto diff1 = shiftedVert(2) - shiftedVert(1);
                if constexpr (N == 3) {
                    auto diff2 = shiftedVert(3) - shiftedVert(1);
                    perp = diff1.cross(diff2);
                }
                else {
                    perp << diff1.y(), -diff1.x();
                }
            }
            else {
                perp << 1.;
            }
            val = perp.dot(shiftedVert(1));
            // Check that our simplex is not degenerate.
            //            assert(val != perp.dot(paddedVerts.col(i)));
            if (i == N && val < perp.dot(paddedVerts.col(N))) {
                // perp causes paddedVerts[N] to be higher than our constant
                // offset, so after the offset it'll still be positive, so
                // we need to negate.
                negateAtN = true;
            }
        }
        else {
            // In 0d, we have one constraint corresponding to the side of the
            // simplex not adjacent to our point, but we never use it for
            // anything, and it has no real meaning since all constraints are
            // only up to a constant positive factor.  So just set it to -1.
            val = -1.;
        }

        // If N is even, rotating the vertices will not change the sign
        // of the determinant that we're effectively comparing to 0 here, 
        // so we negate either everywhere or nowhere.  If N is odd, rotating
        // the vertices by an even amount will not change the sign of the 
        // determinant, but by an odd amount will flip it, so we alternate
        // which ones we flip.
        auto sameAsN = (N % 2 == 0 || i % 2 == 1);
        if (sameAsN == negateAtN) {
            constraints[i].matrix << -perp, val;
        }
        else {
            constraints[i].matrix << perp, -val;
        }
    }
}

template <unsigned N>
SimplexQEF<N>& SimplexQEF<N>::operator+=(const Intersection<N>& intersection) {
    AtA += intersection.AtA;
    AtB += intersection.AtB;
    BtB += intersection.BtB;
    if (intersection.rank > highest_rank) {
        highest_rank = intersection.rank;
        mass_point = intersection.mass_point;
    }
    else if (intersection.rank == highest_rank) {
        mass_point += intersection.mass_point;
    }

    return *this;
}

template <unsigned N>
template <unsigned mask>
const SimplexQEF<N - bitcount(mask)>& SimplexQEF<N>::sub() const
{
    constexpr auto highBit = highestbit(mask);
    static_assert(highBit < N + 1, "Too many constraints");
    static_assert(bitcount(mask) <= N, "Using all constraints, not compatible");

    if constexpr (mask == 0) {
        return *this;
    }
    else if constexpr (bitcount(mask) == 1) {
        auto& out = subs[highBit];
        if (!out.has_value()) {
            const auto& constraint = constraints[highBit];
            auto axis = constraint.strongestDimension();
            assert(axis < N);

            // Get a "reduced" constraint matrix that we will use to adjust
            // the values of outVal's members
            const auto& matrix = constraint.matrix;
            auto normalized = matrix / matrix[axis];
            Eigen::Matrix<double, N, 1> reduced;
            reduced << normalized.head(axis), normalized.tail(N - axis);
            // Generate our restricted set of simplex vertices; we regenerate
            // our lower-dimensional constraints from those rather than by
            // intersecting higher-dimensional constraints in order to avoid
            // floating-point issues when the simplex is nearly degenerate.
            Eigen::Matrix<double, N - 1, N> newVerts;
            if constexpr (N > 1) {
                newVerts <<
                    paddedVerts.topLeftCorner(axis, highBit),
                    paddedVerts.topRightCorner(axis, N - highBit),
                    paddedVerts.bottomLeftCorner(N - 1 - axis, highBit),
                    paddedVerts.bottomRightCorner(N - 1 - axis, N - highBit);
            }
            out.emplace(std::move(newVerts), 0);
            auto& outVal = out.value();
            if constexpr (N > 1) {
                auto reducedHead = reduced.template head<N - 1>();

                // Start by copying all of AtA except the 
                // row and column at axis.
                outVal.AtA <<
                    AtA.topLeftCorner(axis, axis),
                    AtA.topRightCorner(axis, N - 1 - axis),
                    AtA.bottomLeftCorner(N - 1 - axis, axis),
                    AtA.bottomRightCorner(N - 1 - axis, N - 1 - axis);

                // Then use reduced to compensate for the lost dimension.
                Eigen::Matrix<double, 1, N - 1> reducedFromAtA;
                const auto& axisRow = AtA.row(axis);
                reducedFromAtA << 
                    axisRow.head(axis), axisRow.tail(N - 1 - axis);
                auto adjustment1 = reducedHead * reducedFromAtA;
                outVal.AtA -= adjustment1;
                outVal.AtA -= adjustment1.transpose();
                outVal.AtA += reducedHead * reducedHead.transpose() *
                    AtA(axis, axis);
                // Once we're done, use eval to ensure that outVal.AtA doesn't
                // run into problems when this function returns and the 
                // intermediate matrices go out of scope.
                outVal.AtA = outVal.AtA.eval();

                outVal.AtB << AtB.head(axis), AtB.tail(N - 1 - axis);
                outVal.AtB -= AtB[axis] * reducedHead;
                // Note that this is += rather than the -= we used earlier
                // when using 'reduced' only once.  This is because we're
                // applying elements from AtA to AtB, and they are on opposite
                // sides of the equation (which we are trying to solve)
                // AtAx=AtB; if moved to the same side, one would be negated.
                outVal.AtB += reducedFromAtA.transpose() * reduced[N - 1];
                outVal.AtB -= AtA(axis, axis) * reducedHead * reduced[N - 1];
                outVal.AtB = outVal.AtB.eval();
            }

            outVal.BtB = BtB; 
            outVal.BtB += AtA(axis, axis) * reduced[N - 1] * reduced[N - 1];
            outVal.BtB += 2 * AtB[axis] * reduced[N - 1];

            auto massPoint = mass_point.template head<N>();
            auto adjustment = 
                (massPoint.dot(matrix.template head<N>()) + matrix[N]) /
                matrix.template head<N>().squaredNorm();
            auto adjustedMassPoint = (massPoint -
                matrix.template head<N>() * adjustment);
            outVal.mass_point << 
                adjustedMassPoint.head(axis), 
                adjustedMassPoint.tail(N - 1 - axis), 
                mass_point(N);

            // We don't need to copy the rank, since we won't be adding 
            // intersections directly to the subspace QEF.            
        }

        return out.value();
    }
    else {
        constexpr auto fromHighest = 1 << highBit;
        static_assert(fromHighest & mask);
        return sub<fromHighest>().sub<mask & ~fromHighest>();
    }
}

template<unsigned N>
template<unsigned mask>
typename SimplexQEF<N>::Point SimplexQEF<N>::convertFromSub(
    Eigen::Matrix<double, N - bitcount(mask), 1> pt) const
{
    constexpr auto highBit = highestbit(mask);
    static_assert(highBit < N + 1, "Too many constraints");
    static_assert(bitcount(mask) <= N, "Using all constraints, not compatible");
    if constexpr (mask == 0) {
        return pt;
    }
    else if constexpr (bitcount(mask) == 1) {
        const auto& constraint = constraints[highBit];
        auto axis = constraint.strongestDimension();
        assert(axis < N);
        // Get a "reduced" constraint matrix that we will use to get
        // the value of our point on the missing dimension.
        const auto& matrix = constraint.matrix;
        auto normalized = matrix / matrix[axis];
        typename SimplexQEF<N - 1>::Constraint reduced;
        reduced.matrix << normalized.head(axis), normalized.tail(N - axis);
        auto coord = -reduced.val(pt);
        Point out;
        out << pt.head(axis), coord, pt.tail(N - 1 - axis);
        return out;
    }
    else {
        constexpr auto fromHighest = 1 << highBit;
        static_assert(fromHighest & mask);
        auto fromNext = sub<fromHighest>()
            .convertFromSub<mask & ~fromHighest>(pt);
        return convertFromSub<fromHighest>(fromNext);
    }
}

template<unsigned N>
template<unsigned current_mask>
typename SimplexQEF<N>::PointOpt SimplexQEF<N>::solveFromMask(
    unsigned currentMasks, 
    unsigned& nextMasksAnyFailed, 
    unsigned& nextMasksAllFailed) const
{
    if ((1 << current_mask) & currentMasks) {
        auto rawResult = sub<current_mask>().solveRaw();
        auto subFailures = sub<current_mask>().constraintFailures(rawResult);
        // assert(bitcount(subFailures) + bitcount(current_mask) <= N);
        for (unsigned i = 0, subI = 0; i <= N; ++i) {
            assert(i >= subI);
            assert(i <= subI + bitcount(current_mask));
            auto shifted = 1 << i;
            if (shifted & current_mask) {
                // Don't check the constraints that define this subspace.
                continue;
            }
            auto constraintMask = 1 << (current_mask | shifted);
            if (subFailures & (1 << subI)) {
                nextMasksAnyFailed |= constraintMask;
            }
            else {
                nextMasksAllFailed &= ~constraintMask;
            }
            ++subI;
        }
        if (subFailures == 0) {
            auto result = convertFromSub<current_mask>(rawResult);
            return result;
        }
    }
    constexpr auto next_mask = current_mask + 1;
    if constexpr (bitcount(next_mask) > N) {
        return {};
    }
    else {
        return solveFromMask<next_mask>(
            currentMasks, nextMasksAnyFailed, nextMasksAllFailed);
    }
}

template <>
inline typename SimplexQEF<0>::Point SimplexQEF<0>::solve() const
{
    // We're returning a 0d point, there's no data there.
    return {};
}

template <unsigned N>
inline typename SimplexQEF<N>::Point SimplexQEF<N>::solve() const
{
    // Our masks here use bits from 0 to 2^(N+1)-2; any other bits are ignored.
    // Each such bit represents a number that is itself a mask of N+1 bits, 
    // corresponding to which constraints are set to 0 to form a subspace.
    // (Thus, there are 2^(N+1)-2 subspaces).
    static_assert(sizeof(unsigned) >= 2);
    unsigned currentMasks = 1;
    unsigned nextMasksAnyFailed = 0;
    unsigned nextMasksAllFailed = unsigned(-1);
    auto maskBitCount = 0; // For debugging only.
    while (currentMasks != 0) {
        for (auto i = 0; i < 8 * sizeof(unsigned); ++i) {
            if (1 << i & currentMasks) {
                assert(bitcount(i) == maskBitCount);
                if (i != 0) {
                    assert(highestbit(i) <= N + 1);
                }
            }
        }
        auto point = solveFromMask<0>(
            currentMasks, nextMasksAnyFailed, nextMasksAllFailed);
        if (point.has_value()) {
            return point.value();
        }

        currentMasks = nextMasksAnyFailed & nextMasksAllFailed;
        nextMasksAnyFailed = 0;
        nextMasksAllFailed = unsigned(-1);
        ++maskBitCount;
//        assert(maskBitCount <= N);
    }

    // We should not reach this point unless the geometric algorithm has
    // gone awry.
//    assert(false);
    return massPoint();
}

template<unsigned N>
unsigned SimplexQEF<N>::constraintFailures(const Point& pt) const
{
    auto out = 0;
    for (auto i = 0; i <= N; ++i) {
        const auto& constraint = constraints[i];
        auto constraintResult = constraints[i].val(pt);
        if (constraintResult >= 0) {
            out |= (1 << i);
        }
    }
    return out;
}

template <unsigned N>
typename SimplexQEF<N>::Vector SimplexQEF<N>::solveRaw() const
{
    // When the matrix rank is less than full so we have freedom in terms of
    // which point to use, we want the one closest to our "mass point" (average
    // of max-rank intersections).
    auto target = massPoint();

    // Our high-level goal here is to find the pseduo-inverse of AtA,
    // with special handling for when it isn't of full rank.
    Eigen::SelfAdjointEigenSolver<Matrix> es(AtA);
    auto eigenvalues = es.eigenvalues().real();

    // Build the SVD's diagonal matrix; since we'll constrain, there is no need
    // to truncate near-singular eigenvalues.
    Matrix D = Matrix::Zero();

    for (unsigned i = 0; i < N; ++i) {
        if (eigenvalues[i] != 0)
        {
            D.diagonal()[i] = 1 / eigenvalues[i];
        }
    }

    // SVD matrices
    auto U = es.eigenvectors().real().eval(); // = V

#ifdef LIBFIVE_VERBOSE_QEF_DEBUG
    std::cout << "Eigenvalues: [\n" << eigenvalues << "]\n";
    std::cout << "Eigenvectors: [\n" << U << "]\n";
#endif

    // Pseudo-inverse of A
    auto AtA_inv = (U * D * U.transpose()).eval();

    // Solve for vertex position (minimizing distance to target)
    return AtA_inv * (AtB - (AtA * target)) + target;
}

template <>
typename SimplexQEF<0>::Vector SimplexQEF<0>::solveRaw() const
{
    // There's no actual data here, we can just return a default-constructed
    // object.
    return SimplexQEF<0>::Vector{};
}

template class SimplexQEF<3>;

}   // namespace libfive
