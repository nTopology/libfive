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
SimplexQEF<N>::SimplexQEF(Eigen::Matrix<double, N, N + 1>&& vertices) :
    AtA(Matrix::Zero()), AtB(Vector::Zero()), BtB(0.)
{
    if constexpr (N == 2 || N == 3) {
        auto center = vertices.rowwise().mean();
        auto diff = (vertices.colwise() - center).eval();
        vertices -= constraintPadding * diff;

        // Our formula to produce a constraint only gives us where it is 0; to
        // determine which side is negative, we need to test one of them.  We
        // choose to test constraint N + 1 (the one associated with vertices
        // from 0 to N), since that will ideally correspond to an axis-aligned
        // side of the simplex and therefore minimize the effect of numeric
        // error.
        auto negateAtN = false;
        for (auto i = N; i >= 0; ++i) {
            auto shiftedVert = [&](int idx)
            {return vertices.col((i + idx) % (N + 1)); };
            auto diff1 = shiftedVert(2) - shiftedVert(1);
            Vector perp;
            if constexpr (N == 2) {
                perp << diff1.y(), -diff1.x();
            }
            else {
                auto diff2 = shiftedVert(3) - shiftedVert(1);
                perp = diff1.cross(diff2);
            }
            auto val = perp.dot(shiftedVert(1));
            assert(val != perp.dot(vertices.col(i))); // If it is, our simplex
                                                      // is degenerate.
            if (i == N && val < perp.dot(vertices.col(N))) {
                // perp causes vertices[N] to be higher than our constant
                // offset, so after the offset it'll still be positive, so
                // we need to negate.
                negateAtN = true;
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
    else {
        // QEFs with N of 0 or 1 should be occurring only as part of a higher-
        // dimensional QEF, and thus should not be using this constructor at 
        // all.
        assert(false);
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
            out.emplace();
            auto& outVal = out.value();
            const auto& constraint = constraints[highBit];
            auto axis = constraint.strongestDimension();
            assert(axis < N);

            // Get a "reduced" constraint matrix that we will use to adjust
            // the values of outVal's members
            const auto& matrix = constraint.matrix;
            auto normalized = matrix / matrix[axis];
            Eigen::Matrix<double, N, 1> reduced;
            reduced << normalized.head(axis), normalized.tail(N - axis);
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
                outVal.AtB += reducedFromAtA.transpose() * reduced[N];
                outVal.AtB -= AtA(axis, axis) * reducedHead * reduced[N];
                outVal.AtB = outVal.AtB.eval();
            }

            outVal.BtB = BtB; 
            outVal.BtB += AtA(axis, axis) * reduced[N] * reduced[N];
            outVal.BtB += 2 * AtB[axis] * reduced[N];

            for (auto idx = 0, newIdx = 0; idx < constraints.size(); ++idx) {
                if (idx == highBit) {
                    continue;
                }
                auto& matrix = outVal.constraints[newIdx].matrix;
                const auto& oldMatrix = constraints[idx].matrix;
                matrix << oldMatrix.head(axis), oldMatrix.tail(N - axis);
                matrix -= oldMatrix[N] * reduced;
                ++newIdx;
            }

            outVal.mass_point << 
                mass_point.head(axis), mass_point.tail(N - axis);
            outVal.mass_point -= mass_point[N] * reduced;

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
    unsigned& nextMasksAnyValid, 
    unsigned& nextMasksAllValid) const
{
    if ((1 << current_mask) & currentMasks) {
        auto rawResult = sub<current_mask>().solveRaw();
        auto result = convertFromSub<current_mask>(rawResult);
        auto failedConstraint = false;
        for (auto i = 0; i <= N; ++i) {
            auto shifted = 1 << i;
            if (shifted & current_mask) {
                // Don't check the constraints that define this subspace.
                continue;
            }
            const auto& constraint = constraints[i];
            auto constraintResult = constraints[i].val(result);
            auto constraintMask = 1 << (current_mask | shifted);
            if (constraintResult > 0) {
                assert(bitcount(current_mask) < N);
                failedConstraint = true;
                nextMasksAllValid &= ~constraintMask;
            }
            else {
                nextMasksAnyValid |= constraintMask;
            }
        }
        if (!failedConstraint) {
            return result;
        }
    }
    constexpr auto next_mask = current_mask + 1;
    if constexpr (bitcount(next_mask) > N) {
        return {};
    }
    else {
        return solveFromMask<next_mask>(
            currentMasks, nextMasksAnyValid, nextMasksAllValid);
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
    unsigned nextMasksAnyValid = 0;
    unsigned nextMasksAllValid = unsigned(-1);
    auto maskBitCount = 1; // For debugging only.
    while (currentMasks != 0) {
        for (auto i = 0; i < 8 * sizeof(unsigned); ++i) {
            if (1 << i & currentMasks) {
                assert(bitcount(i) == maskBitCount);
                assert(highestbit(i) <= N + 1);
            }
        }
        auto point = solveFromMask<0>(
            currentMasks, nextMasksAnyValid, nextMasksAllValid);

        if (point.has_value()) {
            return point.value();
        }
        currentMasks = nextMasksAnyValid & nextMasksAllValid;
        nextMasksAnyValid = 0;
        nextMasksAllValid = unsigned(-1);
        ++maskBitCount;
        assert(maskBitCount <= N);
    }

    // We should not reach this point unless the geometric algorithm has
    // gone awry.
    assert(false);
    return massPoint();
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

    for (unsigned i = 0; i < N + 1; ++i) {
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
