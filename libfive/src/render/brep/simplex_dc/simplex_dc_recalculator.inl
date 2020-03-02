/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "libfive/render/brep/simplex_dc/simplex_dc_recalculator.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "libfive/render/brep/indexes.hpp"

namespace libfive {

template<unsigned N, unsigned K>
SimplexDCRecalculator<N, K>::SimplexDCRecalculator(
    Evaluator* eval, SimplexDCTree<N>* tree, int targetSub, double boundCutoff)
    : eval(eval), tape(tree->leaf->tape), 
    is_ambiguous(tree->type == Interval::AMBIGUOUS),
    bound_cutoff(boundCutoff)
{
    assert(NeighborIndex(targetSub).dimension() == K);

    QEF<K>(*getQEFSub)(const QEF<N>&);
    Region<K>(*getRegionSub)(const Region<N>&);
    const auto& rawRegion = tree->region;
    if constexpr (N == K) {
        missing_dim = N == 3 ? -1 : 2;
        if constexpr (N == 2) {
            missing_dim_value = tree->region.perp[0];
        }
        getQEFSub = [](const QEF<N>& qef) {return qef; };
        getRegionSub = [](const Region<N>& region) {return region; };
    }
    else {
        Eigen::Vector3d regionCorner;
        switch (targetSub) {
        case 8:
        case 20:
        case 24:
            regionCorner = rawRegion.lower;
            break;
        case 17:
        case 23:
        case 25:
            regionCorner = rawRegion.upper;
            break;
        default:
            assert(false);
            return;
        }
        switch (targetSub) {
        case 8:
        case 17:
            missing_dim = 2;
            missing_dim_value = regionCorner[2];
            getQEFSub = [](const QEF<N>& qef) {return qef.template sub<3>(); };
            getRegionSub = [](const Region<N>& region) 
                {return region.template subspace<3>(); };
            break;
        case 20:
        case 23:
            missing_dim = 1;
            missing_dim_value = regionCorner[1];
            getQEFSub = [](const QEF<N>& qef) {return qef.template sub<5>(); };
            getRegionSub = [](const Region<N>& region) 
                {return region.template subspace<5>(); };
            break;
        case 24:
        case 25:
            missing_dim = 0;
            missing_dim_value = regionCorner[0];
            getQEFSub = [](const QEF<N>& qef) {return qef.template sub<6>(); };
            getRegionSub = [](const Region<N>& region) 
                {return region.template subspace<6>(); };
            break;
        default:
            assert(false);
            return;
        }
    }

    total_region = getRegionSub(rawRegion);

    for (auto i = 0; i < ipow(3, N); ++i) {
        if (NeighborIndex(targetSub).contains(NeighborIndex(i))) {
            if (i != targetSub) {
                assert(NeighborIndex(i).dimension() < K);
                addBorderSubs(tree, i);
            }
            summed_qef += getQEFSub(tree->leaf->sub[i].load()->qef);
        }
    }
}

template<unsigned N, unsigned K>
Eigen::Matrix<double, N, 1> SimplexDCRecalculator<N, K>::getNewVertex()
{
    current_vertex = total_region.center();
    // We only need to check border subs if this cell is ambiguous; otherwise,
    // nothing's going to happen in the tets inside it anyway.
    if (is_ambiguous) {
        Region<K> region(current_vertex, current_vertex);

        auto toContinue = true;
        while (toContinue) {
            toContinue = false;
            auto nextRegion = region;
            for (auto sub : border_subs) {
                Eigen::Matrix<double, K, 1> extremum;
                if (getLineExtremum(sub, extremum) &&
                    !region.contains(extremum, 0)) {
                    toContinue = true;
                    nextRegion.lower = nextRegion.lower.min(extremum.array());
                    nextRegion.upper = nextRegion.upper.max(extremum.array());
                }
            }
            auto upperDiff = total_region.upper - nextRegion.upper;
            auto lowerDiff = nextRegion.lower - total_region.lower;
            auto upperRatios = upperDiff / (total_region.upper - region.upper);
            auto lowerRatios = lowerDiff / (region.lower - total_region.lower);

            for (auto i = 0; i < K; ++i) {
                assert(upperRatios[i] <= 1.);
                assert(lowerRatios[i] <= 1.);
                assert(upperRatios[i] > 0.);
                assert(lowerRatios[i] > 0.);
                if (upperRatios[i] < 1. && upperRatios[i] >= 0.99) {
                    // It's converging too slowly, let's speed it up a bit.
                    nextRegion.upper(i) += 0.1 * upperDiff[i];
                }
                if (lowerRatios[i] < 1. && lowerRatios[i] >= 0.99) {
                    nextRegion.lower(i) -= 0.1 * lowerDiff[i];
                }
            }
            region = nextRegion;
            current_vertex = summed_qef
              .solveBounded(region, 1., bound_cutoff).position;
        }
    }
    if constexpr (N == K) {
        return current_vertex;
    }
    else {
        static_assert(N == 3);
        static_assert(K == 2);
        Eigen::Matrix<double, N, 1> out;
        out.head(missing_dim) = current_vertex.head(missing_dim);
        out[missing_dim] = missing_dim_value;
        out.tail(2 - missing_dim) = current_vertex.tail(2 - missing_dim);
        return out;
    }
}

template<unsigned N, unsigned K>
inline bool SimplexDCRecalculator<N, K>::getLineExtremum(
    SimplexLeafSubspace<N>* borderSub, 
    Eigen::Matrix<double, K, 1>& result) const
{
    const auto& borderPtFull = borderSub->vert;
    Eigen::Matrix<double, K, 1> borderPtReduced;
    auto projectK = [this](const Eigen::Vector3d& vec) {
        if constexpr (K == 3) {
            return vec;
        }
        else {
            static_assert(K == 2);
            assert(vec[missing_dim] == missing_dim_value);
            Eigen::Matrix<double, K, 1> out;
            out.head(missing_dim) = vec.head(missing_dim);
            out.tail(2 - missing_dim) = vec.tail(2 - missing_dim);
            return out;
        }
    };

    if constexpr (N == 3) {
        borderPtReduced = projectK(borderPtFull);
    }
    else {
        borderPtReduced = borderPtFull;
    }
    auto diff = (borderPtReduced - current_vertex).template cast<float>();

    auto project1 = [&diff, this](const Eigen::Vector4f& grad)
    {
        Eigen::Matrix<float, K, 1> reduced;
        if constexpr (N == K) {
            reduced = grad.template head<N>();
        }
        else {
            auto reduced3 = grad.template head<3>();
            reduced.head(missing_dim) = reduced3.head(missing_dim);
            reduced.tail(2 - missing_dim) = reduced3.tail(2 - missing_dim);
        }
        return reduced.dot(diff);
    };

    auto extend = [this](const Eigen::Matrix<double, K, 1>& reduced) {
        if constexpr (K == 3) {
            return reduced;
        }
        else {
            Eigen::Vector3d out;
            out.head(missing_dim) = reduced.head(missing_dim);
            out[missing_dim] = missing_dim_value;
            out.tail(2 - missing_dim) = reduced.tail(2 - missing_dim);
            return out;
        }
    };
    auto extendF =
        [&extend](const Eigen::Matrix<double, K, 1>& reduced)->Eigen::Vector3f {
        return extend(reduced).template cast<float>();
    };
    Eigen::Vector3d borderPoint3;
    if constexpr (N == 3) {
        assert(extend(borderPtReduced) == borderPtFull.transpose());
        borderPoint3 = borderPtFull;
    }
    else {
        borderPoint3 = extend(borderPtFull);
    }

    // We are currently evaluating only one intersection at a time.  Later, we
    // will likely want to do more (as there are multiple subs we can compare 
    // to), and will want to loop over them.

    // First, we evaluate at our endpoints; if they result in the same sign or,
    // either is 0, we're done.
    eval->set(extendF(current_vertex), 0);
    eval->set(borderPoint3.template cast<float>(), 1);
    auto out = eval->derivs(2, *tape);
    auto currentVal = project1(out.col(0));
    auto borderVal = project1(out.col(1));

    // The calculation of our result is now very similar to that used elsewhere
    // (e.g. simplex_dc_intersecter.inl) for bisection, except based on the
    // derivative rather than the original.

    // We are currently evaluating only one boundary subspace vertex at a time.
    // Later, we will likely want to do more, and will want to loop over them;
    // a bit of the loop code is preserved from the DC code, but a lot will 
    // still need to be changed to make it functional.

    constexpr int SEARCH_COUNT = 4;
    constexpr int POINTS_PER_SEARCH = 16;
    static_assert(POINTS_PER_SEARCH <= DerivArrayEvaluator::N,
        "Overflowing DerivArrayEvaluator data array");

    constexpr auto eval_count = 1;

    std::array<std::pair<Eigen::Vector3d, Eigen::Vector3d>, eval_count> targets;
    if (currentVal < 0 && borderVal > 0) {
        targets[0].first = extend(current_vertex);
        targets[0].second = borderPoint3;
    }
    else if (currentVal > 0 && borderVal < 0) {
        targets[0].first = borderPoint3;
        targets[0].second = extend(current_vertex);
    }
    else {
        // early return
        return false;
    }

    // Multi-stage binary search for intersection
    for (int s = 0; s < SEARCH_COUNT; ++s)
    {
        // Load search points into the evaluator
        Eigen::Array<double, 3, POINTS_PER_SEARCH> ps;
        for (int j = 0; j < POINTS_PER_SEARCH; ++j)
        {
            const double frac = j / (POINTS_PER_SEARCH - 1.0);
            ps.col(j) = (targets[0].first * (1 - frac)) +
                (targets[0].second * frac);
            eval->set(ps.col(j).template cast<float>(), j);
        }

        auto out = eval->derivs(POINTS_PER_SEARCH, *tape);

        // Skip one point, because the very first point is
        // already known to have a negative projection (but
        // sometimes, due to numerical issues, it registers
        // as positive!)
        for (unsigned j = 1; j < POINTS_PER_SEARCH; ++j)
        {
            // We're searching for the first point that has a positive (or 0)
            // value.  There's a special case for the final point in the
            // search, working around  numerical issues where different
            // evaluators disagree with whether points are inside or outside.
            auto projected = project1(out.col(j));
            auto validateResult = [this, &result]() {
                for (auto i = 0; i < K; ++i) {
                    if (result(i) == total_region.lower(i) ||
                        result(i) == total_region.upper(i)) {
                        return false;
                    }
                }
                return true;
            };
            if (projected == 0) {
                // We found our result, and can return early.
                result = projectK(ps.col(j).matrix().eval());
                return validateResult();
            }
            else if (projected > 0 || j == POINTS_PER_SEARCH - 1)
            {
                if (s == SEARCH_COUNT - 1) {
                    // This is our best result.
                    auto rawResult = ps.col(j).matrix().eval();
                    if (projected > 0) {
                        rawResult += ps.col(j - 1).matrix();
                    }
                    else {
                        assert(j == POINTS_PER_SEARCH - 1);
                        rawResult += targets[0].second;
                    }
                    rawResult /= 2;

                    result = projectK(rawResult.eval());
                    return validateResult();
                }
                else {
                    targets[0] = { ps.col(j - 1), ps.col(j) };
                    break;
                }
            }
        }
    }
    // We should never hit this point, as if we reach the last run of the inner
    // loop, it should always either return or hit the line checking if 
    // s == SEARCH_COUNT - 1, and on the last run of the outer loop that returns
    // true, so we return regardless.
    assert(false);
    return false;
}

template<unsigned N, unsigned K>
inline void SimplexDCRecalculator<N, K>::addBorderSubs(
    SimplexDCTree<N>* tree, int targetSub)
{
    NeighborIndex targetNeighbor(targetSub);
    if (tree->isBranch()) {
        auto targetFloating = targetNeighbor.floating();
        auto targetPos = targetNeighbor.pos();
        for (auto i = 0; i < ipow(3, N); ++i) {
            NeighborIndex childIndex(i);
            auto childFloating = childIndex.floating();
            if ((childFloating | targetFloating) != (1 << N) - 1) {
                // There is a dimension in which targetNeighbor is fixed,
                // and our child index is fixed (meaning we want it to be
                // on one side or the other); in any dimension where
                // targetNeighbor is fixed, the child must cover the "entire"
                // span of targetNeighbor, so it can't be fixed in that 
                // dimension.
                continue;
            }
            // Determine which child to use; if either the child or the target 
            // neighbor is fixed in a dimension, we want to use it.  If both are
            // floating, we use the lower child.  (It is not possible to reach
            // this point if both are fixed.)
            auto cellChild = childIndex.pos() | targetPos;
            // To calculate the new sub for this child: In each dimension that
            // is not floating, targetNeighbor must be floating to reach this
            // point, and we want to keep it floating.  In each dimension that
            // is floating, if targetNeighbor was fixed, we must keep it (and
            // its pos) so that it remains inside targetNeighbor.  If both were
            // floating, that represents a lower-dimensional interior child of
            // our original subspace, so we want to set our new sub to fixed,
            // with upper position.  Thus (exploiting the fact that in a 
            // dimension that is floating, the pos passed is irrelevant):

            auto newPos = targetFloating | targetPos;
            auto newFloat = targetFloating & ~childFloating;

            auto newSub = NeighborIndex::fromPosAndFloating(newPos, newFloat);
            addBorderSubs(tree->child(cellChild), newSub.i);
        }
    }
    else {
        if (targetNeighbor.dimension() > 0) {
            // The target subspace may branch even if 'tree' does not, if
            // one of the neighbors containing it branches.
            constexpr auto maxBits = (1 << N) - 1;
            for (auto toFloat = 0; toFloat < maxBits - 1; ++toFloat) {
                if (targetNeighbor.floating() & ~toFloat) {
                    // If a dimension is floating in the subspace we're 
                    // checking, any neighbors containing that subspace must 
                    // likewise have that dimension floating in their direction.
                    continue;
                }
                auto testDirection = NeighborIndex::fromPosAndFloating(
                    targetNeighbor.pos(), toFloat);
                auto testNeighbor = tree->neighbor(testDirection);
                if (testNeighbor && testNeighbor->isBranch()) {
                    // Use it instead, adjusting our target subspace 
                    // accordingly.
                    auto newSub = NeighborIndex::fromPosAndFloating(
                        targetNeighbor.pos() ^ toFloat ^ maxBits,
                        targetNeighbor.floating());
                    addBorderSubs(testNeighbor, newSub.i);
                    return;
                }
            }
        }
        // If we've reached this point, our subspace does in fact not branch.
        // Therefore, we want to add the subspace, unless its vertex collapses
        // to another vertex.  How we determine this depends on its dimension.
        auto sub = tree->leaf->sub[targetSub].load();
        switch (targetNeighbor.dimension()) {
        case 0:
            // It's always minimal.
            break;
        case 1:
        {
            auto vert = sub->vert;
            const auto& region = tree->region;
            for (auto i = 0; i < N; ++i) {
                if (!targetNeighbor.isAxisFixed(i)) {
                    if (vert[i] == region.lower[i] ||
                        vert[i] == region.upper[i]) {
                        // It either did collapse, or will collapse; we don't
                        // need to use it.
                        return;
                    }
                    break; // There is exactly one non-fixed axis of a 
                           // dimension-1 neighbor
                }
            }
            break; // From the switch.
        }
        case 2:
            // targetSub should always be of dimensionality lower than K, which
            // should be at most 3.
            assert(K == 3);
            if (sub->collapseRef != 8) {
                // It collapses.
                assert(sub->collapseRef < 8);
                return;
            }
            break;
        default:
            // Again, K must have dimensionality at most 3 and more than the
            // dimensionality of targetSub, so a dimensionality higher than 2
            // should not be possible.
            assert(false);
            return;
        }
        border_subs.push_back(sub);
    }
}


}   // namespace libfive
