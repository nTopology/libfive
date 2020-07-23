/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "simplex_dc_collapser.inl"
#include "simplex_dc_tree3.cpp"

namespace libfive {

template <>
void SimplexDCCollapser<3, 3>::load(Input* t)
{
    const auto& subs = t->leaf->sub;
    auto cellSub = subs[26].load();
    auto& cellVert = cellSub->vert;
    const auto& region = t->region;
    auto pos = 0;
    auto floating = 0;
    auto targetDimension = 0;
    for (auto i = 0; i < 3; ++i) {
        assert(region.lower(i) < region.upper(i));
        assert(region.lower(i) <= cellVert(i));
        assert(region.upper(i) >= cellVert(i));
        if (region.upper(i) == cellVert(i)) {
            pos |= 1 << i;
        }
        else if (region.lower(i) < cellVert(i)) {
            floating |= 1 << i;
            ++targetDimension;
        }
    }

    auto targetNeighbor = NeighborIndex::fromPosAndFloating(pos, floating);
    auto recalcCell = [&, this]()
    {
        // We can't use any intermediates because they all collapse
        // (or are split), but some collapse to something else 
        // (or are split), so we have to us the recalculator.
        SimplexDCRecalculator<3, 3> recalc(eval, t, 26, bound_cutoff);

        cellSub->vert = recalc.getNewVertex();
        switch (t->type) {
        case Interval::AMBIGUOUS:
        {
            std::array<bool, ipow(3, 3)> alreadySolved;
            alreadySolved.fill(true);
            alreadySolved[26] = false;
            t->saveVertexSigns(
                eval, eval->getDeck()->tape, alreadySolved);
            break;
        }
        case Interval::FILLED:
        case Interval::EMPTY:
            cellSub->inside = (t->type == Interval::FILLED);
            break;
        default:
            assert(false);
            cellSub->inside = false;
        }
        targetNeighbor = 26;
        cellSub->collapseRef = 26;
    };

    if (targetDimension == 1 || targetDimension == 2) {
        for (auto toFloat = 0; toFloat < 7; ++toFloat) {
            if (floating & ~toFloat) {
                continue;
            }
            auto testCell = t->neighbor(
                NeighborIndex::fromPosAndFloating(pos, toFloat));
            if (testCell && testCell->isBranch()) {
                // Our subspace splits, and there's no way we can
                // collapse to both (unless all the splits collapse
                // to the center, but that's too much of an edge case
                // to implement yet).
                recalcCell();
                return;
            }
        }
    }

    auto targetSub = subs[targetNeighbor.i].load();
    if (targetDimension < 3) {
        if (targetDimension > 0) {
            auto newTarget = targetNeighbor
                .fromRelativeToThis(targetSub->collapseRef);
            if (newTarget.i != targetNeighbor.i) {
                targetNeighbor = newTarget;
                floating = targetNeighbor.floating();
                pos = targetNeighbor.pos();
                targetDimension = bitcount(floating);
            }
        }
        auto checkIntermediates = [&](auto intermediates) {
            // intermediates is std::array<NeighborIndex, something>.

            // We want to check which (if any) intermediates are uncollapsing, 
            // and if there are any that are collapsing but not to 
            // targetNeighbor.  For these purposes, a "split" (a subspace that
            // is further subdivided due to a branch neighbor) is considered
            // collapsing.
            constexpr auto size = intermediates.size();
            std::array<bool, size> uncollapsing;
            std::array<bool, size> intermediatesDisqualify;
            intermediatesDisqualify.fill(false);
            bool hasUncollapsing = false;
            bool hasUnusableCollapsing = false;
            for (auto i = 0; i < size; ++i) {
                auto intermediate = intermediates[i];
                // First check if there's a neighbor that forces this 
                // intermediate subspace to split; if it splits, it's not
                // usable.
                auto floating = intermediate.floating();
                auto pos = intermediate.pos();
                auto isSplit = false;
                if (intermediate.dimension() != 0) {
                    for (auto toFloat = 0; toFloat < 8; ++toFloat) {
                        if (floating & ~toFloat) {
                            continue; // If the edge has it floating, so must
                                      // the neighbor direction.
                        }
                        auto neighbor = t->neighbor(
                            NeighborIndex::fromPosAndFloating(pos, toFloat));
                        if (neighbor && neighbor->isBranch()) {
                            // The intermediate subspace is split.
                            hasUnusableCollapsing = true;
                            isSplit = true;
                            uncollapsing[i] = false;
                            break;
                        }
                    }
                }
                if (isSplit) {
                    continue;
                }
                auto intermediateSub = subs[intermediate.i].load();
                auto intermediateCollapse = 
                    intermediate.dimension() == 0
                    ? intermediate
                    : intermediate.fromRelativeToThis(
                    intermediateSub->collapseRef);
                if (intermediateCollapse.i == intermediate.i) {
                    uncollapsing[i] = true;
                    hasUncollapsing = true;
                }
                else {
                    uncollapsing[i] = false;
                    if (intermediateCollapse.i != targetNeighbor.i) {
                        hasUnusableCollapsing = true;
                    }
                }
                // Either way, if we've got size of 6 and this is a face, it 
                // disqualifies the edges "lower down" except for any it 
                // matches.
                if constexpr (size == 6) {
                    if (intermediates[i].dimension() == 2) {
                        for (auto j = 0; j < size; ++j) {
                            if (i == j) {
                                continue;
                            }
                            if (intermediates[j].floating() & ~floating) {
                                continue;
                            }
                            assert(intermediates[j].dimension() == 1);
                            if (intermediateCollapse.i == intermediates[j].i) {
                                continue;
                            }
                            // We don't need to worry about a case where 
                            // intermediates[j] is not intermediateCollapse but
                            // does collapse to it; if it collapses, we're not
                            // using it anyway.
                            intermediatesDisqualify[j] = true;
                        }
                    }
                }
            }
            if (hasUncollapsing) {
                // We can't use our current targetNeighbor, since there is
                // an intermediate space that doesn't collapse to it, but
                // we can use an uncollapsing intermediate subspace as our 
                // collapse target.
                int newTargetIdx;
                static_assert(size == 2 || size == 6);
                /// We would like to capture uncollapsing and 
                /// intermediatesDisqualify in the following lambda, but Visual
                /// Studio 2019 throws a compile error when attempting to do so
                /// (presumably a compiler bug), so as a workaround we get their 
                /// addresses, capture those, and dereference them.
                auto uncPtr = &uncollapsing;
                auto disqPtr = &intermediatesDisqualify;
                auto compare = [&, uncPtr, disqPtr](int a, int b) {
                    if (!(*uncPtr)[a] || (*disqPtr)[a]) {
                        return false;
                    }
                    else if (!(*uncPtr)[b] || (*disqPtr)[b]) {
                        return true;
                    }
                    else {
                        return compareCollapseCandidates(
                            t, *cellSub, intermediates[a].i, intermediates[b].i);
                    }
                };
                if constexpr (size == 2) {
                    newTargetIdx = std::min(0, 1, compare);
                    targetNeighbor = intermediates[newTargetIdx];
                }
                else {
                    newTargetIdx = std::min({ 0, 1, 2, 3, 4, 5 }, compare);
                    if (intermediatesDisqualify[newTargetIdx]) {
                        // This can happen if the only uncollapsing elements are
                        // edges disqualified by intermediate faces.
                        for (auto i = 0; i < 6; ++i) {
                            assert(intermediatesDisqualify[i] ||
                                   !uncollapsing[i]);
                        }
                        recalcCell();
                    }
                    else {
                        targetNeighbor = intermediates[newTargetIdx];
                    }
                }
            }
            else if (hasUnusableCollapsing) {
                recalcCell();
            }
            // Otherwise, everything collapses to our target, so we're good.
        };
        if (targetDimension == 1) {
            std::array<NeighborIndex, 2> intermediates;
            auto currentIdx = 0;
            for (auto i = 0; i < 3; ++i) {
                auto faceAxis = 1 << i;
                if (floating & faceAxis) {
                    continue;
                }
                auto toFloat = 7 & ~(1 << i);
                assert(currentIdx < 2);
                intermediates[currentIdx++] = 
                    NeighborIndex::fromPosAndFloating(pos, toFloat);
            }
            checkIntermediates(intermediates);
        }
        else if (targetDimension == 0) {
            std::array<NeighborIndex, 6> intermediates;
            for (auto i = 0; i < 6; ++i) {
                auto toFloat = 1 + i;
                intermediates[i] =
                    NeighborIndex::fromPosAndFloating(pos, toFloat);
            }
            checkIntermediates(intermediates);
        }
    }
    cellSub->collapseRef = targetNeighbor.i;
}

template class SimplexDCCollapser<3, 3>;
}   // namespace libfive
