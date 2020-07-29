/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "libfive/render/brep/simplex_dc/simplex_dc_intersecter.hpp"

#include "libfive/render/brep/simplex_dc/simplex_dc_circulator.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "libfive/render/brep/object_pool.hpp"
#include "libfive/render/brep/indexes.hpp"
#include "libfive/render/brep/per_thread_brep.hpp"
#include <set>

namespace libfive {

template <unsigned N>
SimplexDCIntersecter<N>::SimplexDCIntersecter(
    PerThreadOutput& m, Tree t, ObjectPool<SimplexDCIntersection<N>>& pool, 
    double mergeRatioSquared, Perp perp)
    : parent_pool(pool), eval(new Evaluator(t)), owned(true), 
    perp(perp), m(m), merge_ratio_squared(mergeRatioSquared)
{
    // Nothing to do here
}

template <unsigned N>
SimplexDCIntersecter<N>::SimplexDCIntersecter(
    PerThreadOutput& m, Evaluator* es, ObjectPool<SimplexDCIntersection<N>>& pool, 
    double mergeRatioSquared, Perp perp)
    : parent_pool(pool), eval(es), owned(false), 
    perp(perp), m(m), merge_ratio_squared(mergeRatioSquared)
{
    // Nothing to do here
}

template <unsigned N>
SimplexDCIntersecter<N>::~SimplexDCIntersecter() {
    if (owned) {
        delete eval;
    }
    parent_pool.claim(my_pool);
}

template <unsigned N>
template <Axis::Axis A>
void SimplexDCIntersecter<N>::load(const std::array<Input*, 1 << (N - 2)>& ts)
{
    if constexpr (N == 3) {
        // Skip this if all of the cells are empty / filled.
        if (std::all_of(ts.begin(), ts.end(),
                        [](const Input* t) 
                            { return t->type != Interval::AMBIGUOUS; }))
        {
            return;
        }

        // We want to use the smallest cell adjoining this face.
        const auto index = std::min_element(ts.begin(), ts.end(),
            [](const Input* a, const Input* b)
        { return a->leafLevel() < b->leafLevel(); }) - ts.begin();

        auto faceSubspaceIndex = 
            (ipow(3, N) - 1) - ipow(3, Axis::toIndex(A)) * (index + 1);

        unsigned trueFaceIdx;
        auto faceSub = ts[index]->leaf->collapsedSub(faceSubspaceIndex, 
                                                     &trueFaceIdx);
        if (trueFaceIdx != faceSubspaceIndex) {
            assert(trueFaceIdx < faceSubspaceIndex);
            // If the face sub collapses to a lower dimension, the face
            // cannot be a minimal dimension for any sub containing it;
            // we'll handle this pair (or already did, or are currently
            // doing so) on the minimal-dimension pair.  If it does not,
            // then it is the minimal pair whenever the cell sub does not
            // collapse to the face sub.  (The cell sub is not allowed
            // to collapse to an edge or corner sub in the face sub if
            // the face sub does not collapse as well.)
            return;
        }

        // Now that we have the intersections, we need to write them to every 
        // simplex that contains this face.  Again, we use the smallest of the 
        // associated trees.
        auto& edges = ts[index]->leaf->edges;

        for (auto cell = 0; cell < 2; ++cell) {
            if (ts[cell]->type != Interval::AMBIGUOUS) {
                continue;
            }
            const auto& cellSub = ts[cell]->leaf->collapsedSub(ipow(3, N) - 1);

            auto isSame = faceSub == cellSub;
            if (!isSame) {
                if (faceSub->inside == cellSub->inside) {
                    // Not an active edge.
                    continue;
                }
            }
            auto inside = faceSub->inside ? faceSub : cellSub;
            auto outside = faceSub->inside ? cellSub : faceSub;
            auto intersection = isSame ? &DCSimplex<N>::dupVertIntersection
                : searchEdge(inside, outside, ts[cell]->leaf->tape);
            auto cellSubIdx = (cell == index ) ? ipow(3, N) - 1 
                                               : ipow(3, N) + faceSubspaceIndex;
            auto insideIdx = faceSub->inside ? faceSubspaceIndex : cellSubIdx;
            auto outsideIdx = faceSub->inside ? cellSubIdx : faceSubspaceIndex;
            SimplexDCCirculator<N, false> circulator(
                ts[index], insideIdx, outsideIdx, isSame);
            while (!circulator.isEnd()) {
                auto simplex = circulator.simplex();
                if (!simplex) {
                    assert(false);
                    break;
                }
                circulator.insertIntersection(intersection);
                ++circulator;
            }
        }
    }
    else {
        // We should never be calling this version for N == 2 (though the code
        // still references N where feasible).
        assert(false); 
    }
}

template <unsigned N>
template <Axis::Axis A>
void SimplexDCIntersecter<N>::load(const std::array<Input*, 1 << (N - 1)>& ts)
{
    // Skip this if all of the cells are empty / filled.
    if (std::all_of(ts.begin(), ts.end(),
        [](const Input * t)
    { return t->type != Interval::AMBIGUOUS; }))
    {
        return;
    }

    // We want to use the smallest cell adjoining this face.
    const auto index = std::min_element(ts.begin(), ts.end(),
        [](const Input * a, const Input * b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();

    // Calculate the subspace for the edge.
    std::array<int, 3> reductions{ 0, 0, 0 };
    reductions[(Axis::toIndex(A) + 1) % N] = 1 + (index & 1);
    if (N == 3) {
        reductions[(Axis::toIndex(A) + 2) % 3] = 1 + bool(index & 2);
    }

    auto totalReductions = reductions[0] + 
                           3 * reductions[1] + 
                           9 * bool(N == 3) * reductions[2];

    auto edgeSubspaceIndex = (ipow(3, N) - 1) - totalReductions;

    unsigned trueEdgeIdx;
    auto edgeSub = ts[index]->leaf->collapsedSub(
        edgeSubspaceIndex, &trueEdgeIdx);
    if (trueEdgeIdx != edgeSubspaceIndex) {
        assert(trueEdgeIdx < edgeSubspaceIndex);
        // As in the cell/face case, this is a condition to determine if this
        // is the minimal pair, and the only such condition when the other
        // sub's dimension is 2.
        return;
    }

    auto& edge = ts[index]->leaf->edgeFromReduced(A, ts.size() - 1 - index);

    if constexpr (N == 3) {
        // Handle edges between face vertices and our edge vertex.
        for (auto faceAxisIdx = 0; faceAxisIdx < 2; ++faceAxisIdx) {
            auto faceAxis = (faceAxisIdx == 0) ? Axis::Q(A) : Axis::R(A);
            for (auto facePosition = 0; facePosition < 2; ++facePosition) {
                auto indexA = facePosition == 0 ? 0 : 2 >> faceAxisIdx;
                auto indexB = indexA ^ (1 << faceAxisIdx);
                if (ts[indexA] == ts[indexB]) {
                    // Due to a merged cell, there is no face in that direction
                    // (not even one extending beyond this edge).
                    continue;
                }
                auto faceIndex = std::min(indexA, indexB,
                    [&](int indexA, int indexB) 
                {return ts[indexA]->leafLevel() < ts[indexB]->leafLevel(); });
                if (ts[faceIndex]->type == Interval::UNKNOWN) {
                    // This can happen if the edge was on the border of the
                    // bounding box and the face is toward the outside, so
                    // the only non-Interval::UNKNOWN cell(s) don't border the face
                    // in question.
                    continue;
                }
                auto indexFromFace = bool(faceIndex & (1 << faceAxisIdx));
                auto faceSubIndex = (ipow(3, N) - 1) - 
                    ipow(3, Axis::toIndex(faceAxis)) * (indexFromFace + 1);
                auto faceSub = ts[faceIndex]->leaf->collapsedSub(faceSubIndex);
                auto isSame = edgeSub == faceSub;
                if (!isSame) {
                    if (faceSub->inside == edgeSub->inside) {
                        // Not an active edge.
                        continue;
                    }
                }
                auto inside = edgeSub->inside ? edgeSub : faceSub;
                auto outside = edgeSub->inside ? faceSub : edgeSub;
                auto intersection = isSame ? &DCSimplex<N>::dupVertIntersection
                    : searchEdge(inside, outside, ts[faceIndex]->leaf->tape);
                auto circulator = [&]() {
                    if (faceIndex == index) {
                        auto insideIdx =
                            edgeSub->inside ? edgeSubspaceIndex : faceSubIndex;
                        auto outsideIdx =
                            edgeSub->inside ? faceSubIndex : edgeSubspaceIndex;
                        return SimplexDCCirculator<N, false>(
                            ts[index], insideIdx, outsideIdx, isSame);
                    }
                    else {
                        assert((faceIndex ^ index) != (1 << faceAxisIdx));
                        return SimplexDCCirculator<3, false>(
                            SimplexDCCirculator<3, true>::FarFaceTag{},
                            ts[index], edgeSubspaceIndex, A, faceAxis,
                            edgeSub->inside, isSame);
                    }
                }();
                while (!circulator.isEnd()) {
                    auto simplex = circulator.simplex();
                    if (!simplex) {
                        assert(false);
                        break;
                    }
                    circulator.insertIntersection(intersection);
                    ++circulator;
                }
            }
        }
    }

    // Now handle edges between cell vertices and our edge vertex.
    for (auto cell = 0; cell < ts.size(); ++cell) {
        if (ts[cell]->type == Interval::UNKNOWN) {
            continue;
        }
        if constexpr (N == 3) {
            auto opposite = index ^ 3;
            if (cell != opposite && ts[cell] == ts[opposite]) {
                assert(cell != index);
                // We have a duplicate cell, and we're handling it when 
                // cell == opposite.  The circulator will cover all the 
                // simplices around this pair of sub vertices, so there's
                // nothing more to do.
                continue;
            }
        }
        auto cellSub = ts[cell]->leaf->collapsedSub(ipow(3, N) - 1);
        auto isSame = edgeSub == cellSub;
        if (!isSame) {
            if (cellSub->inside == edgeSub->inside) {
                // Not an active edge.
                continue;
            }
            if (N == 3) {
                // Unlike in the cell/face and face/edge cases, and cell/edge in
                // 2D, checking if the lower-dimensional subspace collapses is 
                // not sufficient to determine that this is the minimal pair for
                // cell/edge in 3D. If the cell and edge share an axis that is
                // not the edge axis, and the edge is on a boundary of the cell
                // in that dimension, that means that the cell collapses to a 
                // face containing the edge, and so that face/edge is the 
                // minimal pair, not this one.

                auto isMinimalPair = true;
                for (auto testAxis = 0; testAxis < N; ++testAxis) {
                    if (testAxis == Axis::toIndex(A)) {
                        continue;
                    }
                    auto vertValue = edgeSub->vert(testAxis);
                    if (vertValue != ts[cell]->region.lower(testAxis) &&
                        vertValue != ts[cell]->region.upper(testAxis)) {
                        continue;
                    }
                    if (vertValue == cellSub->vert(testAxis)) {
                        isMinimalPair = false;
                        break;
                    }
                }
                if (!isMinimalPair) {
                    continue;
                }
            }
        }
        auto inside = edgeSub->inside ? edgeSub : cellSub;
        auto outside = edgeSub->inside ? cellSub : edgeSub;
        auto intersection = isSame ? &DCSimplex<N>::dupVertIntersection :
            searchEdge(inside, outside, ts[cell]->leaf->tape);

        auto cellSubIdx = ipow(3, N) - 1;
        if (cell != index) {
            if (N == 2) {
                cellSubIdx = edgeSubspaceIndex + 9;
            }
            else {
                auto diff = cell ^ index;
                switch (diff) {
                case 1:
                case 2:
                {
                    // For both of these, our cell is off of index by a face 
                    // NeighborIndex.
                    auto faceAxisIdx = (Axis::toIndex(A) + diff) % 3;
                    bool indexIsHigh(index & diff);
                    auto toSubtract = ipow(3, faceAxisIdx) * (indexIsHigh + 1);
                    cellSubIdx = (ipow(3, N) - 1) - toSubtract;
                    break;
                }
                case 3:
                    cellSubIdx = edgeSubspaceIndex;
                    break;
                default:
                    assert(false);
                    continue;
                }
                cellSubIdx += 27;
            }
        }
        auto insideIdx = edgeSub->inside ? edgeSubspaceIndex : cellSubIdx;
        auto outsideIdx = edgeSub->inside ? cellSubIdx : edgeSubspaceIndex;
        SimplexDCCirculator<N, false> circulator(
            ts[index], insideIdx, outsideIdx, isSame);
        while (!circulator.isEnd()) {
            auto simplex = circulator.simplex();
            if (!simplex) {
                assert(false);
                break;
            }
            circulator.insertIntersection(intersection);
            ++circulator;
        }
    }
}

template <unsigned N>
void SimplexDCIntersecter<N>::load(const std::array<Input*, 1 << N>& ts)
{
    // Skip this if all of the cells are empty / filled.
    if (std::all_of(ts.begin(), ts.end(),
        [](const Input* t)
    { return t->type != Interval::AMBIGUOUS; }))
    {
        return;
    }

    // We want to use the smallest cell adjoining this corner.
    const auto index = std::min_element(ts.begin(), ts.end(),
        [](const Input* a, const Input* b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();

    auto cornerSubspaceIndex = CornerIndex((1 << N) - 1 - index).neighbor().i;

    auto cornerSub = ts[index]->leaf->collapsedSub(cornerSubspaceIndex);

    constexpr auto tsIndices = [](){
        static_assert(N == 2 || N == 3);
        if constexpr (N == 2) {
            return std::array{ 0, 1, 2, 3 };
        }
        else {
            return std::array{ 0, 1, 2, 3, 4, 5, 6, 7 };
        }
    }();

    std::array<SimplexDCEdge<N>*, N * 2> edges;

    // Handle edges between edge vertices and our corner vertex.
    for (auto edgeAxisIdx = 0; edgeAxisIdx < N; ++edgeAxisIdx) {
        auto edgeAxis = Axis::toAxis(edgeAxisIdx);
        for (auto edgePosition = 0; edgePosition < 2; ++edgePosition) {
            auto lowerLevel = [&](int indexA, int indexB) {
                if (bool(indexA & edgeAxis) != bool(edgePosition)) {
                    // That index cell does not border the edge we're 
                    // looking at.
                    return false;
                }
                else if (bool(indexB & edgeAxis) != bool(edgePosition)) {
                    // The second cell is the only one of the two that 
                    // does border the edge we're looking at.
                    return true;
                }
                return ts[indexA]->leafLevel() < ts[indexB]->leafLevel();
            };
            auto edgeIndex = *std::min_element(
                tsIndices.begin(), tsIndices.end(), lowerLevel);
            if (ts[edgeIndex]->type == Interval::UNKNOWN) {
                continue;
            }
            auto edgeIsInternal = false;
            for (auto compareAxis : { Axis::Q(edgeAxis), Axis::R(edgeAxis) }) {
                if (N == 2 && compareAxis == Axis::Z) {
                    continue;
                }
                if (ts[edgeIndex] == ts[edgeIndex ^ compareAxis]) {
                    edgeIsInternal = true;
                    break;
                }
            }
            if (edgeIsInternal) {
                // The edge is internal to a merged cell or face, hence does
                // not exist, so nothing more to do for this edge.
                edges[edgeAxisIdx * 2 + edgePosition] = nullptr;
                continue;
            }
            auto lowerCorner = ~edgeIndex & ~edgeAxis & ((1 << N) - 1);
            const auto& edge = ts[edgeIndex]->leaf->edge(edgeAxis, lowerCorner);
            edges[edgeAxisIdx * 2 + edgePosition] = edge;
            auto edgeSubIndex = CornerIndex(lowerCorner).neighbor().i + 
                2 * ipow(3, Axis::toIndex(edgeAxis));
            auto edgeSub = ts[edgeIndex]->leaf->collapsedSub(edgeSubIndex);
            auto isSame = cornerSub == edgeSub;
            if (!isSame && cornerSub->inside == edgeSub->inside) {
                // This is the only condition we need to check in this case, as
                // a corner/edge pair is always the minimal-dimensional one for
                // that pair of subs.
                continue;
            }
            auto inside = cornerSub->inside ? cornerSub : edgeSub;
            auto outside = cornerSub->inside ? edgeSub : cornerSub;
            auto intersection = isSame ? &DCSimplex<N>::dupVertIntersection : 
                searchEdge(inside, outside, ts[edgeIndex]->leaf->tape);

            auto myCornerSub = 
                CornerIndex((1 << N) - 1 - edgeIndex).neighbor().i;
            auto insideIdx = cornerSub->inside ? myCornerSub : edgeSubIndex;
            auto outsideIdx = cornerSub->inside ? edgeSubIndex : myCornerSub;
            SimplexDCCirculator<N, false> circulator(
                ts[edgeIndex], insideIdx, outsideIdx, isSame);
            while (!circulator.isEnd()) {
                auto simplex = circulator.simplex();
                if (!simplex) {
                    assert(false);
                    break;
                }
                circulator.insertIntersection(intersection);
                ++circulator;
            }
        }
    }

    if constexpr (N == 3) {
        // Handle edges between face vertices and our corner vertex.
        for (auto faceAxisIdx = 0; faceAxisIdx < N; ++faceAxisIdx) {
            auto faceAxis = Axis::toAxis(faceAxisIdx);
            for (auto facePosition = 0; facePosition < 4; ++facePosition) {
                auto shift = faceAxisIdx + 1;
                auto lowerFaceCell =
                    (facePosition << shift) & 7 | (facePosition >> (N - shift));
                auto upperFaceCell = lowerFaceCell | faceAxis;
                if (ts[lowerFaceCell] == ts[upperFaceCell]) {
                    // The face is internal to the cell, hence does not exist,
                    // so don't use it.
                    continue;
                }
                auto useUpper = ts[upperFaceCell]->leafLevel() < 
                                ts[lowerFaceCell]->leafLevel();
                auto faceIndex = useUpper ? upperFaceCell : lowerFaceCell;
                if (ts[faceIndex]->type == Interval::UNKNOWN) {
                    continue;
                }
                auto isDuplicate = false;
                for (auto testAxis : { Axis::Q(faceAxis), Axis::R(faceAxis) }) {
                    if (faceIndex & testAxis &&
                        ts[faceIndex] == ts[faceIndex ^ testAxis]) {
                        // We already did this face from its other position.
                        isDuplicate = true;
                        break;
                    }
                }
                if (isDuplicate) {
                    continue;
                }
                auto faceSubIndex = (ipow(3, N) - 1) - 
                    ipow(3, faceAxisIdx) * (int(useUpper) + 1);
                auto faceSub = ts[faceIndex]->leaf->collapsedSub(faceSubIndex);
                auto isSame = cornerSub == faceSub;
                if (!isSame) {
                    if (cornerSub->inside == faceSub->inside) {
                        // Not an active edge.
                        continue;
                    }
                    // In our current (face/corner) case, the only condition required
                    // for this to be a minimal pair is that the face and corner do
                    // not share any axes for which the corner is at an extremum of
                    // the face (other than the face axis itself).

                    auto isMinimalPair = true;
                    for (auto testAxis = 0; testAxis < N; ++testAxis) {
                        if (testAxis == faceAxisIdx) {
                            continue;
                        }
                        auto vertVal = cornerSub->vert(testAxis);
                        if (vertVal != ts[faceIndex]->region.lower(testAxis) &&
                            vertVal != ts[faceIndex]->region.upper(testAxis)) {
                            continue;
                        }
                        if (vertVal == faceSub->vert(testAxis)) {
                            isMinimalPair = false;
                            break;
                        }
                    }
                    if (!isMinimalPair) {
                        continue;
                    }
                }
                auto inside = cornerSub->inside ? cornerSub : faceSub;
                auto outside = cornerSub->inside ? faceSub : cornerSub;
                auto intersection = isSame ? &DCSimplex<N>::dupVertIntersection :
                    searchEdge(inside, outside, ts[faceIndex]->leaf->tape);

                auto circulator = [&]() {
                    auto Q = Axis::Q(faceAxis);
                    auto R = Axis::R(faceAxis);
                    assert(ts[faceIndex] != ts[faceIndex ^ (Q | R)]);
                    auto edgeAxis = Q;
                    // Get the smallest cell bordering our edge.  We don't
                    // need to check faceIndex ^ faceAxis since we already
                    // compared it to faceIndex when selecting faceIndex,but
                    // there are two more cells (unless there are double cells)
                    // adjacent to our edge.
                    auto thirdCell = faceIndex ^ R;
                    auto edgeIndex = std::min(
                        { faceIndex, thirdCell, thirdCell ^ faceAxis },
                        [&](int a, int b)
                    { return ts[a]->leafLevel() < ts[b]->leafLevel(); });
                    if (ts[edgeIndex] == ts[edgeIndex ^ R]) {
                        edgeAxis = R;
                        auto thirdCell = faceIndex ^ Q;
                        edgeIndex = std::min(
                            { faceIndex, thirdCell, thirdCell ^ faceAxis },
                            [&](int a, int b)
                        { return ts[a]->leafLevel() < ts[b]->leafLevel(); });
                        assert(ts[edgeIndex] != ts[edgeIndex ^ Q]);
                    }
                    auto myCornerIdx = CornerIndex(7 - edgeIndex).neighbor().i;
                    if (edgeIndex != faceIndex) {
                        return SimplexDCCirculator<3, false>(
                            SimplexDCCirculator<3, true>::FarFaceTag{},
                            ts[edgeIndex], myCornerIdx,
                            edgeAxis, faceAxis, cornerSub->inside, isSame);
                    }
                    else {
                        auto insideIdx =
                            cornerSub->inside ? myCornerIdx : faceSubIndex;
                        auto outsideIdx =
                            cornerSub->inside ? faceSubIndex : myCornerIdx;
                        return SimplexDCCirculator<N, false>(
                            ts[edgeIndex], insideIdx, outsideIdx, isSame);
                    }
                }();

                while (!circulator.isEnd()) {
                    auto simplex = circulator.simplex();
                    if (!simplex) {
                        assert(false);
                        break;
                    }
                    circulator.insertIntersection(intersection);
                    ++circulator;
                }
            }
        }
    }

    // Now handle edges between cell vertices and our corner vertex.
    for (auto cell = 0; cell < ts.size(); ++cell) {
        if (ts[cell]->type == Interval::UNKNOWN) {
            continue;
        }
        auto isUpperDuplicate = false;
        for (auto testAxis : { Axis::X, Axis::Y, Axis::Z }) {
            if (!(testAxis & cell)) {
                // Even if it's a duplicate, it's the lower one, and we need
                // to handle one of them.
                continue;
            }
            else if (ts[cell] == ts[cell ^ testAxis]) {
                assert(cell != index);
                assert ((cell ^ testAxis) != index);
                isUpperDuplicate = true;
                break;
            }
        }
        if (isUpperDuplicate) {
            continue;
        }
        auto cellSub = ts[cell]->leaf->collapsedSub(ipow(3, N) - 1);
        auto isSame = cornerSub == cellSub;
        if (!isSame) {
            if (cornerSub->inside == cellSub->inside) {
                // Not an active edge.
                continue;
            }
            // In our current (cell/corner) case, the condition required
            // for this to be a minimal pair is that the face and corner do
            // not share any axes (other than possibly axes in which the corner
            // is not on the border of the cell, but rather an edge or face of
            // a merged cell).

            auto isMinimalPair = true;
            for (auto testAxis = 0; testAxis < N; ++testAxis) {
                auto vertValue = cornerSub->vert(testAxis);
                if (vertValue != ts[cell]->region.lower(testAxis) &&
                    vertValue != ts[cell]->region.upper(testAxis)) {
                    continue;
                }
                if (vertValue == cellSub->vert(testAxis)) {
                    isMinimalPair = false;
                    break;
                }
            }
            if (!isMinimalPair) {
                continue;
            }
        }
        auto inside = cornerSub->inside ? cornerSub : cellSub;
        auto outside = cornerSub->inside ? cellSub : cornerSub;
        auto intersection = isSame ? &DCSimplex<N>::dupVertIntersection :
            searchEdge(inside, outside, ts[cell]->leaf->tape);

        auto cellSubIdx = NeighborIndex::fromPosAndFloating(
            cell, ((1 << N) - 1) ^ (cell ^ index)).i + ipow(3, N);
        auto insideIdx = cornerSub->inside ? cornerSubspaceIndex : cellSubIdx;
        auto outsideIdx = cornerSub->inside ? cellSubIdx : cornerSubspaceIndex;
        SimplexDCCirculator<N, false> circulator(
            ts[index], insideIdx, outsideIdx, isSame);
        while (!circulator.isEnd()) {
            auto simplex = circulator.simplex();
            if (!simplex) {
                assert(false);
                break;
            }
            circulator.insertIntersection(intersection);
            ++circulator;
        }
    }
}

template<unsigned N>
SimplexDCIntersection<N>* SimplexDCIntersecter<N>::searchEdge(
    const SimplexLeafSubspace<N>* inside,
    const SimplexLeafSubspace<N>* outside,
    const std::shared_ptr<Tape>& tape)
{
    // This code is based on the simplex mesher, but without the assumption
    // that N == 3, and with the gradient calculated after finding the
    // intersection (as we need it to do the DC part of simplex DC).
    assert(tape.get() != nullptr);

    // There's an interesting question of precision + speed tradeoffs,
    // which mostly depend on how well evaluation scales in the
    // ArrayEvaluator.  for now, we'll use the same value as XTree.
    constexpr int SEARCH_COUNT = 4;
    constexpr int POINTS_PER_SEARCH = 16;
    static_assert(POINTS_PER_SEARCH <= ArrayEvaluator::N,
        "Overflowing ArrayEvaluator data array");

    // We are currently evaluating only one intersection at a time.  Later, we
    // will likely want to do more, and will want to loop over them; the loop
    // code is partially preserved from the DC code, but a lot will still need 
    // to be changed to make it functional.

    constexpr auto eval_count = 1;

    std::array<std::pair<Eigen::Vector3d, Eigen::Vector3d>, eval_count> targets;
    targets[0].first << inside->vert.transpose(), perp;
    targets[0].second << outside->vert.transpose(), perp;

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

        auto out = eval->values(POINTS_PER_SEARCH, *tape);

        // Skip one point, because the very first point is
        // already known to be inside the shape (but
        // sometimes, due to numerical issues, it registers
        // as outside!)
        for (unsigned j = 1; j < POINTS_PER_SEARCH; ++j)
        {
            // We're searching for the first point that's outside of the
            // surface.  There's a special case for the final point in the
            // search, working around  numerical issues where different
            // evaluators disagree with whether points are inside or outside.
            if (out[j] > 0 || j == POINTS_PER_SEARCH - 1 ||
                (out[j] == 0 && !eval->isInside(
                    ps.col(j).template cast<float>(), tape)))

            {
                targets[0] = { ps.col(j - 1), ps.col(j) };
                break;
            }
        }
    }

    for (unsigned i = 0; i < eval_count; ++i)
    {
        eval->set(targets[i].first.template cast<float>(), 2 * i);
        eval->set(targets[i].second.template cast<float>(), 2 * i + 1);
    }

    // Copy the results to a local array, to avoid invalidating
    // the results array when we call features() below.
    Eigen::Array<float, 4, ArrayEvaluator::N> ds;
    ds.leftCols(2 * eval_count) = eval->derivs(
        2 * eval_count, *tape);
    auto ambig = eval->getAmbiguous(2 * eval_count, *tape);


    auto intersection = my_pool.get();

    // Iterate over all inside-outside pairs, storing the number
    // of intersections before each inside node (in prev_size),
    // then checking the rank of the pair after each outside
    // node based on the accumulated intersections.
    for (unsigned i = 0; i < 2 * eval_count; ++i)
    {
        static_assert(eval_count == 1, 
            "Time to change the return value for multiple intersections");

        // This is the position associated with the intersection
        // being investigated.
        const auto& pos = (i & 1) ? targets[i / 2].second : 
                                    targets[i / 2].first;

        // If this position is unambiguous, then we can use the
        // derivatives value calculated and stored in ds.
        if (!std::isfinite(ds.col(i).w())) {
            continue;
        }
        auto pushDeriv = 
            [&pos, &ds, &intersection, i](auto deriv) {

            if (std::isfinite(deriv.x()) && std::isfinite(deriv.y()) &&
                std::isfinite(deriv.z())) {
                Eigen::Matrix<double, N, 1> derivProper =
                    deriv.template cast<double>().template head<N>();
            intersection->push(pos.template head<N>(), 
                               derivProper,
                               ds.col(i).w());
            auto norm = derivProper.norm();
            if (norm > 1e-20) {
                auto normalized = derivProper / norm;
                intersection->AtANormalized +=
                    normalized * normalized.transpose();
            }
          }
        };

        if (!ambig(i))
        {
            pushDeriv(ds.col(i));
        }
        // Otherwise, we need to use the feature-finding special
        // case to find all possible derivatives at this point.
        else
        {
            const auto fs = eval->features(
                pos.template cast<float>(), tape);

            for (auto& f : fs)
            {
                pushDeriv(f);
            }
        }
    }

    auto ptDouble = intersection->normalized_mass_point().template head<N>();
    auto pt = ptDouble.template cast<float>().eval();

    decltype(inside) closeSub = nullptr;
    auto totalDistanceSquared = (outside->vert - inside->vert).squaredNorm();
    if (totalDistanceSquared * merge_ratio_squared >
        (outside->vert - ptDouble.transpose()).squaredNorm()) {
        closeSub = outside;
    }
    else if (totalDistanceSquared * merge_ratio_squared >
        (inside->vert - ptDouble.transpose()).squaredNorm()) {
        closeSub = inside;
    }

    m.insert(closeSub, pt, intersection);

    intersection->orientationChecker.setEndpts(outside->vert, inside->vert);

    // Calculate the rank now, to avoid threading issues.
    intersection->rank = 0;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, N, N>> 
        es(intersection->AtANormalized);
    auto eigenvalues = es.eigenvalues().real();
    for (unsigned j = 0; j < N; ++j) {
        intersection->rank += (fabs(eigenvalues[j]) >= EIGENVALUE_CUTOFF);
    }
    return intersection;
}

}   // namespace libfive
