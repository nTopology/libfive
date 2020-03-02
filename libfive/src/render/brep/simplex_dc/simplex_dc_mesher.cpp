/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <numeric>
#include <boost/container/static_vector.hpp>

#include "libfive/render/brep/simplex_dc/simplex_dc_mesher.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_circulator.hpp"
#include "libfive/render/brep/per_thread_brep.hpp"

namespace libfive {

template<Axis::Axis A>
void SimplexDCMesher::load(const std::array<Input*, 2>& ts)
{
    // Skip this if all of the cells are empty / filled.
    if (std::all_of(ts.begin(), ts.end(),
        [](const Input* t)
    { return t->type != Interval::AMBIGUOUS; }))
    {
        return;
    }

    // We want to use the smallest cell adjoining this face to get the
    // face subspace.
    const auto index = std::min_element(ts.begin(), ts.end(),
        [](const Input* a, const Input* b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();

    auto faceSubspaceIndex = 26 - ipow(3, Axis::toIndex(A)) * (index + 1);

    unsigned trueFaceIdx;
    const auto& faceSub = *ts[index]->leaf->collapsedSub(faceSubspaceIndex,
                                                         &trueFaceIdx);
    if (trueFaceIdx != faceSubspaceIndex) {
        assert(trueFaceIdx < faceSubspaceIndex);
        // As in the case of the intersecter, we want to only handle minimal
        // pairs of dimensionality for these subs.
        return;
    }
    for (auto cell = 0; cell < 2; ++cell) {
        if (ts[cell]->type != Interval::AMBIGUOUS) {
            continue;
        }
        const auto& cellSub = *ts[cell]->leaf->collapsedSub(26);
        if (faceSub.inside == cellSub.inside) {
            continue;
        }
        if (trueFaceIdx != faceSubspaceIndex) {
            assert(trueFaceIdx < faceSubspaceIndex);
            // If the face sub collapses to a lower dimension, these
            // are not the minimal dimensions for this pair of subs;
            // we'll handle this pair (or already did, or are currently
            // doing so) on the minimal-dimension pair.  If it does not,
            // then it is the minimal pair, since the cell sub does not
            // collapse to the face sub.  (The cell sub is not allowed
            // to collapse to an edge or corner sub in the face sub if
            // the face sub does not collapse as well.)
            continue;
        }
        auto cellSubIdx = (cell == index) ? ipow(3, 3) - 1
            : ipow(3, 3) + faceSubspaceIndex;
        auto insideIdx = faceSub.inside ? faceSubspaceIndex : cellSubIdx;
        auto outsideIdx = faceSub.inside ? cellSubIdx : faceSubspaceIndex;
        SimplexDCCirculator<3, true> circulator(
            ts[index], insideIdx, outsideIdx, false);
        addPolygon(circulator);
    }
}

template<Axis::Axis A>
void SimplexDCMesher::load(const std::array<Input*, 4>& ts)
{
    // Skip this if all of the cells are empty / filled.
    if (std::all_of(ts.begin(), ts.end(),
        [](const Input* t)
    { return t->type != Interval::AMBIGUOUS; }))
    {
        return;
    }

    // We want to use the smallest cell adjoining this edge.
    const auto index = std::min_element(ts.begin(), ts.end(),
        [](const Input* a, const Input* b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();

    // Calculate the subspace for the edge.
    std::array<int, 3> reductions{ 0, 0, 0 };
    reductions[(Axis::toIndex(A) + 1) % 3] = 1 + (index & 1);
    reductions[(Axis::toIndex(A) + 2) % 3] = 1 + bool(index & 2);

    auto totalReductions = reductions[0] +
        3 * reductions[1] +
        9 * reductions[2];

    auto edgeSubspaceIndex = 26 - totalReductions;

    unsigned trueEdgeIdx;
    const auto& edgeSub = *ts[index]->leaf->collapsedSub(
        edgeSubspaceIndex, &trueEdgeIdx);
    if (trueEdgeIdx != edgeSubspaceIndex) {
        assert(trueEdgeIdx < edgeSubspaceIndex);
        // Again, we handle the minimal-dimensionality condition just as in
        // the intersecter.
        return;
    }

    auto& edge = ts[index]->leaf->edgeFromReduced(A, ts.size() - 1 - index);

    // Handle edges between face vertices and our edge vertex.
    for (auto faceAxisIdx = 0; faceAxisIdx < 2; ++faceAxisIdx) {
        auto faceAxis = (faceAxisIdx == 0) ? Axis::Q(A) : Axis::R(A);
        for (auto facePosition = 0; facePosition < 2; ++facePosition) {
            auto indexA = facePosition == 0 ? 0 : (2 >> faceAxisIdx);
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
                continue;
            }
            auto indexFromFace = bool(faceIndex & (1 << faceAxisIdx));
            auto faceSubIndex = 
                26 - ipow(3, Axis::toIndex(faceAxis)) * (indexFromFace + 1);
            const auto& faceSub =
                *ts[faceIndex]->leaf->collapsedSub(faceSubIndex);
            if (edgeSub.inside == faceSub.inside) {
                continue;
            }
            auto circulator = [&]() {
                if (faceIndex == index) {
                    auto insideIdx =
                        edgeSub.inside ? edgeSubspaceIndex : faceSubIndex;
                    auto outsideIdx =
                        edgeSub.inside ? faceSubIndex : edgeSubspaceIndex;
                    return SimplexDCCirculator<3, true>(
                        ts[index], insideIdx, outsideIdx, false);
                }
                else {
                    assert((faceIndex ^ index) != (1 << faceAxisIdx));
                    return SimplexDCCirculator<3, true>(
                        SimplexDCCirculator<3, true>::FarFaceTag{},
                        ts[index], edgeSubspaceIndex, A, faceAxis,
                        edgeSub.inside, false);
                }
            }();

            addPolygon(circulator);
        }
    }

    // Now handle edges between cell vertices and our edge vertex.
    for (auto cell = 0; cell < ts.size(); ++cell) {
        if (ts[cell]->type == Interval::UNKNOWN) {
            continue;
        }
        auto opposite = index ^ 3;
        if (cell != opposite && ts[cell] == ts[opposite]) {
            assert(cell != index);
            // We have a duplicate cell, and we're handling it when 
            // cell == opposite.  The circulator will cover all the 
            // simplices around this pair of sub vertices, so there's
            // nothing more to do.
            continue;
        }
        const auto& cellSub = *ts[cell]->leaf->collapsedSub(26);
        if (edgeSub.inside == cellSub.inside) {
            continue;
        }

        // Just as in the case of the intersecter, we have another test to
        // make sure that we have a minimal dimensionality pair.
        auto isMinimalPair = true;
        for (auto testAxis = 0; testAxis < 3; ++testAxis) {
            if (testAxis == Axis::toIndex(A)) {
                continue;
            }
            auto vertValue = edgeSub.vert(testAxis);
            if (vertValue != ts[cell]->region.lower(testAxis) &&
                vertValue != ts[cell]->region.upper(testAxis)) {
                continue;
            }
            if (vertValue == cellSub.vert(testAxis)) {
                isMinimalPair = false;
                break;
            }
        }
        if (!isMinimalPair) {
            continue;
        }
        auto cellSubIdx = 26;
        if (cell != index) {
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
                cellSubIdx = 26 - toSubtract;
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
        auto insideIdx = edgeSub.inside ? edgeSubspaceIndex : cellSubIdx;
        auto outsideIdx = edgeSub.inside ? cellSubIdx : edgeSubspaceIndex;
        SimplexDCCirculator<3, false> circulator(
            ts[index], insideIdx, outsideIdx, false);

        addPolygon(circulator);
    }
}

void SimplexDCMesher::load(const std::array<Input*, 8>& ts) {
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

    auto cornerSubspaceIndex = CornerIndex(7 - index).neighbor().i;

    const auto& cornerSub = *ts[index]->leaf->collapsedSub(cornerSubspaceIndex);
    constexpr std::array tsIndices{ 0, 1, 2, 3, 4, 5, 6, 7 };

    std::array<SimplexDCMinEdge<3>*, 6> edges;

    // Handle edges between edge vertices and our corner vertex.
    for (auto edgeAxisIdx = 0; edgeAxisIdx < 3; ++edgeAxisIdx) {
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
            const auto edgeIndex = *std::min_element(
                tsIndices.begin(), tsIndices.end(), lowerLevel);
            if (ts[edgeIndex]->type == Interval::UNKNOWN) {
                continue;
            }
            auto edgeIsInternal = false;
            for (auto compareAxis : { Axis::Q(edgeAxis), Axis::R(edgeAxis) }) {
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
            auto lowerCorner = ~edgeIndex & ~edgeAxis & 7;
            const auto& edgeVariant =
                ts[edgeIndex]->leaf->edge(edgeAxis, lowerCorner);
            assert(std::holds_alternative<SimplexDCMinEdge<3>*>(edgeVariant));
            const auto& edge = std::get<SimplexDCMinEdge<3>*>(edgeVariant);
            edges[edgeAxisIdx * 2 + edgePosition] = edge;
            auto edgeSubIndex = CornerIndex(lowerCorner).neighbor().i + 
                2 * ipow(3, Axis::toIndex(edgeAxis));
            const auto& edgeSub =
                *ts[edgeIndex]->leaf->collapsedSub(edgeSubIndex);
            if (cornerSub.inside == edgeSub.inside) {
                continue;
            }

            auto myCornerSub = CornerIndex(7 - edgeIndex).neighbor().i;
            auto insideIdx = cornerSub.inside ? myCornerSub : edgeSubIndex;
            auto outsideIdx = cornerSub.inside ? edgeSubIndex : myCornerSub;
            SimplexDCCirculator<3, true> circulator(
                ts[edgeIndex], insideIdx, outsideIdx, false);
            addPolygon(circulator);
        }
    }

    // Next, handle edges between face vertices and our corner vertex.
    for (auto faceAxisIdx = 0; faceAxisIdx < 3; ++faceAxisIdx) {
        auto faceAxis = Axis::toAxis(faceAxisIdx);
        for (auto facePosition = 0; facePosition < 4; ++facePosition) {
            auto shift = faceAxisIdx + 1;
            auto lowerFaceCell =
                (facePosition << shift) & 7 | (facePosition >> (3 - shift));
            auto upperFaceCell = lowerFaceCell | faceAxis;
            if (ts[upperFaceCell] == ts[lowerFaceCell]) {
                // The face is internal to a merged cell.
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


            auto faceSubIndex = 
                26 - ipow(3, Axis::toIndex(faceAxis)) * (int(useUpper) + 1);
            const auto& faceSub =
                *ts[faceIndex]->leaf->collapsedSub(faceSubIndex);
            if (cornerSub.inside == faceSub.inside) {
                continue;
            }

            // Again, we use the same condition for minimal dimensionality as
            // in the intersecter

            auto isMinimalPair = true;
            for (auto testAxis = 0; testAxis < 3; ++testAxis) {
                if (testAxis == faceAxisIdx) {
                    continue;
                }
                if (testAxis == faceAxisIdx) {
                    continue;
                }
                auto vertVal = cornerSub.vert(testAxis);
                if (vertVal != ts[faceIndex]->region.lower(testAxis) &&
                    vertVal != ts[faceIndex]->region.upper(testAxis)) {
                    continue;
                }
                if (vertVal == faceSub.vert(testAxis)) {
                    isMinimalPair = false;
                    break;
                }
            }
            if (!isMinimalPair) {
                continue;
            }
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
                    return SimplexDCCirculator<3, true>(
                        SimplexDCCirculator<3, true>::FarFaceTag{},
                        ts[edgeIndex], myCornerIdx,
                        edgeAxis, faceAxis, cornerSub.inside, false);
                }
                else {
                    auto insideIdx =
                        cornerSub.inside ? myCornerIdx : faceSubIndex;
                    auto outsideIdx =
                        cornerSub.inside ? faceSubIndex : myCornerIdx;
                    return SimplexDCCirculator<3, true>(
                        ts[edgeIndex], insideIdx, outsideIdx, false);
                }
            }();

            addPolygon(circulator);
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
                assert((cell ^ testAxis) != index);
                isUpperDuplicate = true;
                break;
            }
        }
        if (isUpperDuplicate) {
            continue;
        }

        const auto& cellSub = *ts[cell]->leaf->collapsedSub(26);
        if (cornerSub.inside == cellSub.inside) {
            continue;
        }
        auto isMinimalPair = true;
        for (auto testAxis = 0; testAxis < 3; ++testAxis) {
            auto vertValue = cornerSub.vert(testAxis);
            if (vertValue != ts[cell]->region.lower(testAxis) &&
                vertValue != ts[cell]->region.upper(testAxis)) {
                continue;
            }
            if (vertValue == cellSub.vert(testAxis)) {
                isMinimalPair = false;
                break;
            }
        }
        if (!isMinimalPair) {
            continue;
        }
        auto cellSubIdx = NeighborIndex::fromPosAndFloating(
            cell, 7 ^ (cell ^ index)).i
            + ipow(3, 3);
        auto insideIdx = cornerSub.inside ? cornerSubspaceIndex : cellSubIdx;
        auto outsideIdx = cornerSub.inside ? cellSubIdx : cornerSubspaceIndex;
        SimplexDCCirculator<3, true> circulator(
            ts[index], insideIdx, outsideIdx, false);
        addPolygon(circulator);
    }
}

void SimplexDCMesher::addPolygon(SimplexDCCirculator<3, true> circulator)
{
    // Get the intersection mass point for a neater, and easier to prevent 
    // self-intersection in, triangulation.
    auto intersection = circulator.circulatingIntersection();
    if (!intersection) {
        assert(false);
        return;
    }
    assert(DCSimplex<3>::isValid(intersection));
    assert(intersection->index != 0);


    auto getIndex = [](auto simplex) {
        assert(simplex);
        auto out = simplex->index.load();
        assert(out != 0);
        return out;
    };

    auto getIntersection = [](const auto& circulator) {
        auto out = (*iter)->intersection(dim0, dim1);
        assert(DCSimplex<3>::isValid(out));
        return out;
    };

    auto intersectionVert = 
        intersection->normalized_mass_point().template head<3>().eval();
    auto skippedCount = 0; // For debugging check.

    while (!circulator.isEnd()) {
        Eigen::Matrix<uint32_t, 3, 1> triangle;
        auto vert1 = circulator.simplex()->vert;
        auto index1 = getIndex(circulator.simplex());
        // There is a second intersection shared between the two simplices
        // we are currently connecting, and we may need it to handle the
        // possibility that the line between two simplex vertices may pass
        // outside both of them.
        auto otherIntersection = circulator.passingIntersection(true);
        assert(otherIntersection != nullptr);
        assert(otherIntersection != intersection);
        if (otherIntersection->index == intersection->index) {
            // This happens if both are at their shared corner
            // (up to float precision).  In that case, the triangle this would
            // produce and the one that would be produced for this simplex and
            // a circulator around otherIntersection are flipped versions of
            // each other, so we simply remove (do not create) both of them.
            continue;
        }
        ++circulator;
        assert(circulator.circulatingIntersection() == intersection);
        assert(otherIntersection == circulator.passingIntersection(false));
        auto vert2 = circulator.simplex()->vert;
        auto index2 = getIndex(circulator.simplex());
        auto goodForThis = intersection->orientationChecker.check(
            vert1, vert2, intersectionVert);
        auto goodForOther = otherIntersection->orientationChecker.check(
            vert2, vert1, 
            otherIntersection->normalized_mass_point().template head<3>());
        switch (goodForThis + 2 * goodForOther) {
        case 3:
            // The usual case, the orientations of the simplices and their 
            // vertices are the same throughout.  Equivalently, the line 
            // between the vertices of our two simplices goes through the 
            // triangle between those simplices.
            triangle << intersection->index, index1, index2;
            m.branes.push_back(triangle);
            break;
        case 2:
            // The line between the vertices of our simplices goes past the
            // intersection we are currently going around; this would cause
            // our triangle to be "flipped", having geometric orientation
            // opposite what the positions of the simplices would imply.  We
            // can consider the simplex vertices as forming arcs around the
            // common edge; this can only happen if the arc between these two
            // vertices is greater than 180 degrees, and therefore should
            // happen at most once, assuming no degenerate simplices.  However,
            // we currently still have degenerate and near-degenerate simplices,
            // and that can cause there to be more skipped triangles, so we
            // disable this assert for now.
            assert(++skippedCount == 1);
            // We don't make any triangles in this case; the "missing" triangle
            // will be handled in a corresponding case 1, with the same two
            // simplices in the opposite order and the same two intersections
            // in opposite roles.
            break;
        case 0:
            // In this case, the line between the vertices of our simplices
            // goes past both active edges of the triangle between them.
            // Thus, it functions as both case 2 and case 1, and will
            // correspond to another case 0.  We thus want to treat one as
            // case 1 and one as case 2; the triangles created will be the
            // same either way, so we can use the index to determine which to
            // do.  Either way, this does depend on an arc of greater than 180
            // degrees, so we can increment and check skippedCount.
            assert(intersection->index != otherIntersection->index);
            assert(++skippedCount == 1);
            if (intersection->index < otherIntersection->index) {
                break;
            }
            // Otherwise, fallthrough.
        case 1:
            // In this case, the line between simplex vertices goes past
            // only the other active edge.  It therefore needs to be split
            // into three triangles along the intersection of the other
            // active edge; one of those three is the negation of the triangle
            // ignored (not made) by the corresponding case 2, so they cancel
            // out and we make neither.  The other two are between both 
            // intersection vertices and each of the simplex vertices.
            triangle << index2, intersection->index, otherIntersection->index;
            m.branes.push_back(triangle);

            triangle << index1, otherIntersection->index, intersection->index;
            m.branes.push_back(triangle);
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
//  Template initialization
template void SimplexDCMesher::load<Axis::X>(
        const std::array<SimplexDCMesher::Input*, 4>&);
template void SimplexDCMesher::load<Axis::Y>(
        const std::array<SimplexDCMesher::Input*, 4>&);
template void SimplexDCMesher::load<Axis::Z>(
        const std::array<SimplexDCMesher::Input*, 4>&);
template void SimplexDCMesher::load<Axis::X>(
    const std::array<SimplexDCMesher::Input*, 2>&);
template void SimplexDCMesher::load<Axis::Y>(
    const std::array<SimplexDCMesher::Input*, 2>&);
template void SimplexDCMesher::load<Axis::Z>(
    const std::array<SimplexDCMesher::Input*, 2>&);

}   // namespace libfive
