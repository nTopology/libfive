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

    const auto& faceSub = *ts[index]->leaf->sub[faceSubspaceIndex].load();

    for (auto cell = 0; cell < 2; ++cell) {
        if (ts[cell]->type != Interval::AMBIGUOUS) {
            continue;
        }
        const auto& cellSub = *ts[cell]->leaf->sub[26].load();
        if (faceSub.inside == cellSub.inside) {
            continue;
        }
        boost::container::small_vector<const DCSimplex<3>*, 8> simplices;
        // We will add the simplices, going counterclockwise around our simplex
        // edge when looking from the +A direction.  We'll start at corner 0 of
        // the face.  Our edges will thus go in order Q(increasing, at low R),
        // R(increasing, at high Q), Q(decreasing, at high R), 
        // R(decreasing, at low Q).
        for (auto edgeDirection = 1; edgeDirection >= 0; --edgeDirection) {
            for (auto axisIdx = 0; axisIdx < 2; ++axisIdx) {
                auto edgeAxis = axisIdx == 0 ? Axis::Q(A) : Axis::R(A);
                auto edgePos = (axisIdx + edgeDirection + 1) % 2;
                auto& edge = ts[index]->leaf->edgeFromFaceAndIndex(
                    edgeAxis, A, edgePos == 1, index == 0);
                auto pushSimplices = [&](SimplexDCMinEdge<3>* edge) {
                    for (auto cornerIdx = 0; cornerIdx < 2; ++cornerIdx) {
                        auto cornerPos = (cornerIdx + edgeDirection + 1) % 2;
                        auto corner = axisIdx == 0 ? 2 * edgePos + cornerPos
                            : edgePos + 2 * cornerPos;
                        auto& simplex = edge->simplexWithFaceReducedCell(
                            edgeAxis, A, cell == 1, corner);
                        assert(simplex.index.load() != 0);
                        simplices.push_back(&simplex);
                    }
                };
                if (std::holds_alternative<SimplexDCMinEdge<3>*>(edge)) {
                    pushSimplices(std::get<SimplexDCMinEdge<3>*>(edge));
                }
                else {
                    auto& edgeVec =
                        std::get<SimplexDCMinEdge<3>::EdgeVec>(edge);
                    edgeVec.sort();
                    if (edgeDirection == 1) {
                        for (auto& minEdge : edgeVec) {
                            pushSimplices(minEdge);
                        }
                    }
                    else {
                        for (auto iter = edgeVec.rbegin();
                            iter != edgeVec.rend();
                            ++iter) {
                            pushSimplices(*iter);
                        }
                    }
                }
            }
        }

        if (faceSub.inside == (cell == 0)) {
            // We should be looking from -A.
            std::reverse(simplices.begin(), simplices.end());
        }
        addPolygon(simplices, 2, 3);
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

    const auto& edgeSub = *ts[index]->leaf->sub[edgeSubspaceIndex].load();

    auto& edge = ts[index]->leaf->edgeFromReduced(A, ts.size() - 1 - index);

    // Handle edges between face vertices and our edge vertex.
    for (auto faceAxisIdx = 0; faceAxisIdx < 2; ++faceAxisIdx) {
        auto faceAxis = (faceAxisIdx == 0) ? Axis::Q(A) : Axis::R(A);
        for (auto facePosition = 0; facePosition < 2; ++facePosition) {
            auto indexA = facePosition == 0 ? 0 : 3;
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
                *ts[faceIndex]->leaf->sub[faceSubIndex].load();
            if (edgeSub.inside == faceSub.inside) {
                continue;
            }
            std::array<const DCSimplex<3>*, 4> simplices;

            // Again, we want to go counterclockwise, but in this case it is
            // easier if we assume that we're looking from the positive
            // direction (in this case, in the axis that is neither A nor
            // faceAxis, as that is the axis by which our edge and face differ)
            // for faceAxisIdx == 0, and from the negative direction for
            // faceAxisIdx == 1.  That way, if we start with the simplex where
            // both the corner and cell have their lower value (for their
            // respective axes), we will start by moving along the edge 
            // (switching to low cell/high corner next).
            
            for (auto i = 0; i < 4; ++i) {
                auto cell = i >> 1;
                assert(cell == 0 || cell == 1);
                auto corner = (i & 1) ^ cell;
                assert(corner == 0 || corner == 1);
                auto cellFrom4 = (faceAxisIdx ? facePosition + 2 * cell
                                              : 2 * facePosition + cell);
                auto& simplex =
                    std::get<SimplexDCMinEdge<3>*>(edge)
                    ->simplexWithEdgeReducedCell(
                        A, faceAxis, cellFrom4, corner);
                simplices[i] = &simplex;
            }

            if (edgeSub.inside == (facePosition == faceAxisIdx)) {
                std::reverse(simplices.begin(), simplices.end());
            }
            addPolygon(simplices, 1, 2);
        }
    }

    // Now handle edges between cell vertices and our edge vertex.
    for (auto cell = 0; cell < ts.size(); ++cell) {
        if (ts[cell]->type == Interval::UNKNOWN) {
            continue;
        }
        const auto& cellSub = *ts[cell]->leaf->sub[26].load();
        if (edgeSub.inside == cellSub.inside) {
            continue;
        }
        // In this case, we will assume counterclockwise looking from the cell
        // toward the edge; thus, we start with low corner/counterclockwise face
        // (as measured looking from the +A direction), and will begin moving
        // along the edge (switching corners).  To handle duplicate cells, we 
        // will use only the counterclockwise of the two (an edge cannot have
        // quadrupled cells), and instead of its clockwise face will be the 
        // clockwise face of its neighbor.

        std::array<const DCSimplex<3>*, 4> simplices;

        // Based on the cell, we can determine if the axis of the 
        // counterclockwise face is Q or R.
        auto QIsCounterClockwise = (cell == 0 || cell == 3);

        auto cClockwiseCell = cell ^ (QIsCounterClockwise ? 1 : 2);
        if (ts[cell] == ts[cClockwiseCell]) {
            continue;
        }

        auto clockwiseCell = 3 ^ cClockwiseCell;
        auto isDuplicate = ts[cell] == ts[clockwiseCell];

        for (auto i = 0; i < 4; ++i) {
            auto faceIsCounterClockwise = (i < 2);
            auto cellToUse = cell;
            auto faceAxis = (faceIsCounterClockwise == QIsCounterClockwise)
                ? Axis::Q(A) : Axis::R(A);
            if (isDuplicate && !faceIsCounterClockwise) {
                cellToUse = clockwiseCell;
                faceAxis = QIsCounterClockwise ? Axis::Q(A) : Axis::R(A);
            }
            auto corner = bool(i & 1) == faceIsCounterClockwise;
            auto& simplex = std::get<SimplexDCMinEdge<3>*>(edge)
                ->simplexWithEdgeReducedCell(A, faceAxis, cellToUse, corner);
            assert(simplex.index.load() != 0);
            simplices[i] = &simplex;
        }

        if (cellSub.inside) {
            std::reverse(simplices.begin(), simplices.end());
        }
        addPolygon(simplices, 1, 3);
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

    const auto& cornerSub = *ts[index]->leaf->sub[cornerSubspaceIndex].load();
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
                *ts[edgeIndex]->leaf->sub[edgeSubIndex].load();
            if (cornerSub.inside == edgeSub.inside) {
                continue;
            }

            // Again, we will go counterclockwise when looking from positive 
            // edgeAxis, and reverse afterward if that would be looking from
            // the inside out; we start at the low cell, clockwise face.

            boost::container::static_vector<const DCSimplex<3>*, 8> simplices;
            for (auto cellNum = 0; cellNum < 4; ++cellNum) {
                // Our index with respect to the edge is not cellNum,
                // since going counterclockwise means we need the order 
                // 0,1,3,2.
                auto cellIdx = cellNum ^ (cellNum >> 1);
                auto QIsCounterClockwise = (cellIdx == 0 || cellIdx == 3);
                auto fullCellIdx = cellIdx << 1;
                if (edgePosition == 1) {
                    fullCellIdx |= 1;
                }
                fullCellIdx <<= edgeAxisIdx;
                fullCellIdx |= (fullCellIdx >> 3);
                fullCellIdx &= 7;
                if (ts[fullCellIdx]->type == Interval::UNKNOWN) {
                    // This would imply that we have two points on the
                    // boundary of our original region (the edge vertex
                    // and corner vertex), of which one is inside and
                    // one is outside.
                    assert(false);
                    continue;
                }
                for (auto CCFace : { 0, 1 }) {
                    auto faceAxis = (QIsCounterClockwise == bool(CCFace))
                        ? Axis::Q(edgeAxis) : Axis::R(edgeAxis);
                    if (ts[fullCellIdx] == ts[fullCellIdx ^ faceAxis]) {
                        // There is a merged cell, so this face does not exist.
                        continue;
                    }
                    auto& simplex = edge->simplexWithEdgeReducedCell(
                        edgeAxis, faceAxis, cellIdx, !bool(edgePosition));
                    assert(simplex.index.load() != 0);
                    simplices.push_back(&simplex);
                }
            }

            if (edgeSub.inside == bool(edgePosition)) {
                std::reverse(simplices.begin(), simplices.end());
            }
            addPolygon(simplices, 0, 1);
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

            auto faceSubIndex = 
                26 - ipow(3, Axis::toIndex(faceAxis)) * (int(useUpper) + 1);
            const auto& faceSub =
                *ts[faceIndex]->leaf->sub[faceSubIndex].load();
            if (cornerSub.inside == faceSub.inside) {
                continue;
            }

            // This time, we will go counterclockwise when looking from the
            // face vertex toward the corner vertex.  We will start at the
            // low cell and clockwise (when looking from positive faceAxis)
            // edge, so the next will be low cell/counterclockwise edge.

            std::array<const DCSimplex<3>*, 4> simplices;

            // It is possible that our corner vertex is on an edge of the face.
            // In such case, the face must have a duplicate cell on each side,
            // and thus be duplicate itself.  (One may be quadrupled, but that
            // has no effect here, and if both were quadrupled this would not
            // be an active corner).  In the case of such a duplicate face, we 
            // want to use the face position that takes the lower value among
            // the duplicate axis (the other will not produce a polygon), and
            // its edges will both have that axis.  (There is no clockwise/
            // counterclockwise in such case, but there is upper/lower value
            // for the duplicate axis).

            bool isDuplicate = false;
            Axis::Axis duplicateAxis; // Undefined if !isDuplicate.

            // We only need to check the side with faceIndex; that has a level
            // lower than, or equal to, that of the other side, so it is smaller
            // or equal and therefore if only one side is merged it will be
            // the other one.
            if (ts[faceIndex] == ts[faceIndex ^ Axis::Q(faceAxis)]) {
                assert(ts[faceIndex] != ts[faceIndex ^ Axis::R(faceAxis)]);
                assert(ts[faceIndex ^ faceAxis] == 
                    ts[faceIndex ^ faceAxis ^ Axis::Q(faceAxis)]);
                isDuplicate = true;
                duplicateAxis = Axis::Q(faceAxis);
            }
            else if (ts[faceIndex] == ts[faceIndex ^ Axis::R(faceAxis)]) {
                assert(ts[faceIndex ^ faceAxis] ==
                    ts[faceIndex ^ faceAxis ^ Axis::R(faceAxis)]);
                isDuplicate = true;
                duplicateAxis = Axis::R(faceAxis);
            }
            if (isDuplicate && faceIndex & duplicateAxis) {
                // Higher of the two duplicates.
                continue;
            }

            auto QIsCC = (facePosition == 1 || facePosition == 2);

            for (auto cell : { 0, 1 }) {
                for (auto edgeIdx : { 0, 1 }) {
                    auto CCEdge = edgeIdx != cell;
                    auto edgeToUse = edgeIdx;
                    auto edgeAxisIdx = (CCEdge == QIsCC) ? 0 : 1;
                    auto edgeAxis = (edgeAxisIdx == 0) ? Axis::Q(faceAxis)
                                                       : Axis::R(faceAxis);
                    auto facePositionToUse = facePosition;
                    if (isDuplicate && edgeAxis != duplicateAxis) {
                        edgeAxis = duplicateAxis;
                        edgeAxisIdx = 1 - edgeAxisIdx;
                        facePositionToUse ^= (1 << edgeAxisIdx);
                    }
                    bool edgePosition(facePositionToUse & (1 << edgeAxisIdx));
                    auto& edge = 
                        edges[Axis::toIndex(edgeAxis) * 2 + edgePosition];
                    assert(edge != nullptr);

                    auto& simplex = edge->simplexWithFaceReducedCell(
                        edgeAxis, faceAxis, cell, 3 - facePositionToUse);
                    assert(simplex.intersectionCount() != 0);
                    simplices[cell * 2 + edgeIdx] = &simplex;
                }
            }

            if (faceSub.inside) {
                std::reverse(simplices.begin(), simplices.end());
            }
            addPolygon(simplices, 0, 2);
        }
    }

    // Now handle edges between cell vertices and our corner vertex.
    for (auto cell = 0; cell < ts.size(); ++cell) {
        if (ts[cell]->type == Interval::UNKNOWN) {
            continue;
        }
        const auto& cellSub = *ts[cell]->leaf->sub[26].load();
        if (cornerSub.inside == cellSub.inside) {
            continue;
        }
        // We'll go counterclockwise when looking from the cell vertex toward
        // the corner vertex when bitcount(cell) is odd, and from the corner
        // vertex toward the cell vertex when it is even.  That way, we'll
        // always be doing elements with the same edge and face axes; we'll
        // start with edge X face Y and then next will be edge X face Z.

        // It is possible that we have a duplicate, or even quadrupled, cell.
        // If it is duplicate, we want to use only the cell with the lower 
        // index; in such case, it will have 4, 6, or 8 simplices, depending
        // on which axes have duplicate faces to go with the duplicate cell.

        // If it is quadrupled, there will be 6 or 8 simplices, and we will 
        // use only the cell with the lowest index of the four.

        auto duplicateDirections = 0;
        for (auto axisIdx = 0; axisIdx < 3; ++axisIdx) {
            auto axis = Axis::toAxis(axisIdx);
            if (ts[cell] == ts[cell ^ axis]) {
                duplicateDirections |= axis;
            }
        }
        if (duplicateDirections & cell) {
            // There is a lower cell index for the same cell pointer.
            continue;
        }
        switch (duplicateDirections) {
        case 0:
        {
            std::array<const DCSimplex<3>*, 6> simplices;
            for (auto edgeAxisIdx = 0; edgeAxisIdx < 3; ++edgeAxisIdx) {
                auto edgeAxis = Axis::toAxis(edgeAxisIdx);
                auto edgePosition = bool(cell & edgeAxis);
                auto edge = edges[edgeAxisIdx * 2 + edgePosition];
                // This edge must exist, as if it were internal to a face or
                // cell, every cell bordering it would be duplicate if not
                // quadrupled.
                assert(edge != nullptr);                
                for (auto faceAxisIdx : { 0, 1 }) {
                    auto faceAxis = (faceAxisIdx == 0) ? Axis::Q(edgeAxis)
                        : Axis::R(edgeAxis);
                    auto& simplex = edge->simplex(edgeAxis, faceAxis,
                        cell, !edgePosition);
                    assert(simplex.index != 0);
                    simplices[2 * edgeAxisIdx + faceAxisIdx] = &simplex;
                }
            }
            if (cellSub.inside == bool(bitcount(cell) & 1)) {
                std::reverse(simplices.begin(), simplices.end());
            }
            addPolygon(simplices, 0, 3);
            break;
        }
        case 1:
        case 2:
        case 4:
        {
            auto dupAxis = static_cast<Axis::Axis>(duplicateDirections);
            // We have a duplicate cell.  We will go counterclockwise looking
            // from the cell vertex toward the corner vertex if the other two
            // axes are both positive or both negative and vice versa if they
            // are one positive and one negative; we will start going from the
            // edge at positive duplicateAxis.  Thus, our edges will go
            // +duplicateAxis, Q, -duplicateAxis, R, and back to +duplicateAxis;
            // we will omit any that use an edge interior to a duplicate face.
            // Our faces will thus go +R, -R, -Q, +Q.
            boost::container::static_vector<const DCSimplex<3>*, 8> simplices;
            for (auto axisIdx = 0; axisIdx < 2; ++axisIdx) {
                auto faceAxis = 
                    (axisIdx == 0) ? Axis::R(dupAxis) : Axis::Q(dupAxis);
                for (auto facePosIdx = 0; facePosIdx < 2; ++facePosIdx) {
                    auto facePos = facePosIdx ^ axisIdx ^ 1;
                    for (auto edgeIdx = 0; edgeIdx < 2; ++edgeIdx) {
                        Axis::Axis edgeAxis;
                        bool edgePos;
                        auto cellToUse = cell;
                        if (edgeIdx == facePosIdx) {
                            // The edge is along duplicateAxis, and must exist.
                            edgeAxis = dupAxis;
                            edgePos = facePos;
                        }
                        else {
                            // The edge is either Q or R, if it exists at all.
                            edgeAxis = (axisIdx == 0) ? Axis::Q(dupAxis) 
                                                      : Axis::R(dupAxis);
                            edgePos = edgeAxis & cell;
                            if (facePos) {
                                cellToUse |= dupAxis;
                            }
                            else {
                                cellToUse &= ~dupAxis;
                            }
                        }
                        auto edgeAxisIdx = Axis::toIndex(edgeAxis);
                        auto edge = edges[edgeAxisIdx * 2 + edgePos];
                        if (edge == nullptr) {
                            // The edge does not exist.
                            continue;
                        }
                        auto& simplex = edge->simplex(edgeAxis, faceAxis,
                            cellToUse, !edgePos);
                        assert(simplex.index != 0);
                        simplices.push_back(&simplex);
                    }
                }
            }
            if (cornerSub.inside == bool(bitcount(cell & ~dupAxis) & 1)) {
                std::reverse(simplices.begin(), simplices.end());
            }
            assert(simplices.size() >= 4);
            assert(simplices.size() % 2 == 0);
            addPolygon(simplices, 0, 3);
            break;
        }
        case 3:
        case 5:
        case 6:
        {
            auto faceAxis = static_cast<Axis::Axis>(7 ^ duplicateDirections);
            // We have a quadrupled cell.  We will go counterclockwise looking
            // from the +faceAxis direction, starting and ending at the edge
            // at +Q(faceAxis) (if it exists).  Thus, our edge order will be
            // +Q, +R, -Q, -R, and back to -Q, and our cell order (reduced by
            // faceAxis) will be 3, 2, 0, 1.  It is possible for one of our
            // edges to be nonexistent (so we skip those two simplices), but
            // if two were nonexistent we'd have only three cells around this
            // corner, in which case it would be in the middle of a minimal
            // edge and not handled by the dual walker.
            boost::container::static_vector<const DCSimplex<3>*, 8> simplices;
            for (auto reducedCell : { 3, 2, 0, 1 }) {
                for (auto edgeAxisIdx : { 0, 1 }) {
                    auto QIsCC = reducedCell == 1 || reducedCell == 2;
                    auto useQ = bool(edgeAxisIdx) == QIsCC;
                    auto edgeAxis = useQ ? Axis::Q(faceAxis) : Axis::R(faceAxis);
                    auto edgeAxisIdx = Axis::toIndex(edgeAxis);
                    auto faceAxisIdx = Axis::toIndex(faceAxis);
                    bool edgePos(reducedCell & (useQ ? 1 : 2));
                    auto edge = edges[edgeAxisIdx * 2 + edgePos];
                    if (edge == nullptr) {
                        // The edge does not exist.
                        continue;
                    }
                    auto reducedCellAsAbsolute = 
                        reducedCell << (faceAxisIdx + 1);
                    reducedCellAsAbsolute |= (reducedCellAsAbsolute >> 3);
                    reducedCellAsAbsolute &= 7;
                    auto cellFromEdge =
                        (cell & faceAxis) | reducedCellAsAbsolute;
                    auto& simplex = edge->simplex(
                        edgeAxis, faceAxis, cellFromEdge, !edgePos);
                    assert(simplex.intersectionCount() != 0);
                    assert(simplex.index != 0);
                    simplices.push_back(&simplex);
                }
            }
            if (cellSub.inside == bool(cell & faceAxis)) {
                std::reverse(simplices.begin(), simplices.end());
            }
            assert(simplices.size() == 6 || simplices.size() == 8);
            addPolygon(simplices, 0, 3);
            break;
        }
        case 7:
        default:
            assert(false);
        }
    }
}

template<class Container>
void SimplexDCMesher::addPolygon(const Container& simplices, int dim0, int dim1)
{
    assert(simplices.size() >= 3);
    // Get the intersection mass point for a neater, and easier to prevent 
    // self-intersection in, triangulation.
    auto intersection = simplices[0]->intersection(dim0, dim1);
    assert(intersection != nullptr);
    assert(intersection->index != 0);
    for (auto& simplex : simplices) {
        assert(simplex->intersection(dim0, dim1) == intersection);
    }

    auto getIndex = [](auto iter) {
        auto out = (*iter)->index.load();
        assert(out != 0);
        return out;
    };
    auto getVert = [](auto iter) {
        return (*iter)->vert;
    };
    auto getIntersection = [](auto iter, int dim0, int dim1) {
        return (*iter)->intersection(dim0, dim1);
    };
    auto intersectionVert = 
        intersection->normalized_mass_point().template head<3>().eval();
    auto skippedCount = 0; // For debugging check.
    for (auto iter1 = simplices.begin(); iter1 != simplices.end(); ++iter1) {
        Eigen::Matrix<uint32_t, 3, 1> triangle;
        auto iter2 = iter1 + 1;
        if (iter2 == simplices.end()) {
            iter2 = simplices.begin();
        }
        const SimplexDCIntersection<3>* otherSharedIntersection = nullptr;

        for (auto i = 0; i < 4; ++i) {
            for (auto j = i + 1; j < 4; ++j) {
                assert(dim0 < dim1); // Otherwise, the following condition 
                                     // needs to handle both directions.
                if (i == dim0 && j == dim1) {
                    continue;
                }
                auto intersection1 = getIntersection(iter1, i, j);
                if (intersection1 != nullptr &&
                    intersection1 == getIntersection(iter2, i, j)) {
                    assert(otherSharedIntersection == nullptr);
                    otherSharedIntersection = intersection1;
                }
            }
        }
        assert(otherSharedIntersection != nullptr);

        auto goodForThis = intersection->orientationChecker.check(
            getVert(iter1), getVert(iter2), intersectionVert);
        auto goodForOther = otherSharedIntersection->orientationChecker.check(
            getVert(iter2), getVert(iter1), otherSharedIntersection->
            normalized_mass_point().template head<3>());

        switch (goodForThis + 2 * goodForOther) {
        case 3:
            // The usual case, the orientations of the simplices and their 
            // vertices are the same throughout.  Equivalently, the line 
            // between the vertices of our two simplices goes through the 
            // triangle between those simplices.
            triangle << intersection->index,
                getIndex(iter1),
                getIndex(iter2);

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
            // happen at most once.
            assert(++skippedCount == 1);
            // We don't make any triangles in this case; the "missing" triangle
            // will be handled in a corresponding case 1, with the same two
            // simplices in the opposite order and the same two intersections
            // in opposite roles.
            break;
        case 0:
            // In this case, the line between the vertices of our simplices
            // goes past both active edges of the triangle between them.
            // Thus, it functions as both case 2 and case 1 (corresponding
            // to itself), and will correspond to another case 0.  We thus
            // want to treat one as case 1 and one as case 2; the triangles
            // created will be the same either way, so we can use the index
            // to determine which to do.  Either way, this does depend on
            // an arc of greater than 180 degrees, so we can increment
            // and check skippedCount.
            assert(intersection->index != otherSharedIntersection->index);
            assert(++skippedCount == 1);
            if (intersection->index < otherSharedIntersection->index) {
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
            triangle << getIndex(iter2),
                intersection->index,
                otherSharedIntersection->index;

            m.branes.push_back(triangle);

            triangle << getIndex(iter1),
                otherSharedIntersection->index, 
                intersection->index;

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