/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "libfive/render/brep/simplex_dc/simplex_dc_circulator.hpp"
#include "libfive/render/brep/indexes.hpp"

namespace libfive {
template class SimplexDCCirculator<2, false>;

DCSimplex<2>* SimplexDCCirculator<2, true>::forwardSimplex(
    const SimplexDCTree<2>* cell, unsigned insideIdx, unsigned outsideIdx)
{
    auto getDim = [](unsigned index)->unsigned {
        if (index >= 8) {
            return 3;
        }
        else {
            return NeighborIndex(index).dimension();
        }
    };
    auto getSub = [cell](unsigned index)->SimplexLeafSubspace<2>* {
        if (index <= 8) {
            return cell->leaf->collapsedSub(index);
        }
        else {
            auto neighbor = cell->neighbor(index - 8);
            if (neighbor) {
                return neighbor->leaf->collapsedSub(8);
            }
            else {
                assert(false);
                return nullptr;
            }
        }
    };
    assert(getDim(insideIdx) != getDim(outsideIdx));

    auto insideSub = getSub(insideIdx);
    auto outsideSub = getSub(outsideIdx);
    if (insideSub == outsideSub) {
        assert(false);
        return nullptr;
    }
    if (insideSub->inside == outsideSub->inside) {
        assert(false);
        return nullptr;
    }
    // We can't check that insideSub is actually inside, and outsideSub 
    // actually outside, since this is also called with inside and outside
    // swapped to get the backward simplex.

    unsigned thirdIdx;
    auto insideDim = getDim(insideIdx);
    auto outsideDim = getDim(outsideIdx);

    switch (insideDim + 3 * outsideDim) {
    case 1:
    case 3:
    {
        // Corner/edge.
        NeighborIndex inside(insideIdx);
        NeighborIndex outside(outsideIdx);
        auto edgeBit = inside.floating() | outside.floating();
        auto otherBit = 3 ^ edgeBit;
        bool isHighEdge(inside.pos() & otherBit);
        assert(isHighEdge == bool(outside.pos() & otherBit));
        bool isHighCorner(edgeBit & (inside.pos() | outside.pos()));
        auto cornerIsCC = (isHighCorner == isHighEdge == (edgeBit == 2));
        auto insideIsCC = (cornerIsCC == (insideDim == 0));
        if (!insideIsCC) {
            // The cell is the one we're centered on right now.
            thirdIdx = 8;
        }
        else {
            // The cell is the one on the other side of the edge.
            auto neighborIdx = insideDim == 1 ? insideIdx : outsideIdx;
            auto neighbor = cell->neighbor(neighborIdx);
            if (!neighbor) {
                assert(false);
                return nullptr;
            }
            assert(!neighbor->isBranch());
            thirdIdx = 8 + neighborIdx;
        }
        break;
    }
    case 2:
    case 6:
    {
        // Corner/cell.
        auto cornerIdx = (insideDim == 0 ? insideIdx : outsideIdx);
        auto cellIdx = (insideDim == 0 ? outsideIdx : insideIdx);
        assert(cellIdx >= 8);
        auto cellDirectionDim = (cellIdx == 8) ? 2 : getDim(cellIdx - 9);
        if (cellDirectionDim == 0) {
            assert(cellIdx == cornerIdx + 9);
        }
        auto xIsCC = (cornerIdx == 0 || cornerIdx == 4) ==
            (cellDirectionDim % 2 == 0);
        auto ccIsForward = insideDim != 0;
        auto axisToUse = (xIsCC == ccIsForward) ? Axis::X : Axis::Y;
        auto otherAxis = (xIsCC == ccIsForward) ? Axis::Y : Axis::X;
        NeighborIndex corner(cornerIdx);
        if (cellDirectionDim == 1) {
            NeighborIndex cellNeighbor(cellIdx - 9);
            // Check that our cell does in fact adjoin our corner.
            assert(!((corner.pos() ^ cellNeighbor.pos()) 
                    & cellNeighbor.fixed()));
        }
        // Now that we have our axis to use and corner, we need to make
        // sure that we're going to go from a cell that is actually adjacent
        // to the edge we're about to use.
        auto neighborFloats = 
            (cellIdx == 8) ? 3 : NeighborIndex(cellIdx - 9).floating();
        auto nearNeighbor = cell->neighbor(
            NeighborIndex::fromPosAndFloating(
                corner.pos(), neighborFloats | otherAxis));
        auto farNeighbor = cell->neighbor(
            NeighborIndex::fromPosAndFloating(
                corner.pos(), neighborFloats & axisToUse));
        if (nearNeighbor == farNeighbor) {
            // Our corner is on a side of our cell, and the axis we wanted
            // to use is not available.  So we use the non-doubled cell 
            // opposite the one we originally marked as our cell, as that is
            // the one that is adjacent to our edge and has our corner as
            // a corner.
            assert(!(neighborFloats & axisToUse)); // NeighborFloats is either 
            // 0 if our cell index was connected by the corner, or 3 ^ axisToUse
            // if it was connected by an edge in the other direction.
            if (neighborFloats) {
                cell = cell->neighbor(
                    NeighborIndex::fromPosAndFloating(
                        corner.pos(), axisToUse));
                auto newCornerPos = axisToUse ^ corner.pos();
                corner = NeighborIndex::fromPosAndFloating(
                    newCornerPos, 0);
                cornerIdx = corner.i;
                while (cell->isBranch()) {
                    cell = cell->child(newCornerPos);
                }
            }
            cellIdx = 9 + NeighborIndex::fromPosAndFloating(
                3 ^ corner.pos(), otherAxis).i;
            axisToUse = (axisToUse == Axis::X ? Axis::Y : Axis::X);
        }
        else {
            auto movedBase = 3 ^ neighborFloats;
            auto nearCornerPos = corner.pos() ^ (movedBase & ~axisToUse);
            auto farCornerPos = corner.pos() ^ (movedBase | axisToUse);
            while (nearNeighbor->isBranch()) {
                nearNeighbor = nearNeighbor->child(nearCornerPos);
            }
            while (farNeighbor->isBranch()) {
                farNeighbor = farNeighbor->child(farCornerPos);
            }
            auto useNear = 
                nearNeighbor->leafLevel() <= farNeighbor->leafLevel();
            auto newCornerPos = useNear ? nearCornerPos : farCornerPos;
            if (newCornerPos != corner.pos()) {
                cell = useNear ? nearNeighbor : farNeighbor;
                bool cellIsNear(neighborFloats & otherAxis);
                if (useNear == cellIsNear) {
                    cellIdx = 8;
                }
                else {
                    cellIdx = 9 + NeighborIndex::fromPosAndFloating(
                        newCornerPos, axisToUse).i;
                }
                corner = NeighborIndex::fromPosAndFloating(
                    newCornerPos, 0);
                cornerIdx = corner.i;
            }
        }

        // At this point, we should have a cell next to the edge we want,
        // the correct axis associated with it, and the proper cell and
        // corner indices.
        auto edgeNeighbor = NeighborIndex::fromPosAndFloating(
            corner.pos(), axisToUse);
        insideIdx = (insideDim == 0 ? cornerIdx : cellIdx);
        outsideIdx = (insideDim == 0 ? cellIdx : cornerIdx);
        thirdIdx = edgeNeighbor.i;

        // Our inside and outside subs should not have changed, as all
        // index changes should account for the change (if any) to cell.
        // The dimensionalities should remain the same as well.
        assert(getSub(insideIdx) == insideSub);
        assert(getSub(outsideIdx) == outsideSub);
        assert(getDim(insideIdx) == insideDim);
        assert(getDim(outsideIdx) == outsideDim);
        break;
    }
    case 5:
    case 7:
    {
        // Edge/cell.
        auto edgeIdx = (insideDim == 1 ? insideIdx : outsideIdx);
        auto cellIdx = (insideDim == 1 ? outsideIdx : insideIdx);
        assert(cellIdx == 8 || cellIdx == 9 + edgeIdx);
        auto edgeBit = NeighborIndex(edgeIdx).floating();
        auto otherBit = 3 ^ edgeBit;
        bool isHighEdge(NeighborIndex(edgeIdx).pos() & otherBit);
        if (cellIdx != 8) {
            isHighEdge = !isHighEdge;
        }
        auto highIsCC = (isHighEdge == (edgeBit == Axis::Y));
        auto highIsForward = (highIsCC == (insideDim == 2));
        auto cornerPos = NeighborIndex(edgeIdx).pos();
        if (highIsForward) {
            cornerPos |= edgeBit;
        }
        thirdIdx = NeighborIndex::fromPosAndFloating(cornerPos, 0).i;
        break;
    }
    default:
        assert(false);
        return nullptr;
    }
    auto thirdSub = getSub(thirdIdx);
    if (thirdSub == insideSub) {
        return forwardSimplex(cell, thirdIdx, outsideIdx);
    }
    else if (thirdSub == outsideSub) {
        return forwardSimplex(cell, insideIdx, thirdIdx);
    }
    else {
        assert(cell->leaf);
        auto edgeIdx = insideDim == 1 ? insideIdx : 
                       outsideDim == 1 ? outsideIdx :
                       thirdIdx;
        auto cellIdx = insideDim == 2 ? insideIdx :
                       outsideDim == 2 ? outsideIdx :
                       thirdIdx;
        auto cornerIdx = insideDim == 0 ? insideIdx :
                         outsideDim == 0 ? outsideIdx :
                         thirdIdx;
        auto edgeAxis = static_cast<Axis::Axis>(
            NeighborIndex(edgeIdx).floating());
        assert(edgeAxis == Axis::X || edgeAxis == Axis::Y);
        auto edgePos = NeighborIndex(edgeIdx).pos();
        const auto& edgeVar = cell->leaf->edge(edgeAxis, edgePos);
        assert(std::holds_alternative<SimplexDCMinEdge<2>*>(edgeVar));
        const auto& edge = std::get<SimplexDCMinEdge<2>*>(edgeVar);
        bool isUpperCorner(NeighborIndex(cornerIdx).pos() & edgeAxis);
        bool isUpperCell(!edgePos);
        if (cellIdx != 8) {
            assert(cellIdx == 9 + edgeIdx);
            isUpperCell = !isUpperCell;
        }
        auto& simplex = edge->simplex(
            edgeAxis, Axis::Z, isUpperCell, isUpperCorner);
        return &simplex;
    }

    return nullptr;
}

}   // namespace libfive
