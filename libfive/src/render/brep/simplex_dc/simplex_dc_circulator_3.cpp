/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "libfive/render/brep/simplex_dc/simplex_dc_circulator.hpp"

namespace libfive {
template class SimplexDCCirculator<3, false>;

SimplexDCCirculator<3, true>::SimplexDCCirculator(
    FarFaceTag, const SimplexDCTree<3>* cell, unsigned nonFaceIdx, 
    Axis::Axis edgeAxis, Axis::Axis faceAxis,
    bool faceIsOutside, bool allowDegenerate)
{
    assert(edgeAxis != faceAxis);
    this->edgeCell = cell;
    this->edgeAxis = edgeAxis;
    auto edgePos = NeighborIndex(nonFaceIdx).pos();
    NeighborIndex testNeighbor(nonFaceIdx);
    auto otherDim = NeighborIndex(nonFaceIdx).dimension();
    if (otherDim == 1) {
        assert(edgeAxis == NeighborIndex(nonFaceIdx).floating());
    }
    else if (otherDim == 0) {
        testNeighbor = NeighborIndex::fromPosAndFloating(edgePos, edgeAxis);
    }
    else {
        assert(false);
        hitError = true;
        return;
    }
    if (!allowDegenerate && cell->neighbor(testNeighbor) == nullptr) {
        // Our edge or corner and our face are both on the border of our 
        // original bounding box, but one is inside and one outside; this
        // is an invalid bounding box, which we are not currently equipped to
        // handle.
        assert(false);
        hitError = true;
        return;
    }
    this->edgeRelativePos = ((edgePos * 9) / (2 * edgeAxis)) & 3;
    this->faceAxis = faceAxis;
    this->relativeCellPos = edgeRelativePos; // Opposite the edge from "cell",
                                             // since we want the far face.
    this->cornerPos = edgePos & edgeAxis;

    this->insideDimension = faceIsOutside ? otherDim : 2;
    this->outsideDimension = faceIsOutside ? 2 : otherDim;
    this->allowDegenerate = allowDegenerate;
    this->allowCellsPastBBox = allowDegenerate;
    while (!isValid() && !hitError) {
        advance();
    }
    startSimplex = simplex();
}

void SimplexDCCirculator<3, true>::construct(
    const SimplexDCTree<3>* cell,
    std::array<unsigned, 2> indexPr)
{
    NeighborIndex index0(indexPr[0]);
    NeighborIndex index1(indexPr[1]);
    
    auto getDimAndSetNeighbor = [](NeighborIndex& idx) {
        if (idx.i < 27) {
            return idx.dimension();
        }
        else {
            idx.i -= 27;
            return 3u;
        }
    };

    auto dim0 = getDimAndSetNeighbor(index0);
    auto dim1 = getDimAndSetNeighbor(index1);
    assert(dim0 != dim1);
    if (dim0 > dim1) {
        std::swap(index0, index1);
        std::swap(dim0, dim1);
    }
    switch (4 * dim0 + dim1) {
    case 1:
    {
        // Vertex/edge
        edgeCell = cell;
        auto edgeFloating = index1.floating();
        edgeAxis = Axis::Axis(edgeFloating);
        for (auto toFloat = 0; toFloat <= 6; ++toFloat) {
            if (toFloat & edgeAxis) {
                // Use a lambda so that in release 'neighbor' is not calculated.
                auto badNeighbor = [&]() {
                    auto neighbor = cell->neighbor(
                        NeighborIndex::fromPosAndFloating(index0.pos(), toFloat));
                    return neighbor && neighbor->isBranch();
                };
                assert(!badNeighbor());
            }
        }
        auto pos = index1.pos();
        assert(((pos ^ index0.pos()) & ~edgeFloating) == 0);
        pos = pos >> 1 | (pos << 2);
        pos >>= Axis::toIndex(edgeAxis);
        edgeRelativePos = pos & 3;
        faceAxis = Axis::Q(edgeAxis);
        assert(index0.floating() == 0);
        cornerPos = index0.pos() & edgeFloating;
        relativeCellPos = 0;
        if (!allowCellsPastBBox && !getSub<3>()) {
            assert(false);
            hitError = true;
        }
        break;
    }
    case 2:
    {
        // Vertex/face
        auto faceFloating = index1.floating();
        faceAxis = Axis::Axis(7 & ~faceFloating);
        edgeAxis = Axis::Q(faceAxis);
        assert(index0.floating() == 0);
        auto cornerValue = index0.pos() & edgeAxis;
        cornerPos = cornerValue;
        auto RValue = index0.pos() & Axis::R(faceAxis);
        auto faceValue = index0.pos() & faceAxis;
        assert(faceValue == (index1.pos() & faceAxis));
        auto edgePos = RValue | faceValue;

        // Check everything adjacent to our edge to see which cell we actually
        // want to use; start with the current cell.
        bool bestR(RValue);
        bool bestFace(faceValue);
        auto bestEdgeCell = cell;
        for (auto toFloat = 0; toFloat <= 6; ++toFloat) {
            if (!(toFloat & edgeAxis)) {
                // Not next to our edge regardless of which direction we go.
                continue;
            }
            auto direction = 
                NeighborIndex::fromPosAndFloating(edgePos, toFloat);
            auto candidateCell = cell->neighbor(direction);
            if (!candidateCell) {
                continue;
            }
            auto childIdx = 7 ^ toFloat ^ index0.pos();
            while (candidateCell->isBranch()) {
                candidateCell = candidateCell->child(childIdx);
            }
            if (candidateCell->leafLevel() < bestEdgeCell->leafLevel()) {
                bestEdgeCell = candidateCell;
                bestFace = bool(faceValue) == bool(toFloat & faceAxis);
                bestR = bool(RValue) == bool(toFloat & Axis::R(faceAxis));
            }
        }
        edgeCell = bestEdgeCell;
        edgeRelativePos = int(bestR) | 2 * bestFace;
        relativeCellPos = RValue ? 0 : 3;
        break;
    }
    case 3:
    {
        // Vertex/cell.  Index1 was set to refer to the subspace in
        // the direction of the simplex-containing cell (not necessarily
        // edgeCell).
        assert(index0.floating() == 0);
        assert(((index0.pos() ^ index1.pos()) & ~index1.floating()) == 0);
        auto simplexCell = cell->neighbor(index1);
        assert(simplexCell);
        // We need to choose an edge that actually exists, i.e. simplexCell
        // and (if different) its opposite do not span all four neighbors in 
        // that direction.  We don't need to worry about whether the face 
        // exists; as long as we have an edge to reference, the circulator 
        // can determine that it is invalid and will advance automatically.

        auto cornerPosFromCell = index0.pos();
        auto cellFloat = index1.floating();
        auto isValidEdgeAxis = [&](Axis::Axis edge) {
            if (cellFloat & edge) {
                // Since our vertex is a vertex of 'cell', it must be
                // adjacent to an edge of 'cell' in each direction.  So
                // if cellFloat & edge, that edge of 'cell' is also the
                // edge of that axis adjacent to that vertex and to 
                // 'simplexCell', and so the edge of that axis must 
                // actually exist.
                return true;
            }
            auto opposite = cell->neighbor(
                NeighborIndex::fromPosAndFloating(
                    cornerPosFromCell, cellFloat ^ 7 ^ edge));
            if (opposite == simplexCell) {
                return false;
            }
            auto hadANonNull = false;
            for (auto testAxis : { Axis::Q(edge), Axis::R(edge) }) {
                auto testCell = cell->neighbor(
                    NeighborIndex::fromPosAndFloating(
                        cornerPosFromCell, cellFloat ^ testAxis));
                if (testCell != simplexCell && testCell != opposite) {
                    return true;
                }
                else if (testCell) {
                    hadANonNull = true;
                }
                else {
                    assert(!opposite);
                }
            }
            // If we reached here, we had only two cells among our four
            // positions along the edge.  That makes it an invalid edge,
            // unless we had one real cell at simplexCell and the other
            // three were non-null.
            return !hadANonNull;
        };

        if (isValidEdgeAxis(Axis::X)) {
            edgeAxis = Axis::X;
        }
        else if (isValidEdgeAxis(Axis::Y)) {
            edgeAxis = Axis::Y;
        }
        else {
            assert(isValidEdgeAxis(Axis::Z));
            edgeAxis = Axis::Z;
        }

        faceAxis = Axis::Q(edgeAxis);

        // Now that we know which edge axis we're using, the edge we're using 
        // must be the one next to simplexCell of that axis with our corner
        // as an endpoint; if the corner is on an edge or face of simplexCell,
        // there may be two such edges, and either will do.  We thus want to 
        // find the smallest cell next to that edge.
        // We will iterate through neighbors of cell, not of simplexCell, since
        // cell is guaranteed to have our corner as one of its corners and thus
        // makes it easy to find the correct neighbors.


        auto cornerFromEdge = (index0.pos() ^ index1.fixed()) & edgeAxis;
        cornerPos = bool(cornerFromEdge);
        auto simplexCellFromEdge = 
            (index0.pos() ^ index1.floating()) & (~edgeAxis);
        auto bestCellFromEdge = simplexCellFromEdge;
        auto bestEdgeCell = simplexCell;
        bool edgeNextToCell(index1.floating() & edgeAxis);
        auto simplexRelativePosQ =
            (simplexCellFromEdge & Axis::Q(edgeAxis)) ? 1 : 0;
        auto simplexRelativePosR =
            (simplexCellFromEdge & Axis::R(edgeAxis)) ? 1 : 0;
        this->relativeCellPos = simplexRelativePosQ | 2 * simplexRelativePosR;
        for (auto toFloat = 0; toFloat <= 7; ++toFloat) {
            if (bool(toFloat & edgeAxis) != edgeNextToCell) {
                // Not next to our edge regardless of which direction we go.
                continue;
            }
            auto direction =
                NeighborIndex::fromPosAndFloating(index0.pos(), toFloat);
            auto candidateCell = cell->neighbor(direction);
            if (!candidateCell || candidateCell == simplexCell) {
                continue;
            }
            auto childIdx = index0.pos() ^ toFloat ^ 7;
            while (candidateCell->isBranch()) {
                candidateCell = candidateCell->child(childIdx);
            }
            if (candidateCell->leafLevel() < bestEdgeCell->leafLevel()) {
                bestEdgeCell = candidateCell;
                bestCellFromEdge = 
                    (index0.pos() ^ toFloat) & ~edgeAxis;
            }
        }
        edgeCell = bestEdgeCell;
        auto bestEdgeRelativePosQ = 
            (bestCellFromEdge & Axis::Q(edgeAxis)) ? 0 : 1;
        auto bestEdgeRelativePosR =
            (bestCellFromEdge & Axis::R(edgeAxis)) ? 0 : 1;
        edgeRelativePos = bestEdgeRelativePosQ | 2 * bestEdgeRelativePosR;
        break;
    }
    case 6:
    {
        // Edge/face
        edgeCell = cell;
        auto edgeFloating = index0.floating();
        auto faceFloating = index1.floating();
        assert(edgeFloating & faceFloating);
        edgeAxis = Axis::Axis(edgeFloating);
        faceAxis = Axis::Axis(7 & ~faceFloating);
        assert(edgeAxis != faceAxis);

        for (auto toFloat = 0; toFloat <= 6; ++toFloat) {
            if (toFloat & edgeAxis) {
                // Use a lambda so that in release 'neighbor' is not calculated.
                auto badNeighbor = [&]() {
                    auto neighbor = cell->neighbor(
                        NeighborIndex::fromPosAndFloating(index0.pos(), toFloat));
                    return neighbor && neighbor->isBranch();
                };
                assert(!badNeighbor());
            }
        }

        auto pos = index0.pos();
        assert(((pos ^ index1.pos()) & ~faceFloating) == 0);
        pos = pos >> 1 | (pos << 2);
        pos >>= Axis::toIndex(edgeAxis);
        edgeRelativePos = pos & 3;
        cornerPos = 0;
        relativeCellPos = 3 ^ edgeRelativePos;
        break;
    }
    case 7:
    {
        // Edge/cell
        edgeCell = cell;
        auto edgeFloating = index0.floating();
        assert(edgeFloating & index1.floating());
        edgeAxis = Axis::Axis(edgeFloating);
        for (auto toFloat = 0; toFloat <= 6; ++toFloat) {
            if (toFloat & edgeAxis) {
                // Use a lambda so that in release 'neighbor' is not calculated.
                auto badNeighbor = [&]() {
                    auto neighbor = cell->neighbor(
                        NeighborIndex::fromPosAndFloating(index0.pos(), toFloat));
                    return neighbor && neighbor->isBranch();
                };
                assert(!badNeighbor());
            }
        }

        auto pos = index0.pos();
        assert(((pos ^ index1.pos()) & ~index1.floating()) == 0);
        auto reduceByEdge = [this](unsigned char& pos) {
            pos = pos >> 1 | (pos << 2);
            pos >>= Axis::toIndex(edgeAxis);
            pos &= 3;
        };
        reduceByEdge(pos);
        edgeRelativePos = pos;

        auto simplexCellFloating = index1.floating();
        reduceByEdge(simplexCellFloating);
        
        relativeCellPos = pos ^ simplexCellFloating;

        // Again, we don't need to worry about whether our face exists; if it
        // doesn't, we'll advance from there.
        faceAxis = Axis::Q(edgeAxis);

        cornerPos = 0;
        break;
    }
    case 11:
    {
        // Face/cell
        auto faceFloating = index0.floating();
        faceAxis = Axis::Axis(7 & ~faceFloating);
        assert(index1.i == 26 || index1.i == index0.i);
        assert(!cell->neighbor(index1)->isBranch());
        edgeAxis = Axis::R(faceAxis);
        auto thirdAxis = Axis::Q(faceAxis);
        assert(thirdAxis == Axis::R(edgeAxis));
        auto faceValue = index0.pos();

        // Check everything adjacent to our edge to see which cell we actually
        // want to use; start with the current cell.
        auto bestThird = true; // Whether the best cell is in the +third
                               // direction of the edge.
        auto bestFace = !faceValue; // Whether the best cell is in the +face
                                    // direction of the edge.
        auto bestEdgeCell = cell;
        for (auto toFloat = 0; toFloat <= 6; ++toFloat) {
            if (!(toFloat & edgeAxis)) {
                // Not next to our edge regardless of which direction we go.
                continue;
            }
            auto direction =
                NeighborIndex::fromPosAndFloating(faceValue, toFloat);
            auto candidateCell = cell->neighbor(direction);
            if (!candidateCell) {
                continue;
            }
            auto childIdx = 7 ^ faceValue ^ toFloat;
            while (candidateCell->isBranch()) {
                candidateCell = candidateCell->child(childIdx);
            }
            if (candidateCell->leafLevel() < bestEdgeCell->leafLevel()) {
                bestEdgeCell = candidateCell;
                bestFace = bool(faceValue) != bool(toFloat & faceAxis);
                bestThird = bool(toFloat & thirdAxis);
            }
        }
        assert(bestEdgeCell);
        edgeCell = bestEdgeCell;
        edgeRelativePos = int(!bestFace) | 2 * !bestThird;
        cornerPos = 0;
        if (bool(faceValue) == (index1.i == 26)) {
            // Our cell is lower in faceValue as compared to the edge...and must
            // be upper in thirdValue (since the edge is always on the 
            // -thirdValue side of the starting cell).
            relativeCellPos = 2;
        }
        else {
            // It's upper in faceValue, and still must be upper in thirdValue.
            relativeCellPos = 3;
        }
        break;
    }
    default:
        assert(false);
    }
}

void SimplexDCCirculator<3, true>::advance()
{
    if (hitError) {
        return;
    }
    assert(insideDimension != outsideDimension);
    // Our first step is to get the dimensions of the "forward" and "backward"
    // vertices.
    auto iofb = getIOFB();
    assert(iofb[0] == insideDimension);
    assert(iofb[1] == outsideDimension);
    auto forwardDimension = iofb[2];
    auto backwardDimension = iofb[3];

    auto forwardSub = getSub(forwardDimension);
    if (!forwardSub) {
        assert(forwardDimension == 2 || 
               (forwardDimension == 3 && allowCellsPastBBox)); 
                // Otherwise, it must exist.
    }

    auto dimensionToSwap = backwardDimension;
    if (!forwardSub) {
        if (forwardDimension == 2) {
            // This means that we just advanced across a face that doesn't
            // exist; there are thus three possibilities for the simplex edge
            // we are circulating around.  It may be between a cell sub and an
            // edge on its face; in such case, what we actually want to do is 
            // swap the cell, in order to reach the next "version" of our  
            // "duplicate" cell.  Because our face does not exist, swapping the
            // cell across it will not actually change the identity of our cell
            // or its sub, but will change which dimension is forward, so that
            // the next call to advance() will swap it to the next face without
            // changing the corner, and then the following one will move on to 
            // the other corner.  
            // Alternatively, our simplex edge may be between a cell sub and a
            // corner that is on either its edge or its face.  In such case, we
            // again want to swap the cell, not actually changing the cell
            // identity but causing the next advance() to move on to the next
            // face without changing the edge, and the one after that to swap
            // the edge.  
            // Finally, our simplex edge may be between an edge sub and a
            // corner, where that edge is on a face of a larger cell.  In such
            // case, we again want to swap the cell, though in that case this is
            // what would be done anyway.
            dimensionToSwap = 3;
        }
    }
    else if (!allowDegenerate) {
        auto insideSub = getSub(insideDimension);
        auto outsideSub = getSub(outsideDimension);
        if (insideSub == outsideSub || !insideSub || !outsideSub) {
            // We're circulating around an invalid simplex edge.  (Even if we
            // allow subs outside the bounding box which will therefore be
            // nullptr, only forwardSub and backwardSub are allowed to be such,
            // not the two we're cycling around.)
            assert(false);
            hitError = true;
            return;
        }
        auto backSubIsDuplicateOrInvalid = [&]() {
            auto backwardSub = getSub(backwardDimension);
            if (!backwardSub) {
                return true;
            }
            return backwardSub == insideSub || backwardSub == outsideSub;
        };
        if (forwardSub == insideSub) {
            if (insideDimension < forwardDimension || 
                !backSubIsDuplicateOrInvalid()) {
                insideDimension = forwardDimension;
                return;
            };
        }
        else if (forwardSub == outsideSub && 
                 (outsideDimension < forwardDimension || 
                  !backSubIsDuplicateOrInvalid())) {
            outsideDimension = forwardDimension;
            return;
        }
    }

    switch (dimensionToSwap) {
    case 0:
    {
        // We just need to swap to the other corner adjacent to this edge.
        cornerPos = !cornerPos;
        break;
    }
    case 2:
    {
        // We need to swap our face axis, thereby swapping the face around
        // the cell while keeping the same edge.  This may result in an
        // invalid face, but that's not an issue; it'll just be detected
        // as invalid and we'll advance again.
        assert(faceAxis == Axis::Q(edgeAxis) || faceAxis == Axis::R(edgeAxis));
        faceAxis = Axis::Axis(7 & ~faceAxis & ~edgeAxis);
        break;
    }
    case 3:
    {
        // We need to swap the cell across our face, and then we're done.
        auto relativeFaceAxis = (faceAxis == Axis::Q(edgeAxis) ? 1 : 2);
        relativeCellPos ^= relativeFaceAxis;
        if (!allowCellsPastBBox && !getSub<3>()) {
            assert(false);
            hitError = true;
            return;
        }
        break;
    }
    case 1:
    {
        // We need to swap our edge dimension.  This is where things get tricky.
        // We need to swap between the two edges adjacent to our corner 
        // and face (which must exist, as otherwise we shouldn't enter this 
        // switch.)
        auto newEdgeAxis = Axis::Axis(7 & ~faceAxis & ~edgeAxis);
        auto newEdgeRelativeMask = (faceAxis == Axis::Q(edgeAxis) ? 2 : 1);
        bool stillAdjacentToOldEdgeCell(
            (edgeRelativePos^ relativeCellPos)& newEdgeRelativeMask);
        const SimplexDCTree<3>* bestCell = nullptr;
        unsigned bestFloat;
        auto bestCellisDuplicate = false;
        auto absoluteEdgeDir = edgeRelativePos * edgeAxis * 2;
        absoluteEdgeDir |= absoluteEdgeDir >> 3;
        absoluteEdgeDir &= 7;
        auto cornerDirection = (edgeAxis * cornerPos) | absoluteEdgeDir;
        auto simplexCellFromCorner = relativeCellPos * (edgeAxis << 1);
        simplexCellFromCorner |= simplexCellFromCorner >> 3;
        if (!cornerPos) {
            simplexCellFromCorner |= edgeAxis;
        }
        simplexCellFromCorner &= 7;

        for (auto toFloat = 0; toFloat <= 7; ++toFloat) {
            if (bool(toFloat & newEdgeAxis) != stillAdjacentToOldEdgeCell) {
                continue; // Not adjacent to our new edge.
            }
            auto candidateCell = edgeCell->neighbor(
                NeighborIndex::fromPosAndFloating(
                    cornerDirection, toFloat));
            if (!candidateCell) {
                continue;
            }
            while (candidateCell->isBranch()) {
                auto childIdx = cornerDirection ^ (7 & ~toFloat);
                candidateCell = candidateCell->child(childIdx);
            }
            if (candidateCell == bestCell) {
                bestCellisDuplicate = true;
                continue;
            }
            if (bestCell &&
                bestCell->leafLevel() <= candidateCell->leafLevel()) {
                continue;
            }
            bestCell = candidateCell;
            bestFloat = toFloat;
            bestCellisDuplicate = false;
        }
        assert(bestCell);
        if (!bestCellisDuplicate) {
            // This is the normal situation, in which our corner is at a
            // corner of our face.
            auto newCornerDirection = cornerDirection ^ (7 & ~bestFloat);
            edgeCell = bestCell;
            edgeAxis = newEdgeAxis;
            edgeRelativePos = (newCornerDirection * 9) / (edgeAxis << 1);
            edgeRelativePos &= 3;
            cornerPos = newCornerDirection & edgeAxis;
            relativeCellPos = (simplexCellFromCorner * 9) / (edgeAxis << 1);
            relativeCellPos &= 3;
        }
        else {
            // The new edge can't be adjacent to the old edge cell, as
            // if it were, the old edge cell would (due to having our
            // corner as a corner) necessarily be better than any cell
            // that does not have it as a corner (the implication of
            // being a duplicate).
            // And our best cell can't be duplicate in the direction of
            // faceAxis, since if it were then to be the best cell all
            // four candidates next to the edge would need to be duplicate
            // in that direction, but one of those two pairs has our face
            // in between, and we know that our face exists (as we branched
            // on that condition).
            // Thus:
            assert(bestCell == edgeCell->neighbor(
                NeighborIndex::fromPosAndFloating(
                    cornerDirection, bestFloat ^ edgeAxis)));
            // This means that our corner was on an edge of our face,
            // and the other edge adjacent to both is of the same axis
            // but on the other side of our corner with respect to that
            // axis.
            bestCell = nullptr;
            bestCellisDuplicate = false; // For asserts only at this point;
                                         // our situation does not allow for
                                         // the new edge to be duplicate.
            cornerPos = !cornerPos;
            for (auto toFloat = 0; toFloat <= 6; ++toFloat) {
                if (toFloat & edgeAxis) {
                    continue; // Adjacent to our old edge, 
                              // hence not to the new one.
                }
                auto candidateCell = edgeCell->neighbor(
                    NeighborIndex::fromPosAndFloating(
                        cornerDirection, toFloat));
                if (!candidateCell) {
                    continue;
                }
                while (candidateCell->isBranch()) {
                    auto childIdx = cornerDirection ^ (7 & ~toFloat);
                    candidateCell = candidateCell->child(childIdx);
                }
                if (candidateCell == bestCell) {
                    bestCellisDuplicate = true;
                    continue;
                }
                if (bestCell &&
                    bestCell->leafLevel() <= candidateCell->leafLevel()) {
                    continue;
                }
                bestCell = candidateCell;
                bestFloat = toFloat;
                bestCellisDuplicate = false;
            }
            assert(bestCell);
            assert(!bestCellisDuplicate);
            edgeCell = bestCell;
            if (!(bestFloat & Axis::Q(edgeAxis))) {
                edgeRelativePos ^= 1;
            }
            if (!(bestFloat & Axis::R(edgeAxis))) {
                edgeRelativePos ^= 2;
            }
        }
        break;
    }
    default:
        assert(false);
        hitError = true;
        return;
    }
}

bool SimplexDCCirculator<3, true>::isValid() const
{
    if (hitError) {
        return false;
    }
    auto faceSub = getSub<2>();
    if (!faceSub) {
        return false;
    }
    auto cellSub = getSub<3>();
    if (!cellSub) {
        assert(allowCellsPastBBox);
        return false; // It's not an error, but still not a valid simplex.
    }
    else if (allowDegenerate) {
        return true;
    }
    else if (cellSub == faceSub) {
        return false;
    }
    auto cornerSub = getSub<0>();
    assert(cornerSub);
    if (cornerSub == faceSub || cornerSub == cellSub) {
        return false;
    }
    auto edgeSub = getSub<1>();
    assert(edgeSub);
    if (edgeSub == cornerSub || edgeSub == faceSub || edgeSub == cellSub) {
        return false;
    }
    return true;
}

std::array<int, 4> SimplexDCCirculator<3, true>::getIOFB() const
{
    assert(insideDimension != outsideDimension);
    // We do this by determining if inside / outside / forward / backward
    // is an even or odd permutation of 0/1/2/3, since that is dependent only
    // on what simplex we're using.

    auto QFaceIsCCLookingFromPos =
        (relativeCellPos == 0 || relativeCellPos == 3);
    auto faceIsQ = (faceAxis == Axis::Q(edgeAxis));
    auto cornerIsPositiveOnEdge = cornerPos;
    auto QFaceIsCCLookingFromCorner =
        (QFaceIsCCLookingFromPos == cornerIsPositiveOnEdge);
    auto faceIsForwardForCornerOutside =
        (QFaceIsCCLookingFromCorner == faceIsQ);
    auto permutationIsOdd = faceIsForwardForCornerOutside;

    // Now we use swaps in an array to find the actual forward and backward 
    // dimensions, adjusting permutationIsOdd accordingly.  We need to end up
    // with an even permutation that has the correct inside and outside, as
    // that will then be the trivial permutation so we can read our dimensions
    // off of the array.
    std::array iofb{ 0, 1, 2, 3 };
    if (insideDimension != 0) {
        permutationIsOdd = !permutationIsOdd;
        std::swap(iofb[0], iofb[insideDimension]);
    }
    if (outsideDimension != iofb[1]) {
        permutationIsOdd = !permutationIsOdd;
        auto swapWith = (outsideDimension == iofb[2] ? 2 : 3);
        std::swap(iofb[1], iofb[swapWith]);
    }
    assert(insideDimension == iofb[0]);
    assert(outsideDimension == iofb[1]);
    if (permutationIsOdd) {
        permutationIsOdd = false;
        std::swap(iofb[2], iofb[3]);
    }
    return iofb;
}

void SimplexDCCirculator<3, false>::insertIntersection(
    const SimplexDCIntersection<3>* intersection) const
{
    auto mySimplex = simplex();
    if (!mySimplex) {
        assert(false);
        return;
    }
    [[maybe_unused]] auto [res, success] = 
        mySimplex->insertIntersection(
        insideDimension, outsideDimension, intersection);
    assert(res == intersection);
    assert(success);
}

template <>
SimplexLeafSubspace<3>* SimplexDCCirculator<3, true>::getSub<0>() const
{
    if (hitError) {
        return nullptr;
    }
    auto pos = edgeRelativePos * (edgeAxis * 2);
    pos |= (edgeRelativePos * edgeAxis / 4);
    if (cornerPos) {
        pos |= edgeAxis;
    }
    pos &= 7;
    return edgeCell->leaf->collapsedSub(
        NeighborIndex::fromPosAndFloating(pos, 0).i);
}

template <>
SimplexLeafSubspace<3>* SimplexDCCirculator<3, true>::getSub<1>() const
{
    if (hitError) {
        return nullptr;
    }
    auto pos = edgeRelativePos * (edgeAxis * 2);
    pos |= (edgeRelativePos * edgeAxis / 4);
    pos &= 7;
    return edgeCell->leaf->collapsedSub(
        NeighborIndex::fromPosAndFloating(pos, edgeAxis).i);
}

template <>
SimplexLeafSubspace<3>* SimplexDCCirculator<3, true>::getSub<2>() const
{
    if (hitError) {
        return nullptr;
    }
    auto relativeFloatToSimplexCell = edgeRelativePos ^ relativeCellPos;
    auto relativeToAbsolute = [this](int relative) {
        auto absolute = relative * 2 * edgeAxis;
        absolute |= absolute >> 3;
        absolute &= 7;
        return absolute;
    };
    auto floatToSimplexCell = 
        edgeAxis | relativeToAbsolute(relativeFloatToSimplexCell);
    auto floatToNeighborCell = floatToSimplexCell ^ faceAxis;
    auto posToSimplexAndNeighbor = relativeToAbsolute(edgeRelativePos);

    auto simplexCell = edgeCell->neighbor(
        NeighborIndex::fromPosAndFloating(
            posToSimplexAndNeighbor, floatToSimplexCell));
    if (!simplexCell) {
        assert(allowCellsPastBBox);
    }
    auto neighborCell = edgeCell->neighbor(
        NeighborIndex::fromPosAndFloating(
            posToSimplexAndNeighbor, floatToNeighborCell));
    if (simplexCell == neighborCell) {
        // The face does not actually exist.
        return nullptr;
    }
    auto neighborIsHigher = 
        (posToSimplexAndNeighbor ^ floatToNeighborCell) & faceAxis;
    auto faceFromHigher = NeighborIndex::fromPosAndFloating(0, 7 ^ faceAxis);
    auto faceFromLower = NeighborIndex::fromPosAndFloating(7, 7 ^ faceAxis);
    auto faceFromSimplexCell = 
        neighborIsHigher ? faceFromLower : faceFromHigher;
    auto faceFromNeighbor = neighborIsHigher ? faceFromHigher : faceFromLower;
    if (!simplexCell) {
        if (!neighborCell) {
            return nullptr;
        }
        else {
            return neighborCell->leaf->collapsedSub(faceFromNeighbor.i);
        }
    }
    else if (neighborCell && 
        neighborCell->leafLevel() < simplexCell->leafLevel()) {
        return neighborCell->leaf->collapsedSub(faceFromNeighbor.i);
    }
    else {
        return simplexCell->leaf->collapsedSub(faceFromSimplexCell.i);
    }
}

template <>
SimplexLeafSubspace<3>* SimplexDCCirculator<3, true>::getSub<3>() const
{
    if (hitError) {
        return nullptr;
    }
    auto absoluteFromRelative = [this](int relative) {
        auto out = relative * (edgeAxis * 2);
        out |= (relative * edgeAxis / 4);
        out &= 7;
        return out;
    };

    auto simplexCellPos = absoluteFromRelative(edgeRelativePos);
    auto simplexCellFloat =
        absoluteFromRelative(edgeRelativePos ^ relativeCellPos) | edgeAxis;
    auto simplexCell = edgeCell->neighbor(
        NeighborIndex::fromPosAndFloating(simplexCellPos, simplexCellFloat));
    if (!simplexCell) {
        assert(allowCellsPastBBox);
        // It is the caller's responsibility to ensure that hitError is set to
        // true, since this method is const.
        return nullptr;
    }
    return simplexCell->leaf->collapsedSub(26);
}

}   // namespace libfive
