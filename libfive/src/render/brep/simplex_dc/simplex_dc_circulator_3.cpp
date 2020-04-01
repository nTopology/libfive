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
    switch (dim0 + 4 * dim1) {
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
        assert((pos ^ index0.pos()) & edgeFloating == 0);
        pos = pos >> 1 | (pos << 2);
        pos >>= Axis::toIndex(edgeAxis);
        edgeRelativePos = pos & 3;
        faceAxis = Axis::Q(edgeAxis);
        assert(index0.floating() == 0);
        cornerPos = index0.pos() & edgeFloating;
        relativeCellPos = 0;
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
        assert(faceValue == index1.pos() & faceAxis);
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
            auto childIdx = (7 ^ edgePos) | cornerValue;
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
        assert(index0.pos() == 0);
        assert((index0.pos() ^ index1.pos()) & ~index1.floating() == 0);
        auto simplexCell = cell->neighbor(index1);
        assert(simplexCell);
        const SimplexDCTree<3>* neighborCell;
        // We need to choose an edge that actually exists, i.e. simplexCell
        // and (if different) its opposite do not span all four neighbors in 
        // that direction.  We don't need to worry about whether the face 
        // exists; as long as we have an edge to reference, the circulator 
        // can determine that it is invalid and will advance automatically.

        auto relativeCellPos = index1.pos();
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
                    relativeCellPos, cellFloat ^ 7 ^ edge));
            if (opposite == simplexCell) {
                return false;
            }
            auto hadANonNull = false;
            for (auto testAxis : { Axis::Q(edge), Axis::R(edge) }) {
                auto testCell = cell->neighbor(
                    NeighborIndex::fromPosAndFloating(
                        relativeCellPos, cellFloat ^ testAxis));
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
        // must be the one next to simplexCell of that axis with our vertex
        // as an endpoint, so we can find the smallest cell next to that edge.
        auto cornerValue = (index0.pos() ^ index1.fixed()) & edgeAxis;
        cornerPos = cornerValue;
        auto faceValue = (index0.pos() ^ index1.fixed()) & faceAxis;
        auto RValue = (index0.pos() ^ index1.fixed()) & Axis::R(edgeAxis);
        bool bestR(RValue);
        bool bestFace(faceValue);
        auto bestEdgeCell = simplexCell;
        auto edgePos = RValue | faceValue;
        relativeCellPos = !bestFace | 2 * !bestR;
        for (auto toFloat = 0; toFloat <= 6; ++toFloat) {
            if (!(toFloat & edgeAxis)) {
                // Not next to our edge regardless of which direction we go.
                continue;
            }
            auto direction =
                NeighborIndex::fromPosAndFloating(edgePos, toFloat);
            auto candidateCell = simplexCell->neighbor(direction);
            auto childIdx = (7 ^ edgePos) | cornerValue;
            while (candidateCell->isBranch()) {
                candidateCell = candidateCell->child(childIdx);
            }
            if (candidateCell->leafLevel() < bestEdgeCell->leafLevel()) {
                bestEdgeCell = candidateCell;
                bestFace = bool(faceValue) == bool(toFloat & faceAxis);
                bestR = bool(RValue) == bool(toFloat & Axis::R(edgeAxis));
            }
        }
        edgeCell = bestEdgeCell;
        edgeRelativePos = int(bestFace) | 2 * bestR;
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
        assert((pos ^ index1.pos()) & ~faceFloating == 0);
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
        assert((pos ^ index1.pos()) & ~index1.floating() == 0);
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
        auto bestThird = false;
        bool bestFace(faceValue);
        auto bestEdgeCell = cell;
        auto bestUnreducedCellPos = index1.pos();
        for (auto toFloat = 0; toFloat <= 6; ++toFloat) {
            if (!(toFloat & edgeAxis)) {
                // Not next to our edge regardless of which direction we go.
                continue;
            }
            auto direction =
                NeighborIndex::fromPosAndFloating(faceValue, toFloat);
            auto candidateCell = cell->neighbor(direction);
            auto childIdx = 7 ^ faceValue ^ edgeAxis;
            while (candidateCell->isBranch()) {
                candidateCell = candidateCell->child(childIdx);
            }
            if (candidateCell->leafLevel() < bestEdgeCell->leafLevel()) {
                bestEdgeCell = candidateCell;
                bestFace = bool(faceValue) == bool(toFloat & faceAxis);
                bestThird = !(bool(toFloat & thirdAxis));
                bestUnreducedCellPos = index1.pos() ^ (7 & ~toFloat);
            }
        }
        assert(bestEdgeCell);
        edgeCell = bestEdgeCell;
        edgeRelativePos = int(bestFace) | 2 * bestThird;
        cornerPos = 0;
        relativeCellPos = int(bool(bestUnreducedCellPos & faceAxis)) | 
                  2 * bool(bestUnreducedCellPos & thirdAxis);
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
    // vertices.  We do this by determining if inside/outside/forward/backward
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
    auto forwardDimension = iofb[2];
    auto backwardDimension = iofb[3];

    auto insideSub = getSub(insideDimension);
    auto outsideSub = getSub(outsideDimension);
    if (insideSub == outsideSub || !insideSub || !outsideSub) {
        // We're circulating around an invalid simplex edge.
        hitError = true;
        return;
    }
    auto forwardSub = getSub(forwardDimension);
    if (forwardSub == insideSub) {
        insideDimension = forwardDimension;
        return;
    }
    else if (forwardSub == outsideSub) {
        outsideDimension = forwardDimension;
        return;
    }
    else if (!forwardSub) {
        assert(forwardDimension == 2); // All others must always exist.
    }

    switch (backwardDimension) {
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
        if (!getSub<3>()) {
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
        // and face, but only if our face actually exists.  If it does not
        // exist, this means that the face cannot be one of the two subspaces 
        // we're circulating around, and must be our forward subspace.

        if (forwardDimension == 2 && !forwardSub) {
            assert(forwardDimension == 2);
            // This can thus only happen if we're circulating around a simplex 
            // edge between a cell sub and a corner that is either on its edge
            // or on its face.  In either case, what we actually want to do
            // is swap the cell; because our face does not exist, swapping
            // the cell across it will not actually change the identity of our
            // cell or its sub, but will change which dimension is forward,
            // so that the next call to advance() will swap it to the
            // next face around the same edge, and then the following one will
            // move on to the next edge.
            auto relativeFaceAxis = (faceAxis == Axis::Q(edgeAxis) ? 1 : 2);
            relativeCellPos ^= relativeFaceAxis;
            if (!getSub<3>()) {
                assert(false);
                hitError = true;
                return;
            }
        }
        else {
            auto newEdgeAxis = Axis::Axis(7 & ~faceAxis & ~edgeAxis);
            auto newEdgeRelativePos = (faceAxis == Axis::Q(edgeAxis) ? 2 : 1);
            bool stillAdjacentToOldEdgeCell(
                (edgeRelativePos ^ relativeCellPos) & newEdgeRelativePos);
            const SimplexDCTree<3>* bestCell = nullptr;
            unsigned bestFloat;
            auto bestCellisDuplicate = false;
            auto absoluteEdgeDir = ((edgeRelativePos * 9) / (2 * edgeAxis)) & 7;
            auto cornerDirection = (edgeAxis * cornerPos) | absoluteEdgeDir;

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
                edgeRelativePos = newCornerDirection * (edgeAxis << 1);
                edgeRelativePos |= edgeRelativePos >> 3;
                edgeRelativePos &= 3;
            }
            else {
                // The new edge can't be adjacent to the old edge cell, as
                // if it were, the old edge cell would (due to having our
                // corner as a corner) necessarily be better than any cell
                // that does not have it as a corner (the implication of
                // being a duplicate).
                // And it  can't be in the direction of faceAxis, since if 
                // it were then to be the best cell all four candidates next 
                // to the edge would need to be duplicate in that direction, but
                // one of those two pairs has our face in between, and we
                // know that our face exists (as we branched on that condition).
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
    auto faceSub = getSub<2>();
    if (!faceSub) {
        return false;
    }
    auto cornerSub = getSub<0>();
    assert(cornerSub);
    if (cornerSub == faceSub) {
        return false;
    }
    auto edgeSub = getSub<1>();
    assert(edgeSub);
    if (edgeSub == cornerSub || edgeSub == faceSub) {
        return false;
    }
    auto cellSub = getSub<3>();
    assert(cellSub);
    if (cellSub == cornerSub || cellSub == edgeSub || cellSub == faceSub) {
        return false;
    }
    return true;
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
    assert(simplexCell);
    auto neighborCell = edgeCell->neighbor(
        NeighborIndex::fromPosAndFloating(
            posToSimplexAndNeighbor, floatToNeighborCell));
    if (simplexCell == neighborCell) {
        // The face does not actually exist.
        return nullptr;
    }
    auto neighborIsHigher = 
        (posToSimplexAndNeighbor ^ floatToNeighborCell) & faceAxis;
    auto faceFromHigher = NeighborIndex::fromPosAndFloating(0, faceAxis);
    auto faceFromLower = NeighborIndex::fromPosAndFloating(7, faceAxis);
    auto faceFromSimplexCell = 
        neighborIsHigher ? faceFromLower : faceFromHigher;
    auto faceFromNeighbor = neighborIsHigher ? faceFromHigher : faceFromLower;
    if (neighborCell && neighborCell->leafLevel() < simplexCell->leafLevel()) {
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
    assert(simplexCell);
    return simplexCell->leaf->collapsedSub(26);
}

}   // namespace libfive
