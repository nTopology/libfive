/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "libfive/render/brep/simplex_dc/simplex_dc_collapser.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "libfive/render/brep/indexes.hpp"

namespace libfive {

template <unsigned N, unsigned K>
SimplexDCCollapser<N, K>::SimplexDCCollapser(
    PerThreadOutput& m, Tree t)
    : eval(new Evaluator(t)), owned(true)
{
    // Nothing to do here
}

template <unsigned N, unsigned K>
SimplexDCCollapser<N, K>::SimplexDCCollapser(
    PerThreadOutput& m, Evaluator* es)
    : eval(es), owned(false)
{
    // Nothing to do here
}

template <unsigned N, unsigned K>
SimplexDCCollapser<N, K>::~SimplexDCCollapser() {
    if (owned) {
        delete eval;
    }
}

template <unsigned N, unsigned K>
template <Axis::Axis A>
void SimplexDCCollapser<N, K>::load(const std::array<Input*, 1 << (N - 2)>& ts)
{
    static_assert(K == 2); // For <3, 3>, it should have been specialized (as
                           // trivial) in the header.
    static_assert(N == 3 || A == Axis::Z); // Likewise, in the .cpp.
    auto index = 0;
    if (N == 3 && ts[1]->leafLevel() < ts[0]->leafLevel()) {
        index = 1;
    }
    auto cell = ts[index];
    constexpr auto AInd = Axis::toIndex(A);
    auto subOffset = 1 + bool(N == 2) + index;
    auto faceSubIndex = 26 - (ipow(3, AInd)) * subOffset;
    NeighborIndex faceNeighbor(faceSubIndex);
    auto& subs = cell->leaf->sub;
    auto faceSub = subs[faceSubIndex].load();
    const auto& faceSubVert = faceSub->vert;
    const auto& region = cell->region;
    if (N == 3) {
        assert(faceSubVert(AInd) == index ? region.lower(AInd) 
                                          : region.upper(AInd));
    }
    std::array<int, 2> posArray;
    for (auto i = 0; i < 2; ++i) {
        auto iShift = (i + AInd + 1) % 3;
        assert(faceSubVert(iShift) >= region.lower(iShift));
        assert(faceSubVert(iShift) <= region.upper(iShift));
        assert(region.lower(iShift) < region.upper(iShift));
        if (faceSubVert(iShift) == region.lower(iShift)) {
            posArray[i] = 0;
        }
        else if (faceSubVert(iShift) == region.upper(iShift)) {
            posArray[i] = 1;
        }
        else {
            posArray[i] = 2;
        }
    }
    auto pos = posArray[0] + 3 * posArray[1];
    assert(pos <= 8);
    if (pos != 8) {
        NeighborIndex relativeNeighbor(pos);
        assert(relativeNeighbor.dimension() < 2);
        // In several situations, we want to "test" an edge.  This returns 0 or 
        // 1 if that edge will collapse to a corner, 2 if it will not collapse,
        // and 3 if it is not minimal.
        auto testEdge = 
            [&cell, &faceNeighbor, &subs, &region](NeighborIndex relativeEdge) {
            auto absoluteEdge =
                faceNeighbor.fromRelativeToThis(relativeEdge);
            // First, check if the edge is in fact minimal; if not, we cannot
            // collapse to it.
            auto floating = absoluteEdge.floating();
            auto pos = absoluteEdge.pos();
            for (auto toFloat = 0; toFloat < (1 << N); ++toFloat) {
                if (floating & ~toFloat) {
                    continue; // If the edge has it floating, so must 
                              // the neighbor direction.
                }
                auto neighbor = cell->neighbor(
                    NeighborIndex::fromPosAndFloating(pos, toFloat));
                if (neighbor && neighbor->isBranch()) {
                    // The edge for our cell is split.  (This must be
                    // from a cell not by our face, as otherwise "cell" would
                    // not be the lowest-level neighbor of the face.
                    assert(toFloat != faceNeighbor.floating());
                    return 3;
                }
            }

            // Then, check if it in turn will collapse to a corner.
            auto edgeAxis = highestbit(absoluteEdge.floating());
            auto edgeSub = 
                subs[faceNeighbor.fromRelativeToThis(relativeEdge).i].load();
            const auto& edgeVert = edgeSub->vert;
            assert(edgeVert(edgeAxis) >= region.lower(edgeAxis));
            assert(edgeVert(edgeAxis) <= region.upper(edgeAxis));
            if (edgeVert(edgeAxis) == region.lower(edgeAxis)) {
                return 0;
            }
            else if (edgeVert(edgeAxis) == region.upper(edgeAxis)) {
                return 1;
            }
            else {
                return 2;
            }
        };
        auto setToCenter = [&]() {
            for (auto i = 0; i < 3; ++i) {
                if (i != AInd) {
                    faceSub->vert(i) = region.center()(i);
                }
            }
            switch (cell->type) {
            case Interval::AMBIGUOUS:
            {
                std::array<bool, ipow(3, N)> alreadySolved;
                alreadySolved.fill(true);
                alreadySolved[faceSubIndex] = false;
                cell->saveVertexSigns(eval, eval->getDeck()->tape, alreadySolved);
                break;
            }
            case Interval::FILLED:
            case Interval::EMPTY:
                faceSub->inside = (cell->type == Interval::FILLED);
                break;
            default:
                assert(false);
                faceSub->inside = false;
            }
            pos = 8;
        };
        if (relativeNeighbor.dimension() == 1) {
            auto testResult = testEdge(relativeNeighbor);
            faceSub->faceEdgeCheckResult[0] = testResult;
            faceSub->faceEdgeCheckResult[1] = -1;
            if (testResult == 3) {
                setToCenter();
            }
            else {
                pos = relativeNeighbor.fromRelativeToThis(testResult).i;
            }
            if (pos != relativeNeighbor.i) {
                relativeNeighbor = pos;
            }
        }
        // If we're now at a corner (whether originally or from collapsing an
        // edge) we need to look at intermediate subspaces (i.e. edges).
        if (relativeNeighbor.dimension() == 0) 
        {
            std::array<NeighborIndex, 2> intermediates; // Relative to the face.
            std::array<int, 2> intermediatesResults;
            auto cornerPos = relativeNeighbor.pos();
            assert(relativeNeighbor.floating() == 0);
            for (auto i = 0; i < 2; ++i) {
                intermediates[i] = 
                    NeighborIndex::fromPosAndFloating(cornerPos, 1 << i);
                intermediatesResults[i] = testEdge(intermediates[i]);
            }
            faceSub->faceEdgeCheckResult = intermediatesResults;
            // If there is an intermediate space that does not collapse,
            // we want to collapse to it.  If both, compare them.
            auto collapseTarget = -1;
            if (intermediatesResults[0] == 2) {
                if (intermediatesResults[1] == 2) {
                    auto absoluteA = faceNeighbor.fromRelativeToThis(
                        intermediates[0]);
                    auto absoluteB = faceNeighbor.fromRelativeToThis(
                        intermediates[1]);
                    if (compareCollapseCandidates(
                        cell, *faceSub, absoluteA.i, absoluteB.i)) {
                        collapseTarget = 0;
                    }
                    else {
                        collapseTarget = 1;
                    }
                }
                else {
                    collapseTarget = 0;
                }
            }
            else if (intermediatesResults[1] == 2) {
                collapseTarget = 1;
            }
            else {
                // We won't collapse to an edge; we will either use this
                // corner, or set to the center.
                for (auto i = 0; i < 2; ++i) {
                    if (intermediatesResults[i] == 3 ||
                        intermediates[i].fromRelativeToThis(
                            intermediatesResults[i]).i != relativeNeighbor.i) {
                        setToCenter();
                        break;
                    }
                }
            }
            if (collapseTarget >= 0) {
                pos = intermediates[collapseTarget].i;
            }
        }
    }
    faceSub->collapseRef = pos;
}

template <unsigned N, unsigned K>
template <Axis::Axis A>
void SimplexDCCollapser<N, K>::load(const std::array<Input*, 1 << (N - 1)>& ts)
{
    auto index = std::min_element(ts.begin(), ts.end(),
        [](const Input* a, const Input* b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();
    auto cell = ts[index];
    constexpr auto AInd = Axis::toIndex(A);
    auto subCorner = ts.size() - 1 - index;
    subCorner *= (2 * A);
    assert(!(subCorner & subCorner >> N));
    subCorner |= (subCorner >> N);
    subCorner &= ((1 << N) - 1);
    auto edgeNeighbor = NeighborIndex::fromPosAndFloating(subCorner, A);
    auto& subs = cell->leaf->sub;
    auto edgeSub = subs[edgeNeighbor.i].load();
    const auto& edgeVert = edgeSub->vert;
    const auto& region = cell->region;

    assert(edgeVert(AInd) >= region.lower(AInd));
    assert(edgeVert(AInd) <= region.upper(AInd));
    if (edgeVert(AInd) == region.lower(AInd)) {
        edgeSub->collapseRef = 0;
    }
    else if (edgeVert(AInd) == region.upper(AInd)) {
        edgeSub->collapseRef = 1;
    }
    else {
        edgeSub->collapseRef = 2;
    }
}

template<unsigned N, unsigned K>
bool SimplexDCCollapser<N, K>::compareCollapseCandidates(
    const Input* cell, const SimplexLeafSubspace<N>& toCollapse, 
    unsigned candidateA, unsigned candidateB)
{
    std::array<NeighborIndex, 2> neighbors{ candidateA, candidateB };
    std::array<int, 2> containing{ 0, 0 };
    std::array<int, 2> floating;
    std::array<int, 2> dimensions;
    const auto& vert = toCollapse.vert;
    const auto& region = cell->region;
    std::array<int, N> vertSub;
    for (auto axis = 0; axis < N; ++axis) {
        assert(region.lower(axis) < region.upper(axis));
        assert(vert(axis) >= region.lower(axis));
        assert(vert(axis) <= region.upper(axis));
        if (vert(axis) == region.lower(axis)) {
            vertSub[axis] = 0;
        }
        else if (vert(axis) == region.upper(axis)) {
            vertSub[axis] = 1;
        }
        else {
            vertSub[axis] = 2;
        }
    }

    for (auto cand = 0; cand < 2; ++cand) {
        auto pos = neighbors[cand].pos();
        floating[cand] = neighbors[cand].floating();
        for (auto axis = 0; axis < N; ++axis) {
            if (floating[cand] & (1 << axis)) {
                // Our first criterion is in how many axes our candidate
                // contains vert; if the candidate is floating, it does contain
                // vert.
                ++containing[cand];
            }
            else {
                bool high(pos & (1 << axis));
                if (int(high) == vertSub[axis]) {
                    // This means that the candidate contains vert even without
                    // being floating there.
                    ++containing[cand];
                }
            }
        }
    }
    assert(floating[0] != floating[1]); // Since they share a corner 
                                        // but are not the same.
    if (containing[0] != containing[1]) {
        return containing[0] > containing[1];
    }
    // If both candidates contain vert in the same number of dimensions, we want
    // to use the one in which fewer of those dimensions are floating (as 
    // containing vert in a fixed dimension is a better match).
    auto dim0 = bitcount(floating[0]);
    auto dim1 = bitcount(floating[1]);
    if (dim0 != dim1) {
        return dim0 < dim1;
    }
    // If we're still tied, compare the candidates' vertices.  First, see if
    // one is a better match in terms of inside/outside.
    assert(cell->leaf);
    auto subA = cell->leaf->sub[candidateA].load();
    auto subB = cell->leaf->sub[candidateB].load();

    if (subA->inside != subB->inside) {
        return subA->inside == toCollapse.inside;
    }

    // Otherwise, check distance to the original vertex.

    auto distA = (vert - subA->vert).squaredNorm();
    auto distB = (vert - subB->vert).squaredNorm();

    if (distA != distB) {
        return distA < distB;
    }

    // Otherwise, as a final tiebreaker, just use which axes are floating.
    return floating[0] < floating[1];

}

}   // namespace libfive
