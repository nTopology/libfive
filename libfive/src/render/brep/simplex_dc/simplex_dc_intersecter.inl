/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "libfive/render/brep/simplex_dc/simplex_dc_intersecter.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "libfive/render/brep/simplex/simplex_tree.hpp"
#include "libfive/render/brep/object_pool.hpp"
#include "libfive/render/brep/indexes.hpp"
#include "libfive/render/brep/per_thread_brep.hpp"
#include <set>

namespace libfive {

template <unsigned N>
SimplexDCIntersecter<N>::SimplexDCIntersecter(
    PerThreadOutput& m, Tree t, 
    ObjectPool<SimplexDCIntersection<N>>& pool, Perp perp)
    : parent_pool(pool), eval(new Evaluator(t)), owned(true), perp(perp), m(m)
{
    // Nothing to do here
}

template <unsigned N>
SimplexDCIntersecter<N>::SimplexDCIntersecter(
    PerThreadOutput& m, Evaluator* es, 
    ObjectPool<SimplexDCIntersection<N>>& pool, Perp perp)
    : parent_pool(pool), eval(es), owned(false), perp(perp), m(m)
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

        const auto& faceSub = *ts[index]->leaf->sub[faceSubspaceIndex].load();

        // Now that we have the intersections, we need to write them to every 
        // simplex that contains this face.  Again, we use the smallest of the 
        // associated trees.
        auto& edges = ts[index]->leaf->edges;

        for (auto cell = 0; cell < 2; ++cell) {
            if (ts[cell]->type != Interval::AMBIGUOUS) {
                continue;
            }
            const auto& cellSub = *ts[cell]->leaf->sub[ipow(3, N) - 1].load();
            if (faceSub.inside == cellSub.inside) {
                continue;
            }
            auto& inside = faceSub.inside ? faceSub.vert : cellSub.vert;
            auto& outside = faceSub.inside ? cellSub.vert : faceSub.vert;
            auto intersection = searchEdge(inside, outside,
                                           ts[cell]->leaf->tape);
            for (auto axisIdx = 0; axisIdx < 2; ++axisIdx) {
                auto edgeAxis = axisIdx == 0 ? Axis::Q(A) : Axis::R(A);
                for (auto edgePos = 0; edgePos < 2; ++edgePos) {
                    auto& edge = ts[index]->leaf->edgeFromFaceAndIndex(
                        edgeAxis, A, edgePos == 1, index == 0);
                    for (auto cornerIdx = 0; cornerIdx < 2; ++cornerIdx) {
                        // cornerIdx is the position with respect to edgeAxis,
                        // and edgePos is that with respect to the third axis
                        // (neither edge nor face), so how we combine the two
                        // into a single face-reduced corner depends on whether
                        // the edge axis is Q or R of the face axis (as bit 1
                        // of the combined corner represents Q of the face axis,
                        // and bit 2 represents R.)
                        auto corner = axisIdx == 0 ? 2 * edgePos + cornerIdx :
                                                     edgePos + 2 * cornerIdx;
                        auto writeToSimplex = [&](SimplexDCMinEdge<N>* edge) {
                            auto& simplex = edge->simplexWithFaceReducedCell(
                            edgeAxis, A, cell == 1, corner);
                            [[maybe_unused]] auto [res, success] =
                                simplex.insertIntersection(2, 3, intersection);
                            assert(success);
                            assert(res == intersection);
                        };
                        if (std::holds_alternative<SimplexDCMinEdge<N>*>(edge)) {
                            writeToSimplex(std::get<SimplexDCMinEdge<N>*>(edge));
                        }
                        else {
                            auto& edgeVec = std::get<SimplexDCMinEdge<N>::EdgeVec>(edge);
                            for (auto& minEdge : edgeVec) {
                                writeToSimplex(minEdge);
                            }
                        }
                    }
                }
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
    reductions[(Axis::toIndex(A) + 1) % 3] = 1 + (index & 1);
    reductions[(Axis::toIndex(A) + 2) % 3] = 1 + bool(index & 2);

    auto totalReductions = reductions[0] + 
                           3 * reductions[1] + 
                           9 * bool(N == 3) * reductions[2];

    auto edgeSubspaceIndex = (ipow(3, N) - 1) - totalReductions;

    const auto& edgeSub = *ts[index]->leaf->sub[edgeSubspaceIndex].load();

    auto& edge = ts[index]->leaf->edgeFromReduced(A, ts.size() - 1 - index);

    if constexpr (N == 3) {
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
                    // This can happen if the edge was on the border of the
                    // bounding box and the face is toward the outside, so
                    // the only non-Interval::UNKNOWN cell(s) don't border the face
                    // in question.
                    continue;
                }
                auto indexFromFace = bool(faceIndex & (1 << faceAxisIdx));
                auto faceSubIndex = (ipow(3, N) - 1) - 
                    ipow(3, Axis::toIndex(faceAxis)) * (indexFromFace + 1);
                const auto& faceSub = 
                    *ts[faceIndex]->leaf->sub[faceSubIndex].load();
                if (edgeSub.inside == faceSub.inside) {
                    continue;
                }
                auto& inside = edgeSub.inside ? edgeSub.vert : faceSub.vert;
                auto& outside = edgeSub.inside ? faceSub.vert : edgeSub.vert;
                auto intersection = searchEdge(inside, outside,
                    ts[faceIndex]->leaf->tape);
                for (auto cell = 0; cell < 2; ++cell) {
                    auto cellFrom4 = (faceAxisIdx ? facePosition + 2 * cell
                                                  : 2 * facePosition + cell);
                    for (auto corner = 0; corner < 2; ++corner) {
                        // Unlike the face loader, we don't need to handle the
                        // case in which there are multiple edges, since we're
                        // using the minimal cell adjoining the edge in 
                        // question (which is the edge we're loading).
                        auto& simplex = 
                            std::get<SimplexDCMinEdge<N>*>(edge)
                            ->simplexWithEdgeReducedCell(
                            A, faceAxis, cellFrom4, corner);
                        [[maybe_unused]] auto [res, success] =
                            simplex.insertIntersection(1, 2, intersection);
                        assert(success);
                        assert(res == intersection);
                    }
                }
            }
        }
    }

    SimplexDCIntersection<N>* doubledCellIntersection = nullptr;
    // Now handle edges between cell vertices and our edge vertex.
    for (auto cell = 0; cell < ts.size(); ++cell) {
        if (ts[cell]->type == Interval::UNKNOWN) {
            continue;
        }
        const auto& cellSub = *ts[cell]->leaf->sub[ipow(3, N) - 1].load();
        if (edgeSub.inside == cellSub.inside) {
            continue;
        }
        auto isDoubledOrOpposite = false;
        if constexpr (N == 3) {
            auto opposite = index ^ 3;
            if (ts[cell] == ts[opposite]) {
                assert(cell != index);
                isDoubledOrOpposite = true;
            }
        }
        auto& inside = edgeSub.inside ? edgeSub.vert : cellSub.vert;
        auto& outside = edgeSub.inside ? cellSub.vert : edgeSub.vert;
        auto intersection = [&](){
            if (isDoubledOrOpposite && doubledCellIntersection != nullptr) {
                return doubledCellIntersection;
            }
            auto out = searchEdge(inside, outside, ts[cell]->leaf->tape);
            if (isDoubledOrOpposite) {
                doubledCellIntersection = out;
            }
            return out;
        }();
        assert(intersection != nullptr);
        for (auto faceAxisIdx = 0; faceAxisIdx < 2; ++faceAxisIdx) {
            auto faceAxis = (faceAxisIdx == 0) ? Axis::Q(A) : Axis::R(A);
            if constexpr (N == 2) {
                if (faceAxis == Axis::Z) {
                    continue;
                }
            }
            else if (isDoubledOrOpposite) {
                auto acrossFace = cell ^ (1 << faceAxisIdx);
                if (ts[cell] == ts[acrossFace]) {
                    // We'd be looking at a "simplex" that uses a face that
                    // is internal to the merged cell and thus does not exist.
                    continue;
                }
            }
            for (auto corner = 0; corner < 2; ++corner) {
                // Unlike the face loader, we don't need to handle the
                // case in which there are multiple edges, since we're
                // using the minimal cell adjoining the edge in 
                // question (which is the edge we're loading).
                auto& simplex = std::get<SimplexDCMinEdge<N>*>(edge)
                    ->simplexWithEdgeReducedCell(A, faceAxis, cell, corner);
                [[maybe_unused]] auto [res, success] =
                    simplex.insertIntersection(1, N, intersection);
                assert(success);
                assert(res == intersection);
            }
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

    const auto& cornerSub = *ts[index]->leaf->sub[cornerSubspaceIndex].load();

    constexpr auto tsIndices = [](){
        static_assert(N == 2 || N == 3);
        if constexpr (N == 2) {
            return std::array{ 0, 1, 2, 3 };
        }
        else {
            return std::array{ 0, 1, 2, 3, 4, 5, 6, 7 };
        }
    }();

    std::array<SimplexDCMinEdge<N>*, N * 2> edges;

    // Handle edges between edge vertices and our corner vertex.
    for (auto edgeAxisIdx = 0; edgeAxisIdx < N; ++edgeAxisIdx) {
        auto edgeAxis = Axis::toAxis(edgeAxisIdx);
        for (auto edgePosition = 0; edgePosition < 2; ++edgePosition) {
            auto lowerLevel = [&](int indexA, int indexB) {
                if (bool(indexA & edgeAxis) != bool(edgePosition)) {
                    // That index cell does not border the edge we're 
                    // looking at.
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
            const auto& edgeVariant = 
                ts[edgeIndex]->leaf->edge(edgeAxis, lowerCorner);
            assert(std::holds_alternative<SimplexDCMinEdge<N>*>(edgeVariant));
            const auto& edge = std::get<SimplexDCMinEdge<N>*>(edgeVariant);
            edges[edgeAxisIdx * 2 + edgePosition] = edge;
            auto edgeSubIndex = CornerIndex(lowerCorner).neighbor().i + 
                2 * ipow(3, Axis::toIndex(edgeAxis));
            const auto& edgeSub = 
                *ts[edgeIndex]->leaf->sub[edgeSubIndex].load();
            if (cornerSub.inside == edgeSub.inside) {
                continue;
            }
            auto& inside = cornerSub.inside ? cornerSub.vert : edgeSub.vert;
            auto& outside = cornerSub.inside ? edgeSub.vert : cornerSub.vert;
            auto intersection = searchEdge(inside, outside,
                                           ts[edgeIndex]->leaf->tape);
            for (auto faceAxis : { Axis::Q(edgeAxis), Axis::R(edgeAxis) }) {
                if (N == 2 && faceAxis == Axis::Z) {
                    continue;
                }
                for (auto cell = 0; cell < (1 << (N - 1)); ++cell) {
                    if constexpr (N == 3) {
                        auto fullCell = edgeAxis * edgePosition |
                                        Axis::Q(edgeAxis) * (cell & 1) |
                                        Axis::R(edgeAxis) * (cell >> 1);
                        if (ts[fullCell] == ts[fullCell ^ faceAxis]) {
                            // The face that would be in this simplex
                            // does not actually exist.
                            continue;
                        }
                    }
                    auto& simplex = edge->simplexWithEdgeReducedCell(
                        edgeAxis, faceAxis, cell, !bool(edgePosition));
                    [[maybe_unused]] auto [res, success] =
                        simplex.insertIntersection(0, 1, intersection);
                    assert(success);
                    assert(res == intersection);
                }
            }
        }
    }

    if constexpr (N == 3) {
        // Handle edges between face vertices and our corner vertex.
        for (auto faceAxisIdx = 0; faceAxisIdx < N; ++faceAxisIdx) {
            auto faceAxis = Axis::toAxis(faceAxisIdx);
            std::array<SimplexDCIntersection<N>*, 4> intersectForDup;
            intersectForDup.fill(nullptr);
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
                auto faceSubIndex = (ipow(3, N) - 1) - 
                    ipow(3, Axis::toIndex(faceAxis)) * (int(useUpper) + 1);
                const auto& faceSub = 
                    *ts[faceIndex]->leaf->sub[faceSubIndex].load();
                if (cornerSub.inside == faceSub.inside) {
                    continue;
                }
                auto& inside = cornerSub.inside ? cornerSub.vert : faceSub.vert;
                auto& outside = cornerSub.inside ? faceSub.vert : cornerSub.vert;
                auto intersection = [&]() {
                    for (auto i : { 1, 2 }) {
                        auto& fromNeighbor = intersectForDup[facePosition ^ i];
                        if (fromNeighbor == nullptr) {
                            continue;
                        }
                        auto otherAxis = 
                            (i == 1) ? Axis::Q(faceAxis) : Axis::R(faceAxis);
                        if (ts[lowerFaceCell ^ otherAxis] == ts[lowerFaceCell] &&
                            ts[upperFaceCell ^ otherAxis] == ts[upperFaceCell]) {
                            intersectForDup[facePosition] = fromNeighbor;
                            return fromNeighbor;
                        }
                    }
                    auto out = searchEdge(inside, outside,
                                          ts[faceIndex]->leaf->tape);
                    intersectForDup[facePosition] = out;
                    return out;
                }();

                for (auto edgeAxisIdx = 0; edgeAxisIdx < 2; ++edgeAxisIdx) {
                    auto edgeAxis = (edgeAxisIdx == 0) ? Axis::Q(faceAxis) 
                                                       : Axis::R(faceAxis);
                    bool edgePosition(facePosition & (1 << edgeAxisIdx));
                    auto edge = 
                        edges[Axis::toIndex(edgeAxis) * 2 + edgePosition];
                    if (edge == nullptr) {
                        // This occurs if the corresponding edge is internal
                        // to a face or cell; in such case, we don't want
                        // to add any simplices using it.
                        continue;
                    }
                    for (auto cell = 0; cell < 2; ++cell) {
                        auto& simplex = edge->simplexWithFaceReducedCell(
                            edgeAxis, faceAxis, cell, 3 - facePosition);
                        [[maybe_unused]] auto [res, success] =
                            simplex.insertIntersection(0, 2, intersection);
                        assert(success);
                        assert(res == intersection);
                    }
                }
            }
        }
    }

    // Now handle edges between cell vertices and our corner vertex.
    std::array<SimplexDCIntersection<N>*, 1 << N> intersectForDup;
    intersectForDup.fill(nullptr);
    for (auto cell = 0; cell < ts.size(); ++cell) {
        if (ts[cell]->type == Interval::UNKNOWN) {
            continue;
        }
        const auto& cellSub = *ts[cell]->leaf->sub[ipow(3, N) - 1].load();
        if (cornerSub.inside == cellSub.inside) {
            continue;
        }
        auto& inside = cornerSub.inside ? cornerSub.vert : cellSub.vert;
        auto& outside = cornerSub.inside ? cellSub.vert : cornerSub.vert;
        auto intersection = [&]() {
            for (auto i = 0; i < N; ++i) {
                auto axis = Axis::toAxis(i);
                if (ts[cell] == ts[cell ^ axis] && 
                    intersectForDup[cell ^ axis] != nullptr) {
                    intersectForDup[cell] = intersectForDup[cell ^ axis];
                    return intersectForDup[cell];
                }
            }
            auto out = searchEdge(inside, outside, 
                                  ts[cell]->leaf->tape);
            intersectForDup[cell] = out;
            return out; 
        }();
        for (auto edgeAxisIdx = 0; edgeAxisIdx < N; ++edgeAxisIdx) {
            auto edgeAxis = Axis::toAxis(edgeAxisIdx);
            auto edgePosition = bool(cell & edgeAxis);
            auto edge = edges[edgeAxisIdx * 2 + edgePosition];
            if (edge == nullptr) {
                // All simplices from this iteration would use a 
                // nonexistent edge.
                continue;
            }
            for (auto faceAxisIdx = 0; faceAxisIdx < N - 1; ++faceAxisIdx) {
                auto faceAxis = (faceAxisIdx == 0) ? Axis::Q(edgeAxis) : 
                                                     Axis::R(edgeAxis);
                if (N == 2 && faceAxis == Axis::Z) {
                    continue;
                }
                if (N == 3 && ts[cell] == ts[cell ^ faceAxis]) {
                    // The simplex would use a nonexistent face.
                    continue;
                }
                auto& simplex = edge->simplex(edgeAxis, faceAxis, 
                                              cell, !edgePosition);
                [[maybe_unused]] auto [res, success] = 
                    simplex.insertIntersection(0, N, intersection);
                assert(success);
                assert(res == intersection);
            }
        }
    }
}

template<unsigned N>
SimplexDCIntersection<N>* SimplexDCIntersecter<N>::searchEdge(
    Eigen::Matrix<double, N, 1> inside, 
    Eigen::Matrix<double, N, 1> outside, 
    const std::shared_ptr<Tape>& tape)
{
    // This code is based on the simplex mesher, but without the assumption
    // that N == 3, and with the gradient calculated after finding the
    // intersection (as we need it to do the DC part of simplex DC).
    assert(tape.get() != nullptr);

    // There's an interesting question of precision + speed tradeoffs,
    // which mostly depend on how well evaluation scales in the
    // ArrayEaluator.  for now, we'll use the same value as XTree.
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
    targets[0].first << inside, perp;
    targets[0].second << outside, perp;

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
        if (!ambig(i))
        {
            intersection->push(pos.template head<N>(),
                               ds.col(i).template cast<double>()
                               .template head<N>(),
                               ds.col(i).w());
        }
        // Otherwise, we need to use the feature-finding special
        // case to find all possible derivatives at this point.
        else
        {
            const auto fs = eval->features(
                pos.template cast<float>(), tape);

            for (auto& f : fs)
            {
                intersection->push(pos.template head<N>(),
                                   f.template head<N>()
                                   .template cast<double>(),
                                   ds.col(i).w());
            }
        }
    }
    intersection->index = 
        m.pushVertex(intersection->normalized_mass_point().head<N>().eval());
    intersection->orientationChecker.setEndpts(outside, inside);
    return intersection;
}

}   // namespace libfive
