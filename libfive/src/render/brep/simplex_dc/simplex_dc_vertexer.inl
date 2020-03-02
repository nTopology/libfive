/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "libfive/render/brep/simplex_dc/simplex_dc_vertexer.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "libfive/render/brep/simplex_dc/simplex_qef.hpp"
#include "libfive/render/brep/per_thread_brep.hpp"
#include "libfive/render/brep/settings.hpp"

namespace libfive {

template <unsigned N>
template <Axis::Axis A>
void SimplexDCVertexer<N>::load(const std::array<Input*, 1 << (N - 1)> & ts)
{
    const auto index = std::min_element(ts.begin(), ts.end(),
        [](const Input* a, const Input* b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();

    auto& edgeVariant = ts[index]->leaf->
                        edgeFromReduced(A, ts.size() - 1 - index);
    assert(std::holds_alternative<SimplexDCMinEdge<N>*>(edgeVariant));

    auto& edge = *std::get<SimplexDCMinEdge<N>*>(edgeVariant);

    // Calculate the subspace for the edge.
    std::array<int, 3> reductions{ 0, 0, 0 };
    reductions[(Axis::toIndex(A) + 1) % 3] = 1 + (index & 1);
    reductions[(Axis::toIndex(A) + 2) % 3] = 1 + bool(index & 2);

    auto totalReductions = reductions[0] +
        3 * reductions[1] +
        9 * bool(N == 3) * reductions[2];

    auto edgeSubspaceIndex = (ipow(3, N) - 1) - totalReductions;

    auto edgeSub = ts[index]->leaf->collapsedSub(edgeSubspaceIndex);

    auto lowCornerSubIndex = edgeSubspaceIndex - 2 * ipow(3, Axis::toIndex(A));
    auto highCornerSubIndex = edgeSubspaceIndex - ipow(3, Axis::toIndex(A));

    std::array<const SimplexLeafSubspace<N>*, 2> cornerSubs = {
        ts[index]->leaf->collapsedSub(lowCornerSubIndex),
        ts[index]->leaf->collapsedSub(highCornerSubIndex)
    };

    // In 3d, faceSubs follows the order axis=Q(A)/low R(A), 
    // axis=Q(A)/high R(A), axis=R(A)/low Q(A), axis=R(A)/high Q(A).
    // In 2d, of course, it is ignored entirely (though sometimes
    // copied to another variable that is then ignored; that should
    // get optimized out).
    std::array<const SimplexLeafSubspace<N>*, 4> faceSubs;

    if constexpr (N == 3) {
        for (auto i = 0; i < 4; ++i) {
            auto cellIdxA = (i & 1) ? 3 : 0;
            auto cellIdxB = cellIdxA ^ (i & 2 ? 2 : 1);
            if (ts[cellIdxA] == ts[cellIdxB]) {
                // The face is interior to a merged cell and hence
                // does not exist.
                faceSubs[i] = nullptr;
                continue;
            }
            auto bestCellIdx = std::min(cellIdxA, cellIdxB,
                [&ts](const int a, const int b)
            { return ts[a]->leafLevel() < ts[b]->leafLevel(); });
            if (ts[bestCellIdx]->type == Interval::UNKNOWN) {
                continue;
            }
            auto faceAxis = (i & 2) ? R(A) : Q(A);
            bool isUpperToFace(bestCellIdx & (i & 2 ? 2 : 1));
            auto faceSubspaceIndex = 
                26 - ipow(3, Axis::toIndex(faceAxis)) * (isUpperToFace + 1);
            faceSubs[i] = ts[bestCellIdx]->leaf->collapsedSub(faceSubspaceIndex);
        }
    }

    for (auto cellIdx = 0; cellIdx < ts.size(); ++cellIdx) {
        if (ts[cellIdx]->type == Interval::UNKNOWN) {
            continue;
        }
        auto cellSub = ts[cellIdx]->leaf->collapsedSub(ipow(3, N) - 1);
        if (cellSub == edgeSub) {
            // Degenerate simplices do not get vertices.
            continue;
        }
        for (auto faceAxisIdx = 0; faceAxisIdx < 2; ++faceAxisIdx) {
            auto faceAxis = faceAxisIdx ? R(A) : Q(A);
            if (N == 2 && faceAxis == Axis::Z) {
                continue;
            }
            auto faceSubIdx = 
                faceAxisIdx << 1 | int((cellIdx & (2 >> faceAxisIdx)) != 0);
            auto faceSub = faceSubs[faceSubIdx];
            if (N == 3 && faceSub == nullptr) {
                // Face does not exist, so neither does the simplex we'd be
                // handling.
                continue;
            }
            if (N == 3 && (faceSub == edgeSub || faceSub == cellSub)) {
                // Degenerate simplices do not get vertices.
                continue;
            }
            for (auto corner = 0; corner < 2; ++corner) {
                auto cornerSub = cornerSubs[corner];
                if (cornerSub == edgeSub || cornerSub == cellSub ||
                    (N == 3 && cornerSub == faceSub)) {
                    // Degenerate simplices do not get vertices.
                    continue;
                }
                auto& simplex = edge.simplexWithEdgeReducedCell(
                    A, faceAxis, cellIdx, corner);
                assert(simplex.intersectionCount() == 3 ||
                    simplex.intersectionCount() == 4 ||
                    simplex.intersectionCount() == 0);
                auto getVerts = [&]()->SubspaceVertArray {
                    if constexpr (N == 3) {
                        return { cornerSub, edgeSub, faceSub, cellSub };
                    }
                    else {
                        return { cornerSub, edgeSub, cellSub };
                    }
                };
                calcAndStoreVert(simplex, getVerts());
            }
        }
    }
}

template<unsigned N>
void SimplexDCVertexer<N>::calcAndStoreVert(
    DCSimplex<N>& simplex, SubspaceVertArray vertsFromSubspaces)
{
    Eigen::Matrix<double, N, N + 1> vertices;
    auto state = 0; // Bit 1 for inside, bit 2 for empty.
    for (auto i = 0; i <= N; ++i) {
        if (vertsFromSubspaces[i]->inside) {
            state |= 1;
        }
        else {
            state |= 2;
        }
        vertices.col(i) = vertsFromSubspaces[i]->vert;
    }
    if (state != 3) {
        assert(state == 1 || state == 2);
        assert(simplex.intersectionCount() == 0);
        return; // No point in making a vertex if the simplex is all
                // filled or all empty.
    }

    SimplexQEF<N> qef(std::move(vertices), 
                      settings.simplex_dc_padding_rate, 
                      settings.simplex_bounding_eigenvalue_cutoff);

    for (auto a = 0; a < N; ++a) {
        for (auto b = a + 1; b <= N; ++b) {
            auto intersection = simplex.intersection(a, b);
            assert(intersection != &DCSimplex<N>::dupVertIntersection);
            if (intersection != nullptr) {
                qef += *intersection;
            }
        }
    }
    simplex.vert = qef.solve();
    assert(simplex.index == 0);
    simplex.index = m.pushVertex(simplex.vert) + offset;
}


}   // namespace libfive
