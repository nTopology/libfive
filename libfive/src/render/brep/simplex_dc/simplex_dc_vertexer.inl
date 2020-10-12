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

template <unsigned N, bool indexing>
template <Axis::Axis A>
void SimplexDCVertexer<N, indexing>::load(
    const std::array<Input*, 1 << (N - 1)> & ts)
{
    const auto index = std::min_element(ts.begin(), ts.end(),
        [](const Input* a, const Input* b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();

    auto& edge = ts[index]->leaf->edgeFromReduced(A, ts.size() - 1 - index);

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
            faceSubs[i] = 
                ts[bestCellIdx]->leaf->collapsedSub(faceSubspaceIndex);
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
                auto& simplex = edge->simplexWithEdgeReducedCell(
                    A, faceAxis, cellIdx, corner);
                assert(simplex.intersectionCount() == 3 ||
                    simplex.intersectionCount() == 4 ||
                    simplex.intersectionCount() == 0);
                // We actually only need getVerts() for the non-indexing version
                // and only need orientation for the indexing version, but 
                // they'll likely be optimized away anyway in the versions where
                // they're not used.
                auto getVerts = [&]()->SubspaceVertArray {
                    if constexpr (N == 3) {
                        return { cornerSub, edgeSub, faceSub, cellSub };
                    }
                    else {
                        return { cornerSub, edgeSub, cellSub };
                    }
                };
                auto directionOrientation =
                    (corner == 1) == (cellIdx == 0 || cellIdx == 3) == (N == 2);
                auto axisOrientation = (N == 2) || (faceAxisIdx == 1);
                calcAndStoreVert(simplex, getVerts(), 
                                 directionOrientation == axisOrientation);
            }
        }
    }
}

template<unsigned N, bool indexing>
void SimplexDCVertexer<N, indexing>::calcAndStoreVert(
    DCSimplex<N>& simplex, 
    SubspaceVertArray vertsFromSubspaces, 
    bool orientation)
{
    Eigen::Matrix<double, N, N + 1> vertices;
    auto state = 0; // Bit 1 << N for bit N being inside.
    for (auto i = 0; i <= N; ++i) {
        if (vertsFromSubspaces[i]->inside) {
            state |= (1 << i);
        }
        vertices.col(i) = vertsFromSubspaces[i]->vert;
    }
    if constexpr (!indexing) {
        const auto& settings = data;
        if (state == 0 || state == (1 << (N + 1)) - 1) {
            assert(simplifiedState == 1 || simplifiedState == 2);
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
        auto vertex = qef.solve();
        auto AtANormalized = Eigen::Matrix<double, N, N>::Zero().eval();
        for (const auto& intersection : simplex.intersections) {
            if (DCSimplex<N>::isValid(intersection.load())) {
                AtANormalized += intersection.load()->AtANormalized;
            }
        }
        auto rank = 0;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, N, N>>
            es(AtANormalized);
        auto eigenvalues = es.eigenvalues().real();
        for (unsigned j = 0; j < N; ++j) {
            rank += (fabs(eigenvalues[j]) >= EIGENVALUE_CUTOFF);
        }
        assert(simplex.collapseTarget == nullptr);
        assert(simplex.hasVert == false);
        simplex.hasVert = true;
        for (const auto& intersection : simplex.intersections) {
            if (DCSimplex<N>::isValid(intersection.load()) &&
                intersection.load()->get_rank() >= rank) {
                if (simplex.collapseTarget) {
                    simplex.collapseTarget = nullptr;
                    break;
                }
                else {
                    simplex.hasVert = false;
                    simplex.collapseTarget = intersection.load();
                }
            }
        }
        if (simplex.hasVert) {
            simplex.vert = vertex;
        }
    }
    else {
        // Indexing version.
        if (simplex.collapseTarget) {
            while (simplex.collapseTarget->collapseTarget != nullptr) {
                simplex.collapseTarget = simplex.collapseTarget->collapseTarget;
            }
            simplex.index.store(simplex.collapseTarget->index);
            simplex.vert = simplex.collapseTarget->
                normalized_mass_point().template head<N>();
        }
        else if (simplex.hasVert) {
            auto offset = data;
            assert(simplex.index == 0);
            simplex.index = m.pushVertex(simplex.vert) + offset;
        }
        else if constexpr (N == 3) {
            assert(simplex.index == 0);
            switch (state) {
            case 7:
            case 11:
            case 13:
            case 14:
                // This is the same as its complement state with the opposite
                // orientation.
                orientation = !orientation;
                state = 15 - state;
                // fallthrough
            case 1:
            case 2:
            case 4:
            case 8:
                // There is exactly one corner different than the others,
                // so we will be making one triangle.
            {
                std::array<std::pair<int, int>, 3> indices;
                switch (state) {
                case 1:
                    indices = { std::pair{0, 1}, {0, 2}, {0, 3} };
                    break;
                case 2:
                    indices = { std::pair{1, 0}, {1, 3}, {1, 2} };
                    break;
                case 4:
                    indices = { std::pair{2, 0}, {2, 1}, {2, 3} };
                    break;
                case 8:
                    indices = { std::pair{3, 0}, {3, 2}, {3, 1} };
                    break;
                default:
                    assert(false);
                    return;
                }
                if (!orientation) {
                    std::swap(indices[1], indices[2]);
                }
                Eigen::Matrix<uint32_t, 3, 1> triangle;
                for (auto i = 0; i < 3; ++i) {
                    auto intersection = simplex.intersection(
                        indices[i].first, indices[i].second);
                    assert(DCSimplex<N>::isValid(intersection));
                    assert(intersection->index != 0);
                    triangle[i] = intersection->index;
                }
                m.branes.push_back(triangle);
            }
            break;
            case 9:
            case 10:
            case 12:
                // Again, reduce cases by flipping orientation.
                orientation = !orientation;
                state = 15 - state;
                // fallthrough
            case 3:
            case 5:
            case 6:
                // The corner states split 2 and 2, so we will be making two
                // different triangles. 
            {
                std::array<std::pair<int, int>, 4> indices;
                switch (state) {
                case 3:
                    indices = { std::pair{0, 2}, {0, 3}, {1, 3}, {1, 2} };
                    break;
                case 5:
                    indices = { std::pair{0, 1}, {2, 1}, {2, 3}, {0, 3} };
                    break;
                case 6:
                    indices = { std::pair{0, 1}, {3, 1}, {3, 2}, {0, 2} };
                    break;
                default:
                    assert(false);
                    return;
                }
                if (!orientation) {
                    std::swap(indices[1], indices[3]);
                }
                std::array<const SimplexDCIntersection<N>*, 4> intersections;
                for (auto i = 0; i < 4; ++i) {
                    auto intersection = simplex.intersection(
                        indices[i].first, indices[i].second);
                    assert(DCSimplex<N>::isValid(intersection));
                    assert(intersection->index != 0);
                    intersections[i] = intersection;
                }
                // Now to determine whether to do {0, 1, 3} and {2, 3, 1}, or
                // {1, 2, 0} and {3, 0, 2}
                auto double1And3 = false;
                std::array<int, 2> ranks02{ intersections[0]->get_rank(),
                                            intersections[2]->get_rank() };
                // Sort it in reverse order, so that the higher rank is first.
                std::sort(ranks02.rbegin(), ranks02.rend());
                std::array<int, 2> ranks13{ intersections[1]->get_rank(),
                                            intersections[3]->get_rank() };
                std::sort(ranks13.rbegin(), ranks13.rend());
                if (ranks02 < ranks13) {
                    // We want to double the higher-rank pair (prioritizing the
                    // highest rank), to better match sharp features.
                    double1And3 = true;
                }
                else if (ranks02 == ranks13) {
                    // We have a tie, so go on to the next approach:
                    auto getPt = [&](int i) {
                        return intersections[i]->normalized_mass_point()
                            .template head<3>().eval();
                    };
                    auto distanceSq02 = (getPt(0) - getPt(2)).squaredNorm();
                    auto distanceSq13 = (getPt(1) - getPt(3)).squaredNorm();
                    // Use the approach that minimizes the distance between
                    // the two doubled points, as that should minimize the
                    // folding over due to the triangulation.
                    double1And3 = distanceSq13 < distanceSq02;
                }
                auto getIdx = [&](int index) {
                    return intersections[index]->index.load();
                };
                Eigen::Matrix<uint32_t, 3, 1> triangle;
                if (double1And3) {
                    triangle << getIdx(0), getIdx(1), getIdx(3);
                    m.branes.push_back(triangle);
                    triangle << getIdx(2), getIdx(3), getIdx(1);
                    m.branes.push_back(triangle);
                }
                else {
                    triangle << getIdx(1), getIdx(2), getIdx(0);
                    m.branes.push_back(triangle);
                    triangle << getIdx(3), getIdx(0), getIdx(2);
                    m.branes.push_back(triangle);
                }
            }
            break;
            default:
                assert(false);
            }
        }
        else {
            static_assert(N == 2);
            // TBD.
        }
    }
}

}   // namespace libfive
