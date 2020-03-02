/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "libfive/render/brep/simplex_dc/simplex_dc_edges.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "libfive/render/brep/simplex/simplex_tree.hpp"
#include "libfive/render/brep/object_pool.hpp"
#include "libfive/render/brep/indexes.hpp"

namespace libfive {

template <unsigned N>
SimplexDCEdges<N>::SimplexDCEdges(PerThreadOutput& m, Pool&pool)
    : parent_pool(pool)
{
    // Nothing to do here
}

template<unsigned N>
inline SimplexDCEdges<N>::~SimplexDCEdges()
{
    parent_pool.claim(my_pool);
}

template <unsigned N>
template <Axis::Axis A>
void SimplexDCEdges<N>::load(const std::array<Input*, 1 << (N - 1)>& ts)
{
    // First, get the minimum leaf level.
    auto minIdx = std::min_element(ts.begin(), ts.end(),
        [](const Input* a, const Input* b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();
    auto minLevel = ts[minIdx]->leafLevel();
    auto ptr = my_pool.get();
    ptr->lowPt = ts[minIdx]->region.lower[Axis::toIndex(A)];
    auto& count = ptr->refcount;
    std::array<bool, 1 << (N - 1)> toUse;
    toUse.fill(true);
    if constexpr (N == 3) {
        auto testAgainst = 3 & ~minIdx;
        for (auto i = 1; i < 4; i = i << 1) {
            auto test = minIdx ^ i;
            if (ts[test] == ts[testAgainst]) {
                toUse[test] = false;
                toUse[testAgainst] = false;
                break; // There can't be 3 matching.
            }
        }
    }
    for (auto index = 0; index < ts.size(); ++index) {
        if (!toUse[index]) {
            continue;
        }
        auto t = ts[index];
        if (t->type == Interval::UNKNOWN) {
            continue;
        }
        auto& edge = t->leaf->edgeFromReduced(A, ts.size() - 1 - index);
        if (t->leafLevel() == minLevel) {
            assert(std::holds_alternative<SimplexDCMinEdge<N>::EdgeVec>(edge));
            assert(std::get<SimplexDCMinEdge<N>::EdgeVec>(edge).empty());
            edge.emplace<SimplexDCMinEdge<N>*>(ptr);
            ++count;
        }
        else {
            assert(std::holds_alternative<SimplexDCMinEdge<N>::EdgeVec>(edge));
            std::get<SimplexDCMinEdge<N>::EdgeVec>(edge).push_back(ptr);
            // Do not increment count, as the EdgeVec is non-owning.
        }
    }
}

}   // namespace libfive
