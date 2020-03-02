/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "libfive/render/brep/simplex_dc/simplex_dc_sub_consolidator.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "libfive/render/brep/indexes.hpp"

namespace libfive {


template <unsigned N>
template <Axis::Axis A>
void SimplexDCSubConsolidator<N>::load(
    const std::array<Input*, 1 << (N - 1)>& ts)
{
    // First, get the minimum leaf level.
    auto minIdx = std::min_element(ts.begin(), ts.end(),
        [](const Input* a, const Input* b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();
    auto minLevel = ts[minIdx]->leafLevel();
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
    auto getSubIdx = [&](int cell) {
        auto pos = ts.size() - 1 - cell;
        pos *= (A * 2);
        pos |= (pos >> N);
        pos &= (1 << N) - 1;
        return NeighborIndex::fromPosAndFloating(pos, A).i;
    };
    auto primary = ts[minIdx]->leaf->sub[getSubIdx(minIdx)].load();
    for (auto cell = 0; cell < ts.size(); ++cell) {
        if (!toUse[cell]) {
            continue;
        }
        if (cell == minIdx) {
            continue;
        }
        if (ts[cell]->leafLevel() != ts[minIdx]->leafLevel()) {
            assert(ts[cell]->leafLevel() > ts[minIdx]->leafLevel());
            continue;
        }
        ts[cell]->leaf->setSub(getSubIdx(cell), primary);
    }
}

template <unsigned N>
template <Axis::Axis A>
void SimplexDCSubConsolidator<N>::load(
    const std::array<Input*, 1 << (N - 2)>& ts)
{
    if constexpr (N == 3) {
        if (ts[0]->leafLevel() != ts[1]->leafLevel()) {
            return;
        }
        auto primary = ts[0]->leaf->sub[26 - ipow(3, Axis::toIndex(A))].load();
        ts[1]->leaf->setSub(26 - 2 * ipow(3, Axis::toIndex(A)), primary);
    }
    else {
        assert(false);
    }
}

template <unsigned N>
void SimplexDCSubConsolidator<N>::load(
    const std::array<Input*, 1 << N>& ts)
{
    // First, get the minimum leaf level.
    auto minIdx = std::min_element(ts.begin(), ts.end(),
        [](const Input* a, const Input* b)
    { return a->leafLevel() < b->leafLevel(); }) - ts.begin();
    std::array<bool, 1 << N> toUse;
    toUse.fill(true);
    for (auto pos = 0; pos < ts.size(); ++pos) {
        for (auto A = 1; A < ts.size(); A <<= 1) {
            if (pos & A) {
                continue;
            }
            else if (ts[pos] == ts[pos | A]) {
                toUse[pos] = false;
                toUse[pos | A] = false;
                if (N == 2) {
                    break; // Max cells sharing a corner in 2d is 2, as if
                           // there were 4 the corner would be interior 
                           // to a cell and this would not be called.
                }
            }
        }
    }
    auto getSubIdx = [&](int cell) {
        auto pos = ts.size() - 1 - cell;
        return NeighborIndex::fromPosAndFloating(pos, 0).i;
    };
    auto primary = ts[minIdx]->leaf->sub[getSubIdx(minIdx)].load();
    for (auto cell = 0; cell < ts.size(); ++cell) {
        if (!toUse[cell]) {
            continue;
        }
        if (cell == minIdx) {
            continue;
        }
        // However, if toUse is true, we want to apply it to this cell
        // even if it does have a higher level; this is the situation
        // where two-differently-sized cells share a corner, and those
        // do need to be identified with each other.
        ts[cell]->leaf->setSub(getSubIdx(cell), primary);
    }
}

}   // namespace libfive
