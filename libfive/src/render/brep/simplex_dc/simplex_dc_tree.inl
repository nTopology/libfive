/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"

namespace libfive {

////////////////////////////////////////////////////////////////////////////////

template<unsigned N>
SimplexDCIntersection<N>::SimplexDCIntersection()
{
    /* No need to reset the Intersection; that was done by its default 
     * constructor. */
}

template<unsigned N>
void SimplexDCIntersection<N>::reset()
{
    assert(refcount == 0);
    refcount = 0;
    orientationChecker.reset();
    Intersection::reset();
}

template<unsigned N>
const SimplexDCIntersection<N>* DCSimplex<N>::intersection(unsigned a, 
                                                           unsigned b) const
{
    if (a == b) {
        assert(false);
        return nullptr;
    }
    if (a > b) {
        std::swap(a, b);
    }
    if (a > N) {
        assert(false);
        return nullptr;
    }
    return intersections[b * (b - 1) / 2 + a];
}

template<unsigned N>
inline std::pair<const SimplexDCIntersection<N>*, bool> 
DCSimplex<N>::insertIntersection(unsigned a, 
                                 unsigned b, 
                                 SimplexDCIntersection<N>* intersection)
{
    if (a == b) {
        assert(false);
        return { nullptr, false };
    }
    if (a > b) {
        std::swap(a, b);
    }
    if (a > N) {
        assert(false);
        return { nullptr, false };
    }
    auto& target = intersections[b * (b - 1) / 2 + a];
    SimplexDCIntersection<N>* ptr(nullptr);
    if (target.compare_exchange_strong(ptr, intersection)) {
        ++intersection->refcount;
        return { intersection, true };
    }
    else {
        return { ptr, false };
    }
}

template<unsigned N>
SimplexDCMinEdge<N>::SimplexDCMinEdge()
{
    reset();
}

template<unsigned N>
void SimplexDCMinEdge<N>::reset()
{
    for (auto& simplex : simplices) {
        std::fill(simplex.intersections.begin(),
            simplex.intersections.end(),
            nullptr);
        simplex.vert.array() = 0.0;
        simplex.index = 0;
    }
    refcount = 0;
}

template<unsigned N>
void SimplexDCMinEdge<N>::releaseTo(Pool& object_pool)
{
    for (auto& simplex : simplices) {
        for (auto& i : simplex.intersections) {
            auto intersection = i.load();
            if (intersection == nullptr) {
                continue;
            }
            if (--intersection->refcount == 0) {
                object_pool.next().put(intersection);
            }
            i = nullptr;
        }
    }
    object_pool.put(this);
}

template<unsigned N>
inline void SimplexDCMinEdge<N>::EdgeVec::push_back(SimplexDCMinEdge* ptr)
{
    std::unique_lock lock(mut);
    vec.push_back(ptr);
}

template<unsigned N>
inline void SimplexDCMinEdge<N>::EdgeVec::sort()
{
    std::unique_lock lock(mut);
    if (isSorted) {
        return;
    }
    std::sort(vec.begin(), vec.end(), 
        [](auto a, auto b) {return a->lowPt < b->lowPt; });
    isSorted = true;
}

////////////////////////////////////////////////////////////////////////////////

template<unsigned N>
SimplexDCLeaf<N>::SimplexDCLeaf()
{
    reset();
}

template<unsigned N>
void SimplexDCLeaf<N>::reset()
{
    level = 0;
    tape.reset();
    std::fill(sub.begin(), sub.end(), nullptr);
    for (auto& edge : edges) {
        edge.emplace<SimplexDCMinEdge<N>::EdgeVec>();
    }
}

template<unsigned N>
void SimplexDCLeaf<N>::releaseTo(Pool& object_pool)
{
    for (auto& s : sub) {
        auto subspace = s.load();
        if (--subspace->refcount == 0) {
            object_pool.next().put(subspace);
        }
        s = nullptr;
    }
    for (auto& edge : edges) {
        if (std::holds_alternative<SimplexDCMinEdge<N>*>(edge)) {
            auto& minEdge = std::get<SimplexDCMinEdge<N>*>(edge);
            if (--minEdge->refcount == 0) {
                minEdge->releaseTo(object_pool.next().next());
            }
        }
        else {
            // Clear any vector or stack to free memory, by turning it
            // into the pointer version.
            edge.emplace<SimplexDCMinEdge<N>*>(nullptr);
        }
    }
    object_pool.put(this);
}

}   // namespace libfive
