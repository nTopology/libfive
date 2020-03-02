/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "libfive/render/brep/indexes.hpp"

namespace libfive {

/*  Free function for debugging; throws if there is a problem with the subs.
 *  Currently only tests that collapsing properly avoids collapsing "past"
 *  an intermediate subspace that does not collapse to the same thing. */
template <unsigned N>
void testSubs(const SimplexDCTree<N>& tree) {
    if (tree.isBranch()) {
        for (auto i = 0; i < (1 << N); ++i) {
            testSubs(*tree.child(i));
        }
    }
    else {
        const auto& subs = tree.leaf->sub;
        auto cellCollapse = (N == 3) 
            ? NeighborIndex(subs[26].load()->collapseRef) 
            : NeighborIndex(0);
        // Check that faces do not collapse past edges.
        for (auto faceAxisIdx = 0; faceAxisIdx < N; ++faceAxisIdx) {
            if (N == 2) {
                faceAxisIdx = 2;
            }
            auto faceAxis = Axis::toAxis(faceAxisIdx);
            auto floating = 7 ^ faceAxis;
            for (auto pos : { 0, static_cast<int>(faceAxis) }) {
                auto index = NeighborIndex::fromPosAndFloating(pos, floating).i;
                auto neighbor = tree.neighbor(index);
                if (neighbor && neighbor->isBranch()) {
                    // The face is not minimal for this cell; almost nothing 
                    // can be established here.
                    if (cellCollapse.i == index) {
                        throw &tree;
                    }
                    else {
                        continue;
                    }
                }
                const auto& sub = subs[index];
                NeighborIndex collapse(sub.load()->collapseRef);
                if (N == 3 && !(faceAxis & cellCollapse.floating()) &&
                    !(faceAxis & (pos ^ cellCollapse.pos()))) {
                    for (auto compAxisIdx = 0; compAxisIdx < 2; ++compAxisIdx) {
                        auto compAxis = (compAxisIdx == 0) ? Axis::Q(faceAxis)
                                                           : Axis::R(faceAxis);
                        if (collapse.i >= 27) {
                            throw &tree;
                        }
                        if (bool(compAxis & cellCollapse.pos()) !=
                            bool((1 << compAxisIdx) & collapse.pos())) {
                            throw &tree;
                        }
                        if (bool(compAxis & cellCollapse.floating()) !=
                            bool((1 << compAxisIdx) & collapse.floating())) {
                            throw &tree;
                        }
                    }
                }
                if (collapse.floating() != 0) {
                    assert(collapse.dimension() == 1 || 
                           collapse.dimension() == 2);
                    continue;
                }
                for (auto edgeAxisIdx = 0; edgeAxisIdx < 2; ++edgeAxisIdx) {
                    auto edgeAxis = edgeAxisIdx == 0 ? Axis::Q(faceAxis) 
                                                     : Axis::R(faceAxis);
                    auto thirdAxis = edgeAxisIdx == 0 ? Axis::R(faceAxis) 
                                                      : Axis::Q(faceAxis);
                    auto edgeRelativePos = edgeAxisIdx == 0 ? collapse.i / 3 
                                                            : collapse.i % 3;
                    auto posOnEdge = edgeAxisIdx == 0 ? collapse.i % 3
                                                      : collapse.i / 3;
                    assert(edgeRelativePos == 0 || edgeRelativePos == 1);
                    assert(posOnEdge == 0 || posOnEdge == 1);
                    auto edgePos = pos;
                    if (edgeRelativePos) {
                        edgePos |= thirdAxis;
                    }
                    else {
                        edgePos &= ~thirdAxis;
                    }
                    auto edgeIdx = NeighborIndex::fromPosAndFloating(
                        edgePos, edgeAxis).i;
                    auto edgeCollapse = subs[edgeIdx].load()->collapseRef;
                    if (edgeCollapse >= 2) {
                        assert(edgeCollapse == 2);
                        throw &tree;
                    }
                    if (edgeCollapse != posOnEdge) {
                        throw &tree;
                    }
                }
            }
        }
    }
}

}   // namespace libfive
