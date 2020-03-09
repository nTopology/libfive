/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "../simplex/simplex_tree.inl"
#include "simplex_dc_tree.inl"
#include "../object_pool.inl"

namespace libfive {

int OrientationChecker<2>::check(
    const Eigen::Vector2d& a, const Eigen::Vector2d& b, 
    const Eigen::Vector2d& endpt1) const
{
    assert(endpt1 == endpoints.col(0) || endpt1 == endpoints.col(1));
    std::array<double, 2> xProds;
    for (auto i = 0; i < 2; ++i) {
        auto diff1 = a - endpoints.col(i);
        auto diff2 = b - endpoints.col(i);
        xProds[i] = diff1.x() * diff2.y() - diff2.x() * diff1.y();
        assert(!std::isnan(xProds[i]));
    }
    if (xProds[0] * xProds[1] < 0) {
        return 0;
    }
    // Each xProd has absolute value equal to twice the area of the triangle 
    // formed by a, b, and the corresponding endpoint, meaning that it is
    // the distance of that endpoint from the line between a and be multiplied
    // by half the length of ab, so a lower xProd value means the endpoint is
    // closer to that line, and thus the endpoint that we're past.
    assert(xProds[0] != xProds[1]);
    int out = (xProds[0] < xProds[1]) ? 1 : 2;
    if (endpt1 != endpoints.col(0)) {
        out = 3 - out;
    }
    return out;
}

template class SimplexTree<2, SimplexDCLeaf<2>>;
template struct SimplexDCLeaf<2>;
template struct SimplexDCIntersection<2>;
template struct SimplexDCMinEdge<2>;
template struct DCSimplex<2>;
template class SimplexDCMinEdge<2>::EdgeVec;
}   // namespace libfive
