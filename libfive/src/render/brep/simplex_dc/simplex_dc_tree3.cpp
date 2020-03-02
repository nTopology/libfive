/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "simplex_dc_tree.inl"
#include "../simplex/simplex_tree.inl"

namespace libfive {

bool OrientationChecker<3>::check(
    const Eigen::Vector3d& clockwisePt, 
    const Eigen::Vector3d& counterclockwisePt,
    const Eigen::Vector3d& center) const
{
    auto relativeClockwise = clockwisePt - center;
    auto relativeCC = counterclockwisePt - center;
    auto xProd = relativeClockwise.cross(relativeCC);
    // If center was 0, clockwisePt is X, counterclockwisePt is Y,
    // and the outward direction is Z, we want this to return true,
    // and in that case, xProd.dot(direction) is positive.
    return xProd.dot(direction) > 0.;
}

const SimplexDCIntersection<3> DCSimplex<3>::dupVertIntersection;

template class SimplexTree<3, SimplexDCLeaf<3>>;
template struct SimplexDCLeaf<3>;
template struct SimplexDCIntersection<3>;
template struct SimplexDCMinEdge<3>;
template struct DCSimplex<3>;
template class SimplexDCMinEdge<3>::EdgeVec;
}   // namespace libfive
