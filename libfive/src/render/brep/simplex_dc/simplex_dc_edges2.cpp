/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "simplex_dc_edges.inl"

namespace libfive {
template class SimplexDCEdges<2>;

template void SimplexDCEdges<2>::load<Axis::X>(
    const std::array<Input*, 2>&);
template void SimplexDCEdges<2>::load<Axis::Y>(
    const std::array<Input*, 2>&);
template void SimplexDCEdges<2>::load<Axis::Z>(
    const std::array<Input*, 2>&);
}   // namespace libfive
