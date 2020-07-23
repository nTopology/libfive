/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "simplex_dc_collapser.inl"
#include "simplex_dc_tree2.cpp"

namespace libfive {
template class SimplexDCCollapser<2, 2>;

// Face loader in 2d should be used only for axis Z (as cell loader).
template <>
template <>
void SimplexDCCollapser<2, 2>::load<Axis::X>(
    const std::array<Input*, 1>&) { assert(false); }

template <>
template <>
void SimplexDCCollapser<2, 2>::load<Axis::Y>(
    const std::array<Input*, 1>&) { assert(false); }

template void SimplexDCCollapser<2, 2>::load<Axis::X>(
    const std::array<Input*, 2>&);
template void SimplexDCCollapser<2, 2>::load<Axis::Y>(
    const std::array<Input*, 2>&);
}   // namespace libfive
