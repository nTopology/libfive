/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "simplex_dc_vertexer.inl"

namespace libfive {
template class SimplexDCVertexer<3, false>;
template class SimplexDCVertexer<3, true>;


template void SimplexDCVertexer<3, false>::load<Axis::X>(
    const std::array<Input*, 4>&);
template void SimplexDCVertexer<3, false>::load<Axis::Y>(
    const std::array<Input*, 4>&);
template void SimplexDCVertexer<3, false>::load<Axis::Z>(
    const std::array<Input*, 4>&);
template void SimplexDCVertexer<3, true>::load<Axis::X>(
    const std::array<Input*, 4>&);
template void SimplexDCVertexer<3, true>::load<Axis::Y>(
    const std::array<Input*, 4>&);
template void SimplexDCVertexer<3, true>::load<Axis::Z>(
    const std::array<Input*, 4>&);
}   // namespace libfive
