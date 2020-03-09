/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "../simplex/simplex_tree.inl"
#include "simplex_dc_tree.inl"

namespace libfive {
template class XTree<3, SimplexTree<3, SimplexDCLeaf<3>>, SimplexDCLeaf<3>>;
}   // namespace libfive
