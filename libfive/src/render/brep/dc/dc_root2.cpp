/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "../root.inl"
#include "libfive/render/brep/dc/dc_neighbors.hpp"
#include "libfive/render/brep/dc/dc_tree.hpp"

namespace libfive {
template class Root<DCTree<2>>;
template class Root<PseudoDCTree2>;
}   // namespace libfive
