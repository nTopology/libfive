/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <memory>

#include "libfive/render/brep/simplex/simplex_neighbors.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "libfive/render/brep/worker_pool.hpp"

namespace libfive {

template <unsigned N>
using SimplexDCWorkerPool = WorkerPool<SimplexTree<N, SimplexDCLeaf<N>>, 
                                       SimplexNeighbors<N, SimplexDCLeaf<N>>, 
                                       N,
                                       false>;

}   // namespace libfive

