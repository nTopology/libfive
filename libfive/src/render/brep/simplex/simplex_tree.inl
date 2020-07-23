/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <future>
#include <numeric>
#include <functional>
#include <limits>

#include <cmath>

#include <Eigen/StdVector>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/stack.hpp>

#include "libfive/render/brep/simplex/simplex_tree.hpp"
#include "libfive/render/brep/simplex/simplex_neighbors.hpp"

#include "libfive/render/brep/region.hpp"
#include "libfive/render/brep/settings.hpp"
#include "libfive/render/brep/neighbor_tables.hpp"

#include "libfive/render/axes.hpp"
#include "libfive/eval/tape.hpp"

#include "../xtree.inl"

namespace libfive {



}   // namespace libfive
