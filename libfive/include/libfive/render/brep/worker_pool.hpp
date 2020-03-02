/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <atomic>
#include <boost/lockfree/stack.hpp>

#include "libfive/render/brep/root.hpp"
#include "libfive/tree/tree.hpp"

namespace libfive {

// Forward declarations
class Evaluator;
template <unsigned N> class Region;
class Tape;
struct BRepSettings;
class VolTree;

/*
 *  A WorkerPool is used to construct a recursive tree (quadtree / octree)
 *  by sharing the work among a pool of threads.
 */
template <typename T, typename Neighbors, unsigned N, bool constOut = true>
class WorkerPool
{
public:
    using OutType = Root<std::conditional_t<constOut, const T, T>>;
    /*
     *  Evaluation function that builds a local evaluator array
     *  (based on settings.workers)
     */
    static OutType build(Tree t, const Region<N>& region,
                         const BRepSettings& settings);
    /*
     *  General-purpose evaluation function
     *
     *  eval must be an array of at least [settings.workers] evaluators
     */
    static OutType build(Evaluator* eval, const Region<N>& region,
                         const BRepSettings& settings);

protected:
    struct Task {
        T* target;
        std::shared_ptr<Tape> tape;
        Neighbors parent_neighbors;
        const VolTree* vol;
    };

    using LockFreeStack =
        boost::lockfree::stack<Task, boost::lockfree::fixed_sized<true>>;

    static void run(Evaluator* eval, LockFreeStack& tasks,
                    OutType& root, std::mutex& root_lock,
                    const BRepSettings& settings,
                    std::atomic_bool& done);
};

}   // namespace libfive
