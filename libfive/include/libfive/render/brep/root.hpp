/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include "libfive/render/brep/object_pool.hpp"
#include "libfive/render/brep/progress.hpp"

#define LIBFIVE_DEBUGDATA

#ifdef LIBFIVE_DEBUGDATA
#include <set>
#endif

namespace Kernel {

template <typename T>
class Root
{
public:
    Root() : ptr(nullptr) {}
    Root(T* ptr) : ptr(ptr) {}
    Root(Root&& other) { *this = std::move(other); }

    Root& operator=(Root&& other) {
        ptr = other.ptr;
        other.ptr = nullptr;
        object_pool = std::move(other.object_pool);
        tree_count = other.tree_count;
        return *this;
    }

    ~Root() { reset(); }

    void reset(unsigned workers=8,
               ProgressCallback progress_callback=EMPTY_PROGRESS_CALLBACK)
    {
        ptr = nullptr;

        std::atomic_bool done(false);
        std::atomic_bool cancel(false);
        auto progress_watcher = ProgressWatcher::build(
                object_pool.total_size(), 2.0f,
                progress_callback, done, cancel);

        object_pool.reset(workers, progress_watcher);
        done.store(true);
        delete progress_watcher;
    }

    const T* operator->() const { return ptr; }
    const T* get() const { return ptr; }

    void claim(typename T::Pool& pool) {
        tree_count += pool.size();
        object_pool.claim(pool);
    }

    int64_t size() const { return tree_count; }

    // Because the algorithm assigns trees to their parents only once
    // they are done, it can be difficult to determine what is
    // going wrong when they fail to ever finish.  This provides an
    // easier way to access them via a debugger, at a small performance cost.
#ifdef LIBFIVE_DEBUGDATA
    std::set<T*> debug_access_undone;
    std::set<T*> debug_access_done;

#endif
protected:
    T* ptr;
    typename T::Pool object_pool;

    // Used for progress tracking.  We use a signed value here because,
    // as we claim Pools of XTrees, it's possible for the intermediate
    // result to go negative (if one pool has claimed many trees from
    // another Pool, so it owns more trees than it has allocated)..
    int64_t tree_count=0;
};

}   // namespace Kernel
