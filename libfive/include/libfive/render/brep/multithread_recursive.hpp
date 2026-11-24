/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <cassert>
#include <type_traits>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/task_group.h>

namespace libfive {

namespace detail {

  // CallerCtx: cancellation shim for TBB 2019 (TBB_INTERFACE_VERSION < 12000).
  // DELETE THIS WHOLE STRUCT once turbo is on oneTBB (nTopology/turbo#19518,
  // the oneTBB switch); the recursion then relies on tg.wait() alone.
  //
  // Why it exists: the old tbb::task implementation shared one
  // task_group_context between every node of the recursion, so cancelling
  // the caller's group/context stopped the whole tree.  After the port to
  // tbb::task_group each branch node owns a task_group.  On TBB 2019U6, a
  // task_group constructed inside a task that the master thread runs
  // directly from its outermost wait gets an isolated context, so cancelling
  // an enclosing group/context never reaches these nested groups and
  // tg.wait() would not report it.
  //
  // What it does: multithreadRecursive() (the public entry point) constructs
  // one CallerCtx, which captures the task_group_context of the task that is
  // executing at that moment (tbb::task::self().group()).  The same object
  // is copied down to every node, and each node checks cancelled() in
  // addition to `tg.wait() == tbb::canceled` before calling post().  It is
  // captured once at entry, so it only reflects the context seen there; it
  // does not track contexts of groups created further down.
  //
  // tbb::task reaches this header via <tbb/task_group.h> on 2019U6.
#if TBB_INTERFACE_VERSION < 12000
  struct CallerCtx {
    tbb::task_group_context* ctx = tbb::task::self().group();
    bool cancelled() const { return ctx->is_group_execution_cancelled(); }
  };
#else
  // oneTBB: a task_group_context binds to its parent context by default, so
  // cancelling an enclosing group propagates to the nested groups and
  // tg.wait() reports it.  Nothing extra to check, hence always false.
  struct CallerCtx {
    bool cancelled() const { return false; }
  };
#endif

  template <class Divisible, class Out, class Local,
    class PreCallable, class PostCallable>
    Out multithreadRecursiveImpl(
      Divisible root, tbb::enumerable_thread_specific<Local>& locals,
      const PreCallable& pre, const PostCallable& post,
      Out* parent, unsigned childNo, const CallerCtx& callerCtx);

}   // namespace detail


  // Calls Callable on Divisible, recursing over its children.  

  // Divisible should be either a (possibly const) class or a pointer to one; in
  // either case, the class should have a method getChildren(), returning an 
  // iterable container of Divisibles with a size() method.

  // If Out is not void, PreCallable should have the signature 
  // Out(Divisible, Out* par, unsigned childNo, bool& recurse, Local&); 
  // par will be the parent of the current node, or nullptr if the current node 
  // is the root, and childNo will be which number child of that parent is
  // currently being handled (or undefined if the current node is the root).  
  // "recurse" must be set by PreCallable to instruct this algorithm whether to
  // continue recursing or to make the current node into a leaf.  While 
  // PreCallable may modify the Out that is passed in, it must do so atomically
  // if at all.
  // If Out is void, the signature should instead be
  // void(Divisible, bool& recurse, Local&).  In either case, the call must be
  // thread-safe with respect to any captured variables (i.e. only writing to 
  // atomics, if at all), and its return value is also the return value of this
  // function.

  // PostCallable should have the signature void(Divisible, Out& child, Local&)
  // if Out is not void, or void(Divisible, Local&) if it is.  PreCallable will be
  // called before recursing to children, while PostCallable will be called
  // afterward; there is no guarantee that the two will be called on the same 
  // thread (or with the same Local&).  If "recurse" is set to false by 
  // PreCallable, or if the algorithm was cancelled via tbb, PostCallable will 
  // not be called.

  // The final two function arguments (parent and childNo) should not be passed in
  // except via recursive calling.

  template <class Divisible, class Out, class Local,
    class PreCallable, class PostCallable>
    Out multithreadRecursive(
      Divisible root, tbb::enumerable_thread_specific<Local>& locals,
      const PreCallable& pre, const PostCallable& post,
      Out* parent = nullptr, unsigned childNo = -1) {
    return detail::multithreadRecursiveImpl(
      root, locals, pre, post, parent, childNo, detail::CallerCtx());
  }

namespace detail {

  template <class Divisible, class Out, class Local,
    class PreCallable, class PostCallable>
    Out multithreadRecursiveImpl(
      Divisible root, tbb::enumerable_thread_specific<Local>& locals,
      const PreCallable& pre, const PostCallable& post,
      Out* parent, unsigned childNo, const CallerCtx& callerCtx) {

    struct Task {
      Divisible root;
      tbb::enumerable_thread_specific<Local>& locals;
      const PreCallable& pre;
      const PostCallable& post;
      Out* parent;
      unsigned childNo;
      CallerCtx callerCtx;

      Task(Divisible root, tbb::enumerable_thread_specific<Local>& locals,
        const PreCallable& pre, const PostCallable& post,
        Out* parent, unsigned childNo, const CallerCtx& callerCtx)
        : root(root), locals(locals), pre(pre), post(post),
        parent(parent), childNo(childNo), callerCtx(callerCtx) {}

      void operator()() const {
        multithreadRecursiveImpl(
          root, locals, pre, post, parent, childNo, callerCtx);
      }
    };

    bool toRecurse;
    auto getChildren = [&root]() {
      if constexpr (std::is_pointer_v<Divisible>) {
        return root->getChildren();
      }
      else {
        return root.getChildren();
      }
    };

    if constexpr (std::is_same_v<Out, void>) {
      pre(root, toRecurse, locals.local());
      if (toRecurse) {
        auto children = getChildren();
        tbb::task_group tg;
        for (auto& child : children) {
          // We don't need to set parent and childNo, since when Out is
          // void they are not used.
          tg.run(Task(child, locals, pre, post, nullptr, -1, callerCtx));
        }
        // Skip post() if this group was cancelled, or (TBB 2019 only, see
        // CallerCtx) if the caller's context was.
        if (tg.wait() != tbb::canceled && !callerCtx.cancelled()) {
          post(std::move(root), locals.local());
        }
      }
    }
    else {
      auto out = pre(root, parent, childNo, toRecurse, locals.local());
      if (toRecurse) {
        auto children = getChildren();
        tbb::task_group tg;
        auto i = 0;
        for (auto& child : children) {
          tg.run(Task(child, locals, pre, post, &out, i++, callerCtx));
        }
        assert(i == children.size());
        if (tg.wait() != tbb::canceled && !callerCtx.cancelled()) {
          post(std::move(root), out, locals.local());
        }
      }
      return out;
    }
  }

}   // namespace detail

}   // namespace libfive
