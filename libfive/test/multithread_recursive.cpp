/*
libfive: a CAD kernel for modeling with implicit functions

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <array>
#include <cassert>
#include <memory>
#include <type_traits>
#include <atomic>
#include <vector>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/task_group.h>

#include "catch.hpp"

#include "libfive/render/brep/multithread_recursive.hpp"

using namespace libfive;

namespace {

// Complete binary tree in heap numbering: node n has children 2n and 2n+1;
// nodes with n >= kLeaves are leaves.  Node 1 is the root.
constexpr int kLeaves = 8;
constexpr int kNodes = 2 * kLeaves;

// Divisible: a lightweight handle with getChildren().
struct Handle {
  int id;
  std::vector<Handle> getChildren() const {
    if (id >= kLeaves) return {};
    return {Handle{2 * id}, Handle{2 * id + 1}};
  }
};

struct Counters {
  std::array<std::atomic<int>, kNodes> pre{};
  std::array<std::atomic<int>, kNodes> post{};
  // Sequence stamp of post calls, to check ordering
  std::atomic<int> clock{0};
  std::array<std::atomic<int>, kNodes> postTime{};
  std::array<std::atomic<int>, kNodes> preTime{};
};

bool isAncestor(int a, int n) {  // a is an ancestor of (or equal to) n
  for (; n >= 1; n /= 2) if (n == a) return true;
  return false;
}

}  // namespace

TEST_CASE("multithreadRecursive (void): pre once per node, post after children") {
  Counters c;
  tbb::enumerable_thread_specific<int> locals;
  multithreadRecursive<Handle, void, int>(
      Handle{1}, locals,
      [&](Handle h, bool& recurse, int&) {
        c.pre[h.id]++;
        c.preTime[h.id] = ++c.clock;
        recurse = true;
      },
      [&](Handle h, int&) {
        c.post[h.id]++;
        c.postTime[h.id] = ++c.clock;
      });

  for (int n = 1; n < kNodes; ++n) {
    CAPTURE(n);
    REQUIRE(c.pre[n] == 1);
    if (n < kLeaves) {
      REQUIRE(c.post[n] == 1);
      REQUIRE(c.postTime[n] > c.postTime[2 * n]);
      REQUIRE(c.postTime[n] > c.postTime[2 * n + 1]);
    }
  }
  // Leaves have recurse == true but no children; post still runs for them
  // as well, once, after pre.
  for (int n = kLeaves; n < kNodes; ++n) {
    REQUIRE(c.post[n] == 1);
    REQUIRE(c.postTime[n] > c.preTime[n]);
  }
}

TEST_CASE("multithreadRecursive (void): recurse=false makes a leaf without post") {
  Counters c;
  tbb::enumerable_thread_specific<int> locals;
  multithreadRecursive<Handle, void, int>(
      Handle{1}, locals,
      [&](Handle h, bool& recurse, int&) {
        c.pre[h.id]++;
        recurse = (h.id != 2);
      },
      [&](Handle h, int&) { c.post[h.id]++; });
  REQUIRE(c.pre[2] == 1);
  REQUIRE(c.post[2] == 0);
  REQUIRE(c.pre[4] == 0);
  REQUIRE(c.post[1] == 1);
}

TEST_CASE("multithreadRecursive (Out): parent, childNo and root result") {
  struct Out { int id; int parentId; unsigned childNo; int sum; };
  std::array<std::atomic<int>, kNodes> parentOf{};
  std::array<std::atomic<int>, kNodes> childNoOf{};
  std::array<std::atomic<int>, kNodes> posted{};
  std::atomic<int> rootParentNull{0};
  tbb::enumerable_thread_specific<int> locals;

  Out result = multithreadRecursive<Handle, Out, int>(
      Handle{1}, locals,
      [&](Handle h, Out* par, unsigned childNo, bool& recurse, int&) {
        if (par == nullptr) {
          rootParentNull++;
        } else {
          parentOf[h.id] = par->id;
          childNoOf[h.id] = childNo;
        }
        recurse = true;
        return Out{h.id, par ? par->id : -1, childNo, 0};
      },
      [&](Handle h, Out& out, int&) {
        posted[h.id]++;
        out.sum = 1;  // marks that post saw this node's own Out
      });

  REQUIRE(rootParentNull == 1);
  REQUIRE(result.id == 1);
  REQUIRE(result.parentId == -1);
  REQUIRE(result.sum == 1);
  for (int n = 2; n < kNodes; ++n) {
    CAPTURE(n);
    REQUIRE(parentOf[n] == n / 2);
    REQUIRE(childNoOf[n] == unsigned(n % 2));
  }
  for (int n = 1; n < kNodes; ++n) REQUIRE(posted[n] == 1);
}

TEST_CASE("multithreadRecursive: tbb cancellation skips post on ancestors") {
  // pre() at node 5 (a leaf's parent, path 1 -> 2 -> 5) cancels the outer
  // task_group.  Node 5's wait, node 2's wait and the root's wait all
  // complete after the cancel, so none of them may run post.
  constexpr int kCancelAt = 5;
  Counters c;
  tbb::enumerable_thread_specific<int> locals;
  tbb::task_group outer;
  // Not run_and_wait(): in TBB 2019U6 that calls the functor inline on the
  // calling thread, outside any group task, so nested contexts would bind to
  // the default context and outer.cancel() would not reach them.
  outer.run([&] {
    multithreadRecursive<Handle, void, int>(
        Handle{1}, locals,
        [&](Handle h, bool& recurse, int&) {
          c.pre[h.id]++;
          recurse = true;
          if (h.id == kCancelAt) outer.cancel();
        },
        [&](Handle h, int&) { c.post[h.id]++; });
  });
  outer.wait();

  REQUIRE(c.pre[1] == 1);
  REQUIRE(c.pre[kCancelAt] == 1);
  for (int n = 1; n < kNodes; ++n) {
    if (isAncestor(n, kCancelAt)) {
      CAPTURE(n);
      REQUIRE(c.post[n] == 0);
    }
  }
}

TEST_CASE("multithreadRecursive (Out): tbb cancellation skips post on ancestors") {
  constexpr int kCancelAt = 5;
  std::array<std::atomic<int>, kNodes> posted{};
  tbb::enumerable_thread_specific<int> locals;
  tbb::task_group outer;
  // Not run_and_wait(): in TBB 2019U6 that calls the functor inline on the
  // calling thread, outside any group task, so nested contexts would bind to
  // the default context and outer.cancel() would not reach them.
  outer.run([&] {
    multithreadRecursive<Handle, int, int>(
        Handle{1}, locals,
        [&](Handle h, int*, unsigned, bool& recurse, int&) {
          recurse = true;
          if (h.id == kCancelAt) outer.cancel();
          return h.id;
        },
        [&](Handle h, int&, int&) { posted[h.id]++; });
  });
  outer.wait();
  for (int n = 1; n < kNodes; ++n) {
    if (isAncestor(n, kCancelAt)) {
      CAPTURE(n);
      REQUIRE(posted[n] == 0);
    }
  }
}
