/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <stack>
#include <boost/lockfree/stack.hpp>

#include "libfive/render/brep/free_thread_handler.hpp"
#include "libfive/render/brep/settings.hpp"
#include "libfive/render/brep/mesh.hpp"
#include "libfive/render/brep/root.hpp"

#include "libfive/render/axes.hpp"
#include "libfive/eval/interval.hpp"

namespace libfive {

/*
 *  Class to walk a dual grid for a quad or octree
 */
template <unsigned N>
class Dual
{
public:
     /*
      *  Basic dual-walking function
      *
      *  The walker type W (often, but not necessarily, a mesher) needs
      *     load<Axis::Axis A>(const std::array<T*, 1 << K>& trees) for 
      *         0 < K < N, and similar untemplated methods for K == 0 and 
      *         k == N; the one for K == 0 should take a T* rather than a const
      *         array of size 1.  A "true" dual-walker will have these trivial 
      *         for K != N - 1, but this class allows for other values of K in 
      *         order to walk faces (in 3d) and corners as well as edges.  T 
      *         may be constant, and should be if the load method will not 
      *         output directly to the input object.  The indices in the arrays 
      *         will correspond to dimensions in the positions of the cells, 
      *         following the order (from lowest to highest bit) Q(A), R(A) for
      *         the templated versions and x, y, z (or just x, y in 2d) for
      *         K == N.
      *     Input (typename, the tree to be walked)
      *     PerThreadOutput (typename)
      *     Output (typename)
      *  and must have a constructor of the form
      *     W(PerThreadOutput&, A...).
      *  PerThreadOutput must be constructible from std::atomic<uint32_t>&, and
      *  Output must have a method collect(std::vector<PerThreadOutput>).
      */
    template<typename W, typename ... A>
    static std::unique_ptr<typename W::Output> walk(
            const Root<typename W::Input>& t,
            const BRepSettings& settings,
            A... args);

     /*
      *  Flexible dual-walking function
      *
      *  The walker type W needs
      *     load(const std::array<T*, N>& trees)
      *     Input (typename)
      *     PerThreadOutput (typename)
      *     Output (typename)
      *
      *  The factory can be anything that spits out valid W objects,
      *  given a PerThreadOutput and worker index.
      */
    template<typename W>
    static std::unique_ptr<typename W::Output> walk_(
            const Root<typename W::Input>& t,
            const BRepSettings& settings,
            std::function<W(typename W::PerThreadOutput&, int)> WalkerFactory);

protected:
    template<typename T, typename Walker>
    static void run(Walker& m,
                    boost::lockfree::stack<T*,
                                           boost::lockfree::fixed_sized<true>>& tasks,
                    const BRepSettings& settings,
                    std::atomic_bool& done);

    template <typename T, typename Walker>
    static void work(T* t, Walker& m);

    template <typename T, typename Walker>
    static void handleTopEdges(T* t, Walker& m);
};

////////////////////////////////////////////////////////////////////////////////
// 2D Implementation

template <typename T, typename V>
void corner2(const std::array<T*, 4> & ts, V& v)
{

    if (std::any_of(ts.begin(), ts.end(),
        [](T* t) { return t->isBranch(); }))
    {
        corner2<T, V>({{ts[0]->child(3), ts[1]->child(2), 
                        ts[2]->child(1), ts[3]->child(0)}}, v);
    }
    else
    {
        v.load(ts);
    }
}

template <typename T, typename V, Axis::Axis A>
void edge2(const std::array<T*, 2>& ts, V& v)
{
    constexpr uint8_t perp = (Axis::X | Axis::Y) ^ A;

    if (std::any_of(ts.begin(), ts.end(),
        [](T* t){ return t->isBranch(); }))
    {
        edge2<T, V, A>({{ts[0]->child(perp), ts[1]->child(0)}}, v);
        edge2<T, V, A>({{ts[0]->child(A|perp), ts[1]->child(A)}}, v);
        std::array<T*, 4> cornerTs{{ts[0]->child(perp), ts[0]->child(A | perp),
                                    ts[1]->child(0), ts[1]->child(A)}};
        if constexpr (A == Axis::Y) {
            std::swap(cornerTs[1], cornerTs[2]);
        }
        corner2(cornerTs, v);
    }
    else
    {
        v.template load<A>(ts);
    }
}

template <>
template <typename T, typename V>
void Dual<2>::work(T* t, V& v)
{
    if (!t->isBranch()) {
        v.load(t);
    }
    edge2<T, V, Axis::Y>({{t->child(0), t->child(Axis::X)}}, v);
    edge2<T, V, Axis::Y>({{t->child(Axis::Y), t->child(Axis::Y | Axis::X)}}, v);
    edge2<T, V, Axis::X>({{t->child(0), t->child(Axis::Y)}}, v);
    edge2<T, V, Axis::X>({{t->child(Axis::X), t->child(Axis::X | Axis::Y)}}, v);
    corner2<T, V>({{t->child(0), t->child(Axis::X), 
                    t->child(Axis::X), t->child(Axis::X | Axis::Y)}}, v);
}

template <>
template <typename T, typename Walker>
void Dual<2>::handleTopEdges(T* t, Walker& m)
{
    (void)t;
    (void)m;

    // TODO
    // No one should be calling this yet, because simplex meshing
    // isn't implemented in 2D.
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////
// 3D Implementation

template <typename T, typename V>
void corner3(const std::array<T*, 8>& ts, V& v)
{
    if (std::any_of(ts.begin(), ts.end(),
        [](T* t) { return t->isBranch(); }))
    {
        corner3<T, V>({{ts[0]->child(7), ts[1]->child(6),
                        ts[2]->child(5), ts[3]->child(4),
                        ts[4]->child(3), ts[5]->child(2),
                        ts[6]->child(1), ts[7]->child(0)}}, v);
    }
    else
    {
        v.load(ts);
    }
}

template <typename T, typename V, Axis::Axis A>
void edge3(const std::array<T*, 4> ts, V& v)
{
    constexpr auto Q = Axis::Q(A);
    constexpr auto R = Axis::R(A);

    if (std::any_of(ts.begin(), ts.end(),
        [](T* t){ return t->isBranch(); }))
    {
        edge3<T, V, A>({{ts[0]->child(Q|R), ts[1]->child(R), ts[2]->child(Q), ts[3]->child(0)}}, v);
        edge3<T, V, A>({{ts[0]->child(Q|R|A), ts[1]->child(R|A), ts[2]->child(Q|A), ts[3]->child(A)}}, v);
        std::array<T*, 8> cornerTs{{ts[0]->child(Q|R), ts[0]->child(Q|R|A),
                                    ts[1]->child(R), ts[1]->child(R|A),
                                    ts[2]->child(Q), ts[2]->child(Q|A),
                                    ts[3]->child(0), ts[3]->child(A)}};
        // The above order works when A is X, but for Y or Z we need to rotate to make it be in the
        // correct order.
        auto rotate = [&cornerTs]() {
            std::swap(cornerTs[1], cornerTs[2]);
            std::swap(cornerTs[1], cornerTs[4]);
            std::swap(cornerTs[6], cornerTs[5]);
            std::swap(cornerTs[6], cornerTs[3]);
        };
        if constexpr (A == Axis::Y) {
            rotate();
        }
        else if constexpr (A == Axis::Z) {
            // Perhaps not that efficient, but it should be optimized away.
            rotate();
            rotate();
        }
        corner3(cornerTs, v);
    }
    else
    {
        v.template load<A>(ts);
    }
}

template <typename T, typename V, Axis::Axis A>
void face3(const std::array<T*, 2> ts, V& v)
{
    if (std::any_of(ts.begin(), ts.end(),
        [](T* t){ return t->isBranch(); }))
    {
        constexpr auto Q = Axis::Q(A);
        constexpr auto R = Axis::R(A);

        for (unsigned k : {0, (int)Q, (int)R, Q|R})
        {
            face3<T, V, A>({{ts[0]->child(k|A), ts[1]->child(k)}}, v);
        }

        edge3<T, V, Q>({{ts[0]->child(A), ts[0]->child(R|A), ts[1]->child(0), ts[1]->child(R)}}, v);
        edge3<T, V, Q>({{ts[0]->child(Q|A), ts[0]->child(Q|R|A), ts[1]->child(Q), ts[1]->child(Q|R)}}, v);

        edge3<T, V, R>({{ts[0]->child(A), ts[1]->child(0), ts[0]->child(A|Q), ts[1]->child(Q)}}, v);
        edge3<T, V, R>({{ts[0]->child(R|A), ts[1]->child(R), ts[0]->child(R|A|Q), ts[1]->child(R|Q)}}, v);

        std::array<T*, 8> cornerTs{{ts[0]->child(A), ts[1]->child(0),
                                    ts[0]->child(Q|A), ts[1]->child(Q),
                                    ts[0]->child(R|A), ts[1]->child(R),
                                    ts[0]->child(Q|R|A), ts[1]->child(Q|R)}};
        // The above order works when A is X, but for Y or Z we need to rotate to make it be in the
        // correct order.
        auto rotate = [&cornerTs]() {
            std::swap(cornerTs[1], cornerTs[2]);
            std::swap(cornerTs[1], cornerTs[4]);
            std::swap(cornerTs[6], cornerTs[5]);
            std::swap(cornerTs[6], cornerTs[3]);
        };
        if constexpr (A == Axis::Y) {
            rotate();
        }
        else if constexpr (A == Axis::Z) {
            // Perhaps not that efficient, but it should be optimized away.
            rotate();
            rotate();
        }
        corner3(cornerTs, v);
    }
    else
    {
        v.template load<A>(ts);
    }
}

template <typename T, typename V, Axis::Axis A>
void call_edge3(T* t, V& v)
{
    for (auto a : {Axis::Axis(0), A})
    {
        edge3<T, V, A>({{t->child(a),
             t->child(Axis::Q(A) | a),
             t->child(Axis::R(A) | a),
             t->child(Axis::Q(A) | Axis::R(A) | a)}}, v);
    }
}

template <typename T, typename V, Axis::Axis A>
void call_face3(T* t, V& v)
{
    constexpr auto q = Axis::Q(A);
    constexpr auto r = Axis::R(A);

    face3<T, V, A>({{t->child(0), t->child(A)}}, v);
    face3<T, V, A>({{t->child(q), t->child(q|A)}}, v);
    face3<T, V, A>({{t->child(r), t->child(r|A)}}, v);
    face3<T, V, A>({{t->child(q|r), t->child(q|r|A)}}, v);
}

template <>
template <typename T, typename V>
void Dual<3>::work(T* t, V& v)
{
    // Handle the cell itself, if it is a leaf.
    if (!t->isBranch()) {
        v.load(t);
    }

    // Call the face procedure on every pair of cells (4x per axis)
    call_face3<T, V, Axis::X>(t, v);
    call_face3<T, V, Axis::Y>(t, v);
    call_face3<T, V, Axis::Z>(t, v);

    // Call the edge function 6 times (2x per axis)
    call_edge3<T, V, Axis::X>(t, v);
    call_edge3<T, V, Axis::Y>(t, v);
    call_edge3<T, V, Axis::Z>(t, v);

    // Handle the center corner.
    corner3<T, V>({{t->child(0), t->child(1), t->child(2), t->child(3),
                    t->child(4), t->child(5), t->child(6), t->child(7)}}, v);
}

template <>
template <typename T, typename V>
void Dual<3>::handleTopEdges(T* t, V& v)
{
    auto e = T::empty();

    for (unsigned i=0; i < 8; ++i)
    {
        std::array<T*, 8> ts = {{e.get(), e.get(), e.get(), e.get(),
                                 e.get(), e.get(), e.get(), e.get()}};
        ts[i] = t;
        corner3<T, V>(ts, v);
    }

    for (unsigned i=0; i < 4; ++i)
    {
        std::array<T*, 4> ts = {{e.get(), e.get(), e.get(), e.get()}};
        ts[i] = t;
        edge3<T, V, Axis::X>(ts, v);
        edge3<T, V, Axis::Y>(ts, v);
        edge3<T, V, Axis::Z>(ts, v);
    }

    for (unsigned i=0; i < 2; ++i)
    {
        std::array<T*, 2> ts = {{e.get(), e.get()}};
        ts[i] = t;
        face3<T, V, Axis::X>(ts, v);
        face3<T, V, Axis::Y>(ts, v);
        face3<T, V, Axis::Z>(ts, v);
    }
}

////////////////////////////////////////////////////////////////////////////////

template <unsigned N>
template<typename W, typename ... A>
std::unique_ptr<typename W::Output> Dual<N>::walk(
            const Root<typename W::Input>& t,
            const BRepSettings& settings,
            A... args)
{
    return walk_<W>(
            t, settings,
            [&args...](typename W::PerThreadOutput& brep, int i) {
                (void)i;
                return W(brep, args...);
                });

}

template <unsigned N>
template<typename W>
std::unique_ptr<typename W::Output> Dual<N>::walk_(
            const Root<typename W::Input>& t,
            const BRepSettings& settings,
            std::function<W(typename W::PerThreadOutput&, int)> WalkerFactory)
{
    boost::lockfree::stack<
        typename W::Input*,
        boost::lockfree::fixed_sized<true>> tasks(settings.workers);
    tasks.push(t.get());
    t->resetPending();

    std::atomic<uint32_t> global_index(1);
    std::vector<typename W::PerThreadOutput> perThread;
    for (unsigned i=0; i < settings.workers; ++i) {
        perThread.emplace_back(typename W::PerThreadOutput(global_index));
    }

    if (settings.progress_handler) {
        settings.progress_handler->nextPhase(t.size() + 1);
    }

    std::vector<std::future<void>> futures;
    futures.resize(settings.workers);
    std::atomic_bool done(false);
    for (unsigned i=0; i < settings.workers; ++i) {
        futures[i] = std::async(std::launch::async,
            [&perThread, &tasks, &WalkerFactory, &settings, &done, i]()
            {
                auto m = WalkerFactory(perThread[i], i);
                Dual<N>::run(m, tasks, settings, done);
            });
    }

    // Wait on all of the futures
    for (auto& f : futures) {
        f.get();
    }

    assert(done.load() || settings.cancel.load());

    // Handle the top tree edges (only used for simplex meshing)
    if (W::needsTopEdges()) {
        auto m = WalkerFactory(perThread[0], 0);
        Dual<N>::handleTopEdges(t.get(), m);
    }

    auto out = std::unique_ptr<typename W::Output>(new typename W::Output);
    out->collect(perThread);
    return out;
}


template <unsigned N>
template <typename T, typename V>
void Dual<N>::run(V& v,
                  boost::lockfree::stack<T*,
                                         boost::lockfree::fixed_sized<true>>& tasks,
                  const BRepSettings& settings,
                  std::atomic_bool& done)

{
    // Tasks to be evaluated by this thread (populated when the
    // MPMC stack is completely full).
    std::stack<T*, std::vector<T*>> local;

    while (!done.load() && !settings.cancel.load())
    {
        // Prioritize picking up a local task before going to
        // the MPMC queue, to keep things in this thread for
        // as long as possible.
        T* t;
        if (local.size())
        {
            t = local.top();
            local.pop();
        }
        else if (!tasks.pop(t))
        {
            t = nullptr;
        }

        // If we failed to get a task, keep looping
        // (so that we terminate when either of the flags are set).
        if (t == nullptr)
        {
            if (settings.free_thread_handler != nullptr) {
                settings.free_thread_handler->offerWait();
            }
            continue;
        }

        if (t->isBranch())
        {
            // Recurse, calling the cell procedure for every child
            for (const auto& c_ : t->children)
            {
                const auto c = c_.load();
                if (!tasks.bounded_push(c)) {
                    local.push(c);
                }
            }
            continue;
        }

        // Special-case for singleton trees, which have null parents
        // (and have already been subtracted from pending)
        if (T::isSingleton(t)) {
            continue;
        }

        if (settings.progress_handler) {
            settings.progress_handler->tick();
        }

        for (t = t->parent; t && t->pending-- == 0; t = t->parent)
        {
            // Do the actual work (specialized for N = 2 or 3)
            Dual<N>::work(t, v);

            // Report trees as completed
            if (settings.progress_handler) {
                settings.progress_handler->tick();
            }
        }

        // Termination condition:  if we've ended up pointing at the parent
        // of the tree's root (which is nullptr), then we're done and break
        if (t == nullptr) {
            break;
        }
    }

    // If we've broken out of the loop, then we should set the done flag
    // so that other worker threads also terminate.
    done.store(true);
}

}   // namespace libfive
