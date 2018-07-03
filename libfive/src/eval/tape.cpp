/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2017  Matt Keeter

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#include <unordered_map>

#include "libfive/eval/tape.hpp"
#include "libfive/eval/deck.hpp"

namespace Kernel {

Clause::Id Tape::rwalk(std::function<void(Opcode::Opcode, Clause::Id,
                                          Clause::Id, Clause::Id)> fn,
                       bool& abort)
{
    for (auto itr = t.rbegin(); itr != t.rend() && !abort; ++itr)
    {
        fn(itr->op, itr->id, itr->a, itr->b);
    }
    return i;
}

void Tape::walk(std::function<void(Opcode::Opcode, Clause::Id,
                                   Clause::Id, Clause::Id)> fn, bool& abort)
{
    for (auto itr = t.begin(); itr != t.end() && !abort; ++itr)
    {
        fn(itr->op, itr->id, itr->a, itr->b);
    }
}

std::shared_ptr<Tape>
Tape::push(const std::shared_ptr<Tape>& tape, Deck& deck,
           std::function<Keep(Opcode::Opcode, Clause::Id,
                              Clause::Id, Clause::Id)> fn,
           Type t, Region<3> r)
{
    // If this tape has no min/max clauses, then return it right away
    if (tape->terminal)
    {
        return tape;
    }

    // Since we'll be figuring out which clauses are disabled and
    // which should be remapped, we reset those arrays here
    std::fill(deck.disabled.begin(), deck.disabled.end(), true);
    std::fill(deck.remap.begin(), deck.remap.end(), 0);

    // Mark the root node as active
    deck.disabled[tape->i] = false;

    bool terminal = true;
    bool changed = false;
    for (const auto& c : tape->t)
    {
        if (!deck.disabled[c.id])
        {
            switch (fn(c.op, c.id, c.a, c.b))
            {
                case KEEP_A:        deck.disabled[c.a] = false;
                                    deck.remap[c.id] = c.a;
                                    changed = true;
                                    break;
                case KEEP_B:        deck.disabled[c.b] = false;
                                    deck.remap[c.id] = c.b;
                                    changed = true;
                                    break;
                case KEEP_BOTH:     terminal = false; // fallthrough
                case KEEP_ALWAYS:   break;
            }

            if (deck.remap[c.id])
            {
                deck.disabled[c.id] = true;
            }
            // Oracle nodes are special-cased here.  They should always
            // return either KEEP_BOTH or KEEP_ALWAYS, but have no children
            // to disable (and c.a is a dummy index into the oracles[]
            // array, so we shouldn't mis-interpret it as a clause index).
            else if (c.op != Opcode::ORACLE)
            {
                deck.disabled[c.a] = false;
                deck.disabled[c.b] = false;
            }
        }
    }

    if (!changed)
    {
        return tape;
    }

    std::shared_ptr<Tape> out;
    if (deck.spares.size())
    {
        out = deck.spares.back();
        deck.spares.pop_back();
    }
    else
    {
        out.reset(new Tape);
    }
    out->t.reserve(tape->t.size());

    out->type = t;
    out->parent = tape;
    out->terminal = terminal;
    out->t.clear(); // preserves capacity

    // Now, use the data in disabled and remap to make the new tape
    for (const auto& c : tape->t)
    {
        if (!deck.disabled[c.id])
        {
            // Oracle nodes use c.a as an index into tape->oracles,
            // rather than the address of an lhs / rhs expression,
            // so we special-case them here to avoid bad remapping.
            if (c.op == Opcode::ORACLE)
            {
                out->t.push_back({c.op, c.id, c.a, c.b});
            }
            else
            {
                Clause::Id ra, rb;
                for (ra = c.a; deck.remap[ra]; ra = deck.remap[ra]);
                for (rb = c.b; deck.remap[rb]; rb = deck.remap[rb]);
                out->t.push_back({c.op, c.id, ra, rb});
            }
        }
    }

    // Remap the tape root index
    for (out->i = tape->i; deck.remap[out->i]; out->i = deck.remap[out->i]);

    // Make sure that the tape got shorter
    assert(out->t.size() <= tape->t.size());

    // Store X / Y / Z bounds (may be irrelevant)
    out->X = {r.lower.x(), r.upper.x()};
    out->Y = {r.lower.y(), r.upper.y()};
    out->Z = {r.lower.z(), r.upper.z()};

    return out;
}

std::shared_ptr<Tape> Tape::getBase(
        std::shared_ptr<Tape> tape,
        const Eigen::Vector3f& p)
{
    // Walk up the tape stack until we find an interval-type tape
    // that contains the given point, or we hit the start of the stack
    while (tape->parent.get())
    {
        if (tape->type == Tape::INTERVAL &&
            p.x() >= tape->X.lower() && p.x() <= tape->X.upper() &&
            p.y() >= tape->Y.lower() && p.y() <= tape->Y.upper() &&
            p.z() >= tape->Z.lower() && p.z() <= tape->Z.upper())
        {
            break;
        }
        else
        {
            tape = tape->parent;
        }
    }

    return tape;
}

}   // namespace Kernel
