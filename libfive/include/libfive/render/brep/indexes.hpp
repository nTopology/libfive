/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <cstdint>
#include <assert.h>
#include "libfive/render/brep/util.hpp"

namespace libfive {

struct CornerIndex; // Forward declaration

/*
 *  A NeighborIndex in an N-dimension space is an N-digit ternary number.
 *  Each digit refers to a particular axis, and is interpreted as
 *      0:  On the lower extreme of that axis
 *      1:  On the upper extreme of that axis
 *      2:  Spanning that index
 *
 *  For example, the subcells in a 2D space are numbered as follows
 *  (left in decimal, right in ternary)
 *
 *    3---5---4      10---12---11
 *    |       |      |          |
 *    6   8   7      20   22   21
 *    |       |      |          |
 *    0---2---1      00---02---01
 */
struct NeighborIndex {
    constexpr NeighborIndex() : i(0) {}
    constexpr NeighborIndex(unsigned i) : i(i) {}

    /*  Calling dimension() on a constexpr NeighborIndex and using it as a 
     *  constexpr value seems to cause an internal compiler error in Visual 
     *  Studio, so we can use this helper static method for a workaround.
     */

    static constexpr unsigned dimension(unsigned i) {
        return i
            ? ((i % 3) == 2) + NeighborIndex(i / 3).dimension()
            : 0;
    }

    constexpr unsigned dimension() const {
        return dimension(i);
    }

    // Defined below, after struct CornerIndex is defined
    constexpr bool contains(const CornerIndex& c) const;

    constexpr bool contains(const NeighborIndex& n) const
    {
        return i
            ?   ((i % 3) == 2 ||
                 (i % 3) == (n.i % 3))
                && NeighborIndex(i / 3).contains(NeighborIndex(n.i / 3))
            : n.i == 0;
    }

    NeighborIndex operator|(const NeighborIndex& other) const {
        return NeighborIndex(
            (i || other.i)
            ? 3 * (NeighborIndex(i / 3) | NeighborIndex(other.i / 3)).i +
                ((i % 3) != (other.i % 3) ? 2 : (i % 3))
            : 0);
    }

    /*
     *  Returns a bitfield representing which axes are floating
     */
    constexpr uint8_t floating() const {
        return i
            ? (NeighborIndex(i / 3).floating() << 1) | ((i % 3) == 2)
            : 0;
    }

    /*
     *  Returns a bitfield representing which axes are fixed
     *
     *  Note that axes above the relevant dimension are all marked
     *  as fixed, because they're 0s in ternary.
     */
    constexpr uint8_t fixed() const { return ~floating(); }

    /*
     *  Checks to see whether a particular axis is fixed
     *  a is 1, 2, or 3 for X, Y, or Z axis respectively.
     */
    constexpr bool isAxisFixed(uint8_t a) const
    { return fixed() & (1 << a); }

    constexpr bool axisPosition(uint8_t a) const
    { return pos() & (1 << a); }

    /*
     *  Returns a bitfield representing which axes are fixed,
     *  with extra-dimension axes masked.
     */
    template <unsigned N>
    constexpr uint8_t fixed() const { return fixed() & ((1 << N) - 1); }

    /*
     *  Returns a bitfield representing low vs high on fixed axes
     *  For floating axes, pos & axis == 0, but that's meaningless.
     */
    constexpr uint8_t pos() const {
        return i
            ? (NeighborIndex(i / 3).pos() << 1) | ((i % 3) == 1)
            : 0;
    }

    /*
     *  Checks whether this is index refers to a corner
     */
    constexpr bool isCorner() const {
        return floating() == 0;
    }

    /*  Determines the absolute neighbor index that is contained in this
     *  and has the specified position relative to it.  Result is undefined
     *  if "relative" would represent an index in a higher dimension than
     *  this, or if floating() has bitcount of at least 2 and the underlying
     *  dimensionality is 4 or higher.*/
    constexpr NeighborIndex fromRelativeToThis(NeighborIndex relative) {
        assert(relative.i < ipow(3, dimension()));
        if (relative.i == ipow(3, dimension())) {
            // relative is the entire space of *this.
            return *this;
        }
        else if (i == ipow(3, dimension())) {
            // *this is the maximal subspace in its dimension, or is maximal
            // for the first K axes and 0 after that. 
            return relative;
        }
        else if (dimension() == 1) {
            // *this is an edge, so the relative neighbor must be a corner.
            auto newPos = pos();
            if (relative.i == 2) {
                newPos |= floating();
            }
            return fromPosAndFloating(newPos, 0);
        }
        else {
            // This is a face.  We assume our dimensionality is 3; we have
            // ruled out 2 or lower, and if we allow 4 or higher there is no
            // way of determining which dimensions should come next.
            auto axis = highestbit(fixed() & 7);
            auto newPos = pos() | (relative.pos() << (axis + 1));
            newPos |= newPos >> 3;
            newPos &= 7;
            auto newFloat = relative.floating() << (axis + 1);
            newFloat |= newFloat >> 3;
            newFloat &= 7;
            return fromPosAndFloating(newPos, newFloat);
        }
    }

    /*
     *  Builds a NeighborIndex from pos and fixed (described above)
     */
    constexpr static NeighborIndex fromPosAndFloating(uint8_t pos, uint8_t floating) {
        return NeighborIndex((floating || pos)
            ? ((floating & 1) ? 2 : (pos & 1)) +
                fromPosAndFloating(pos >> 1, floating >> 1).i * 3
            : 0);
    }

    int i;

};

/*
 *  A CornerIndex in an N-dimensional space is an N-digit binary value.
 *  Each digit represents an axis in that space:
 *    0:  On the lower edge of that axis
 *    1:  On the upper edge of that axis
 *
 *  For example, the corners in a 2D space are numbered as follows
 *  (left in decimal, right in ternary)
 *
 *    2-------3      10-------11
 *    |       |      |         |
 *    |              |         |
 *    |       |      |         |
 *    0-------1      00-------01
 */
struct CornerIndex {
    constexpr CornerIndex() : i(0) {}
    constexpr CornerIndex(unsigned i) : i(i) {}

    /*  Unpacks from a corner index to a neighbor index */
    constexpr NeighborIndex neighbor() const {
        return NeighborIndex(
            i
            ? 3 * CornerIndex(i >> 1).neighbor().i + (i & 1)
            : 0);
    }
    int i;
};

inline constexpr bool NeighborIndex::contains(const CornerIndex& c) const {
    return i
        ? ((i % 3) == 2 ||
           (i % 3) == (c.i & 1))
            && NeighborIndex(i / 3).contains(CornerIndex(c.i >> 1))
        : (c.i == 0);
}

}   // namespace libfive
