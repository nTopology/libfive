/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <array>
#include "libfive/render/axes.hpp"

namespace libfive {

/*  Forward declarations */
template <unsigned N> class DCTree;
template <unsigned N> class PerThreadBRep;
class Contours;

class DCContourer
{
public:
    using Output = Contours;
    using PerThreadOutput = PerThreadBRep<2>;
    using Input = const DCTree<2>;

    DCContourer(PerThreadBRep<2>& m) : m(m) {}

    /* Empty cell loader */
    void load(Input* input) {}

    template <Axis::Axis A>
    void load(const std::array<Input*, 2>& ts);

    /* Empty corner loader */
    void load(const std::array<Input*, 4> & ts) {}

    static bool needsTopEdges() { return false; }

protected:
    template <Axis::Axis A, bool D>
    void load(const std::array<const DCTree<2>*, 2>& ts);

    PerThreadBRep<2>& m;
};

}   // namespace libfive
