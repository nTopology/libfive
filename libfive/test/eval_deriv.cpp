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
#include <cmath>

#include "catch.hpp"

#include "libfive/tree/tree.hpp"
#include "libfive/eval/eval_deriv.hpp"
#include "libfive/eval/deck.hpp"

using namespace Kernel;

TEST_CASE("DerivEvaluator::deriv")
{
    SECTION("Every operator")
    {
        for (unsigned i=7; i < Kernel::Opcode::ORACLE; ++i)
        {
            auto op = (Kernel::Opcode::Opcode)i;
            Tree t = (Opcode::args(op) == 2 ? Tree(op, Tree::X(), Tree(5))
                                            : Tree(op, Tree::X()));
            auto tape = std::make_shared<Deck>(t);
            DerivEvaluator e(tape);
            e.deriv({0, 0, 0});
            REQUIRE(true /* No crash! */ );
        }
    }

    SECTION("var + 2*X")
    {
        auto v = Tree::var();
        auto t = std::make_shared<Deck>(v + 2 * Tree::X());
        DerivEvaluator e(t, {{v.id(), 0}});

        auto out = e.deriv({2, 0, 0});
        REQUIRE(out == Eigen::Vector4f(2, 0, 0, 4));
    }

    SECTION("X^(1/3)")
    {
        auto t = std::make_shared<Deck>(nth_root(Tree::X(), 3));
        DerivEvaluator e(t);

        {
            auto out = e.deriv({0, 0, 0});
            CAPTURE(out);
            REQUIRE(std::isinf(out(0)));
            REQUIRE(out.bottomRows(3).matrix() == Eigen::Vector3f(0, 0, 0));
        }

        {
            auto out = e.deriv({1, 2, 3});
            CAPTURE(out);
            REQUIRE(out(0) == Approx(0.33333));
            REQUIRE(out.bottomRows(3).matrix() == Eigen::Vector3f(0, 0, 1));
        }
    }
}
