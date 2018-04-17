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
#include "catch.hpp"

#include "libfive/tree/tree.hpp"
#include "util/oracles.hpp"

using namespace Kernel;

TEST_CASE("Joining two trees")
{
    auto t = Tree::X() + 1;

    REQUIRE(t->op == Opcode::OP_ADD);
    REQUIRE(t->lhs->op == Opcode::VAR_X);
    REQUIRE(t->rhs->op == Opcode::CONSTANT);
    REQUIRE(t->rhs->value == 1);
}

TEST_CASE("Tree::flags")
{
    auto a = Tree::var();
    REQUIRE(a->flags == Tree::FLAG_LOCATION_AGNOSTIC);

    auto b = Tree::var();
    REQUIRE(b->flags == Tree::FLAG_LOCATION_AGNOSTIC);

    auto c = a + b;
    REQUIRE(c->flags == Tree::FLAG_LOCATION_AGNOSTIC);

    auto d = c + Tree::X();
    REQUIRE(d->flags == 0);
}

TEST_CASE("Tree::serialize")
{
    SECTION("Basic")
    {
        auto a = min(Tree::X(), Tree::Y());
        auto out = a.serialize();
        std::vector<uint8_t> expected =
            {'T', '"', '"', '"', '"', Opcode::VAR_X, Opcode::VAR_Y, Opcode::OP_MIN, 1, 0, 0, 0, 0, 0, 0, 0, 0xFF};
        REQUIRE(out == expected);
    }

    SECTION("With local references")
    {
        auto a = min(Tree::X(), Tree::X());
        auto out = a.serialize();
        std::vector<uint8_t> expected =
            {'T', '"', '"', '"', '"', Opcode::VAR_X, Opcode::OP_MIN, 0, 0, 0, 0, 0, 0, 0, 0, 0xFF};
        REQUIRE(out == expected);
    }
}

TEST_CASE("Tree::deserialize")
{
    SECTION("Simple")
    {
        auto a = Tree::deserialize(min(Tree::X(), Tree::Y()).serialize());
        REQUIRE(a.id() != nullptr);
        REQUIRE(a->op == Opcode::OP_MIN);
        REQUIRE(a->lhs->op == Opcode::VAR_X);
        REQUIRE(a->rhs->op == Opcode::VAR_Y);
    }

    SECTION("With constant")
    {
        auto a = Tree::deserialize(min(Tree::X(), Tree(2.5f)).serialize());
        REQUIRE(a.id() != nullptr);
        REQUIRE(a->op == Opcode::OP_MIN);
        REQUIRE(a->lhs->op == Opcode::VAR_X);
        REQUIRE(a->rhs->op == Opcode::CONSTANT);
        REQUIRE(a->rhs->value == 2.5f);
    }

    SECTION("With variable")
    {
        auto a = Tree::deserialize(min(Tree::X(), Tree::var()).serialize());
        REQUIRE(a.id() != nullptr);
        REQUIRE(a->op == Opcode::OP_MIN);
        REQUIRE(a->lhs->op == Opcode::VAR_X);
        REQUIRE(a->rhs->op == Opcode::VAR_FREE);
    }
}

TEST_CASE("Tree::remap")
{
    SECTION("Simple")
    {
        auto x = Tree::X();
        auto y = x.remap(Tree::Y(), Tree::X(), Tree::X());
        REQUIRE(y == Tree::Y());
    }

    SECTION("Remapping to a constant")
    {
        auto x = Tree::X();
        auto t = x.remap(Tree(12), Tree::X(), Tree::X());
        REQUIRE(t == Tree(12));
    }

    SECTION("Collapsing while remapping")
    {
        auto x = Tree::X() + 5;
        auto t = x.remap(Tree(3), Tree::X(), Tree::X());
        REQUIRE(t == Tree(8));
    }
}

TEST_CASE("Tree: operator<<")
{
    SECTION("Basic")
    {
        std::stringstream ss;
        ss << (Tree::X() + 5);
        REQUIRE(ss.str() == "(+ x 5)");
    }

    SECTION("With oracle")
    {
        std::stringstream ss;
        auto o = Tree(std::unique_ptr<OracleClause>(new CubeOracleClause));
        ss << (Tree::X() + 5 + o);
        REQUIRE(ss.str() == "(+ x 5 'CubeOracle)");
    }
}

TEST_CASE("Tree::makeVarsConstant")
{
    auto v = Tree::var();
    auto w = Tree::var();
    auto a = 2 * v + 5 * w;
    auto b = a.makeVarsConstant();

    {
        std::stringstream ss;
        ss << a;
        REQUIRE(ss.str() == "(+ (* 2 var-free) (* 5 var-free))");

    }
    {
        std::stringstream ss;
        ss << b;
        REQUIRE(ss.str() == "(+ (* 2 (const-var var-free)) (* 5 (const-var var-free)))");
    }
}
