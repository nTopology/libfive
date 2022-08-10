/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2020  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "libfive/tree/data.hpp"
#include "libfive/oracle/oracle.hpp"

namespace libfive {

Opcode::Opcode TreeData::op() const {
    if (auto i = std::get_if<TreeNonaryOp>(this)) {
        return i->op;
    } else if (auto i = std::get_if<TreeUnaryOp>(this)) {
        return i->op;
    } else if (auto i = std::get_if<TreeBinaryOp>(this)) {
        return i->op;
    } else if (std::get_if<TreeConstant>(this)) {
        return Opcode::CONSTANT;
    } else if (std::get_if<TreeOracle>(this)) {
        return Opcode::ORACLE;
    } else if (std::get_if<TreeInvalid>(this)) {
        return Opcode::INVALID;
    } else if (std::get_if<TreeRemap>(this)) {
        return Opcode::INVALID;
    } else if (std::get_if<TreeApply>(this)) {
        return Opcode::INVALID;
    } else {
        assert(false); // All variants should be covered above
        return Opcode::INVALID;
    }
}

const Tree& TreeData::lhs() const {
    if (auto i = std::get_if<TreeUnaryOp>(this)) {
        return i->lhs;
    } else if (auto i = std::get_if<TreeBinaryOp>(this)) {
        return i->lhs;
    } else {
        throw ChildException();
    }
}

const Tree& TreeData::rhs() const {
    if (auto i = std::get_if<TreeBinaryOp>(this)) {
        return i->rhs;
    } else {
        throw ChildException();
    }
}

float TreeData::value() const {
    // Can't use std::get<TreeConstant> because it requires a newer macOS
    // than my main development machine.
    if (auto i = std::get_if<TreeConstant>(this)) {
        return i->value;
    } else {
        throw ValueException();
    }
}

std::unique_ptr<Oracle> TreeData::build_oracle() const {
    if (auto i = std::get_if<TreeOracle>(this)) {
        return i->oracle->getOracle();
    } else {
        throw OracleException();
    }
}

const OracleClause& TreeData::oracle_clause() const {
    if (auto i = std::get_if<TreeOracle>(this)) {
        return *(i->oracle);
    } else {
        throw OracleException();
    }
}

TreeData::Key TreeData::key() const {
    if (auto d = std::get_if<TreeConstant>(this)) {
        if (std::isnan(d->value)) {
            return Key(true);
        } else {
            return Key(d->value);
        }
    } else if (auto d = std::get_if<TreeNonaryOp>(this)) {
        if (d->op == Opcode::VAR_FREE) {
            return Key(std::make_tuple(d->op, this));
        } else {
            return Key(d->op);
        }
    } else if (auto d = std::get_if<TreeUnaryOp>(this)) {
        return Key(std::make_tuple(d->op, d->lhs.get()));
    } else if (auto d = std::get_if<TreeBinaryOp>(this)) {
        return Key(std::make_tuple(d->op, d->lhs.get(), d->rhs.get()));
    } else if (auto d = std::get_if<TreeOracle>(this)) {
        return Key(d->oracle.get());
    } else if (auto d = std::get_if<TreeRemap>(this)) {
        return Key(std::make_tuple(d->x.get(), d->y.get(),
                                   d->z.get(), d->t.get()));
    } else if (auto d = std::get_if<TreeApply>(this)) {
        return Key(std::make_tuple(
                    d->target.get(), d->value.get(), d->t.get()));
    } else {
        return Key(false);
    }
}

uint32_t TreeData::compute_flags() const {
    if (std::get_if<TreeConstant>(this)) {
        return 0;
    } else if (auto d = std::get_if<TreeNonaryOp>(this)) {
        return (d->op == Opcode::VAR_X ||
                d->op == Opcode::VAR_Y ||
                d->op == Opcode::VAR_Z)
            ? TREE_FLAG_HAS_XYZ
            : 0;
    } else if (auto d = std::get_if<TreeUnaryOp>(this)) {
        return d->lhs->flags;
    } else if (auto d = std::get_if<TreeBinaryOp>(this)) {
        return d->lhs->flags | d->rhs->flags;
    } else if (std::get_if<TreeOracle>(this)) {
        return TREE_FLAG_HAS_ORACLE;
    } else if (auto d = std::get_if<TreeRemap>(this)) {
        return d->x->flags | d->y->flags | d->z->flags |
               d->t->flags | TREE_FLAG_HAS_REMAP;
    } else if (auto d = std::get_if<TreeApply>(this)) {
        return d->target->flags | d->value->flags |
               d->t->flags | TREE_FLAG_HAS_REMAP;
    } else {
        return 0;
    }
}

}   // namespace libfive
