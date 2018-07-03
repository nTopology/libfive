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
#include "libfive/eval/eval_jacobian.hpp"
#include "libfive/eval/deck.hpp"
#include "libfive/eval/tape.hpp"

namespace Kernel {

JacobianEvaluator::JacobianEvaluator(std::shared_ptr<Deck> d)
    : JacobianEvaluator(d, std::map<Tree::Id, float>())
{
    // Nothing to do here
}

JacobianEvaluator::JacobianEvaluator(
        std::shared_ptr<Deck> d, const std::map<Tree::Id, float>& vars)
    : DerivEvaluator(d, vars),
      j(Eigen::ArrayXXf::Zero(deck->num_clauses + 1, deck->vars.size()))
{
    // Then drop a 1 at each var's position
    size_t index = 0;
    for (auto& v : deck->vars.left)
    {
        j(v.first, index++) = 1;
    }
}

std::map<Tree::Id, float> JacobianEvaluator::gradient(
        const Eigen::Vector3f& p)
{
    return gradient(p, deck->tape);
}
std::map<Tree::Id, float> JacobianEvaluator::gradient(
        const Eigen::Vector3f& p,
        std::shared_ptr<Tape> tape)
{
    // Perform value evaluation, to make sure the f array is correct
    eval(p, tape);

    // Everybody do the tape walk!
    auto ti = tape->rwalk(*this);

    // Unpack from flat array into map
    // (to allow correlating back to VARs in Tree)
    std::map<Tree::Id, float> out;
    size_t index = 0;
    for (auto v : deck->vars.left)
    {
        out[v.second] = j(ti, index++);
    }
    return out;
}

void JacobianEvaluator::operator()(Opcode::Opcode op, Clause::Id id,
                                   Clause::Id a, Clause::Id b)
{
#define ov f(id)
#define oj j.row(id)

#define av f(a)
#define aj j.row(a)

#define bv f(b)
#define bj j.row(b)

        switch (op) {
            case Opcode::OP_ADD:
                oj = aj + bj;
                break;
            case Opcode::OP_MUL:
                oj = av * bj + bv * aj;
                break;
            case Opcode::OP_MIN:
                oj = (av < bv) ? aj : bj;
                break;
            case Opcode::OP_MAX:
                oj = (av < bv) ? bj : aj;
                break;
            case Opcode::OP_SUB:
                oj = aj - bj;
                break;
            case Opcode::OP_DIV:
                oj = (bv*aj - av*bj) / pow(bv, 2);
                break;
            case Opcode::OP_ATAN2:
                oj = (aj*bv - av*bj) / (pow(av, 2) + pow(bv, 2));
                break;
            case Opcode::OP_POW:
                // The full form of the derivative is
                // oj = m * (bv * aj + av * log(av) * bj))
                // However, log(av) is often NaN and bj is always zero,
                // (since it must be CONST), so we skip that part.
                oj = pow(av, bv - 1) * (bv * aj);
                break;
            case Opcode::OP_NTH_ROOT:
                oj = pow(av, 1.0f/bv - 1) * (1.0f/bv * aj);
                break;
            case Opcode::OP_MOD:
                // This isn't quite how partial derivatives of mod work,
                // but close enough normals rendering.
                oj = aj;
                break;
            case Opcode::OP_NANFILL:
                oj = std::isnan(av) ? bj : aj;
                break;
            case Opcode::OP_COMPARE:
                oj.setZero();
                break;

            case Opcode::OP_SQUARE:
                oj = 2 * av * aj;
                break;
            case Opcode::OP_SQRT:
                if (av < 0) oj.setZero();
                else        oj = (aj / (2 * sqrt(av)));
                break;
            case Opcode::OP_NEG:
                oj = -aj;
                break;
            case Opcode::OP_SIN:
                oj = aj * cos(av);
                break;
            case Opcode::OP_COS:
                oj = aj * -sin(av);
                break;
            case Opcode::OP_TAN:
                oj = aj * pow(1/cos(av), 2);
                break;
            case Opcode::OP_ASIN:
                oj = aj / sqrt(1 - pow(av, 2));
                break;
            case Opcode::OP_ACOS:
                oj = aj / -sqrt(1 - pow(av, 2));
                break;
            case Opcode::OP_ATAN:
                oj = aj / (pow(av, 2) + 1);
                break;
            case Opcode::OP_LOG:
                oj = aj / av;
                break;
            case Opcode::OP_EXP:
                oj = exp(av) * aj;
                break;
            case Opcode::OP_ABS:
                oj = (av > 0 ? 1 : -1) * aj;
                break;
            case Opcode::OP_RECIP:
                oj = -aj / pow(av, 2);
                break;

            case Opcode::CONST_VAR:
                oj.setZero();
                break;

            case Opcode::ORACLE:
                oj.setZero();
                break;

            case Opcode::INVALID:
            case Opcode::CONSTANT:
            case Opcode::VAR_X:
            case Opcode::VAR_Y:
            case Opcode::VAR_Z:
            case Opcode::VAR_FREE:
            case Opcode::LAST_OP: assert(false);
        }
#undef ov
#undef oj

#undef av
#undef aj

#undef bv
#undef bj
}

}   // namespace Kernel


