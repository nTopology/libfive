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
#include "libfive/eval/eval_interval.hpp"
#include "libfive/eval/deck.hpp"
#include "libfive/eval/tape.hpp"

namespace Kernel {

IntervalEvaluator::IntervalEvaluator(std::shared_ptr<Deck> d)
    : IntervalEvaluator(d, std::map<Tree::Id, float>())
{
    // Nothing to do here
}

IntervalEvaluator::IntervalEvaluator(
        std::shared_ptr<Deck> d, const std::map<Tree::Id, float>& vars)
    : BaseEvaluator(d, vars)
{
    i.resize(d->num_clauses + 1);

    // Unpack variables into result array
    for (auto& v : d->vars.right)
    {
        auto var = vars.find(v.first);
        i[v.second] = (var != vars.end()) ? var->second : 0;
    }

    // Unpack constants into result array
    for (auto& c : d->constants)
    {
        i[c.first] = c.second;
    }
}

Interval::I IntervalEvaluator::eval(const Eigen::Vector3f& lower,
                                    const Eigen::Vector3f& upper)
{
    return eval(lower, upper, deck->tape);
}

Interval::I IntervalEvaluator::eval(const Eigen::Vector3f& lower,
                                    const Eigen::Vector3f& upper,
                                    Tape::Handle tape)
{
    i[deck->X] = {lower.x(), upper.x()};
    i[deck->Y] = {lower.y(), upper.y()};
    i[deck->Z] = {lower.z(), upper.z()};

    for (auto& o : deck->oracles)
    {
        o->set(lower, upper);
    }

    // Mark initial safety based on the starting values
    safe = true;
    for (auto a : {deck->X, deck->Y, deck->Z})
    {
        const auto lower = i[a].lower();
        const auto upper = i[a].lower();
        if (std::isnan(lower) || std::isinf(lower) ||
            std::isnan(upper) || std::isinf(upper))
        {
            safe = false;
        }
    }

    return i[tape->rwalk(*this)];
}

std::pair<Interval::I, Tape::Handle> IntervalEvaluator::evalAndPush(
        const Eigen::Vector3f& lower,
        const Eigen::Vector3f& upper)
{
    return evalAndPush(lower, upper, deck->tape);
}

std::pair<Interval::I, Tape::Handle> IntervalEvaluator::evalAndPush(
        const Eigen::Vector3f& lower,
        const Eigen::Vector3f& upper,
        Tape::Handle tape)
{
    auto out = eval(lower, upper, tape);

    auto p = Tape::push(tape, *deck,
        [&](Opcode::Opcode op, Clause::Id /* id */,
            Clause::Id a, Clause::Id b)
    {
        // For min and max operations, we may only need to keep one branch
        // active if it is decisively above or below the other branch.
        if (op == Opcode::OP_MAX)
        {
            if (a == b)
            {
                return Tape::KEEP_A;
            }
            else if (i[a].lower() > i[b].upper() || a == b)
            {
                return Tape::KEEP_A;
            }
            else if (i[b].lower() > i[a].upper())
            {
                return Tape::KEEP_B;
            }
            else
            {
                return Tape::KEEP_BOTH;
            }
        }
        else if (op == Opcode::OP_MIN)
        {
            if (a == b)
            {
                return Tape::KEEP_A;
            }
            else if (i[a].lower() > i[b].upper())
            {
                return Tape::KEEP_B;
            }
            else if (i[b].lower() > i[a].upper())
            {
                return Tape::KEEP_A;
            }
            else
            {
                return Tape::KEEP_BOTH;
            }
        }
        return Tape::KEEP_ALWAYS;
    },
        Tape::INTERVAL,
        {{i[deck->X].lower(), i[deck->Y].lower(), i[deck->Z].lower()},
         {i[deck->X].upper(), i[deck->Y].upper(), i[deck->Z].upper()}});
    return std::make_pair(out, std::move(p));
}

////////////////////////////////////////////////////////////////////////////////

bool IntervalEvaluator::setVar(Tree::Id var, float value)
{
    auto v = deck->vars.right.find(var);
    if (v != deck->vars.right.end())
    {
        bool changed = i[v->second] != value;
        i[v->second] = value;
        return changed;
    }
    else
    {
        return false;
    }
}

////////////////////////////////////////////////////////////////////////////////

void IntervalEvaluator::operator()(Opcode::Opcode op, Clause::Id id,
                                   Clause::Id a_, Clause::Id b_)
{
#define out i[id]
#define a i[a_]
#define b i[b_]
    switch (op) {
        case Opcode::OP_ADD:
            out = a + b;
            break;
        case Opcode::OP_MUL:
            out = a * b;
            break;
        case Opcode::OP_MIN:
            out = boost::numeric::min(a, b);
            break;
        case Opcode::OP_MAX:
            out = boost::numeric::max(a, b);
            break;
        case Opcode::OP_SUB:
            out = a - b;
            break;
        case Opcode::OP_DIV:
            out = a / b;
            break;
        case Opcode::OP_ATAN2:
            out = atan2(a, b);
            break;
        case Opcode::OP_POW:
            out = boost::numeric::pow(a, b.lower());
            break;
        case Opcode::OP_NTH_ROOT:
            out = boost::numeric::nth_root(a, b.lower());
            break;
        case Opcode::OP_MOD:
            out = Interval::I(0.0f, b.upper()); // YOLO
            break;
        case Opcode::OP_NANFILL:
            out = (std::isnan(a.lower()) || std::isnan(a.upper())) ? b : a;
            break;
        case Opcode::OP_COMPARE:
            if      (a.upper() < b.lower()) out = Interval::I(-1, -1);
            else if (a.lower() > b.upper()) out = Interval::I( 1,  1);
            else                            out = Interval::I(-1,  1);
            break;

        case Opcode::OP_SQUARE:
            out = boost::numeric::square(a);
            break;
        case Opcode::OP_SQRT:
            out = boost::numeric::sqrt(a);
            break;
        case Opcode::OP_NEG:
            out = -a;
            break;
        case Opcode::OP_SIN:
            out = boost::numeric::sin(a);
            break;
        case Opcode::OP_COS:
            out = boost::numeric::cos(a);
            break;
        case Opcode::OP_TAN:
            out = boost::numeric::tan(a);
            break;
        case Opcode::OP_ASIN:
            out = boost::numeric::asin(a);
            break;
        case Opcode::OP_ACOS:
            out = boost::numeric::acos(a);
            break;
        case Opcode::OP_ATAN:
            // If the interval has an infinite bound, then return the largest
            // possible output interval (of +/- pi/2).  This rescues us from
            // situations where we do atan(y / x) and the behavior of the
            // interval changes if you're approaching x = 0 from below versus
            // from above.
            out = (std::isinf(a.lower()) || std::isinf(a.upper()))
                ? Interval::I(-M_PI/2, M_PI/2)
                : boost::numeric::atan(a);
            break;
        case Opcode::OP_EXP:
            out = boost::numeric::exp(a);
            break;
        case Opcode::OP_LOG:
            out = boost::numeric::log(a);
            break;
        case Opcode::OP_ABS:
            out = boost::numeric::abs(a);
            break;
        case Opcode::OP_RECIP:
            out = Interval::I(1,1) / a;
            break;

        case Opcode::CONST_VAR:
            out = a;
            break;

        case Opcode::ORACLE:
            deck->oracles[a_]->evalInterval(out);
            break;

        case Opcode::INVALID:
        case Opcode::CONSTANT:
        case Opcode::VAR_X:
        case Opcode::VAR_Y:
        case Opcode::VAR_Z:
        case Opcode::VAR_FREE:
        case Opcode::LAST_OP: assert(false);
    }

    // Track the safety of the evaluation
    if (safe)
    {
        const auto lower = out.lower();
        const auto upper = out.upper();
        if (std::isnan(lower) || std::isinf(lower) ||
            std::isnan(upper) || std::isinf(upper))
        {
            safe = false;
        }
    }

#undef out
#undef a
#undef b
}

}   // namespace Kernel
