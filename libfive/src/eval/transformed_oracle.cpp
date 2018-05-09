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
#include "libfive/eval/transformed_oracle.hpp"

namespace Kernel {

TransformedOracle::TransformedOracle(
        std::unique_ptr<Oracle> underlying, Tree X_, Tree Y_, Tree Z_)
    : underlying(std::move(underlying)),
      xEvaluator(X_), yEvaluator(Y_), zEvaluator(Z_)
{
    //nothing more to do here.
}

void TransformedOracle::set(const Eigen::Vector3f& p, size_t index)
{
    OracleStorage::set(p, index);
    xEvaluator.array.set(p, index);
    yEvaluator.array.set(p, index);
    zEvaluator.array.set(p, index);
}

void TransformedOracle::evalInterval(Interval::I& out)
{
    auto xRange = xEvaluator.interval.eval(lower, upper);
    auto yRange = yEvaluator.interval.eval(lower, upper);
    auto zRange = zEvaluator.interval.eval(lower, upper);

    Eigen::Vector3f rangeLower{
        xRange.lower(), yRange.lower(), zRange.lower() };
    Eigen::Vector3f rangeUpper{
        xRange.upper(), yRange.upper(), zRange.upper() };

    underlying->set(rangeLower, rangeUpper);
    underlying->evalInterval(out);
}

void TransformedOracle::evalArray(
    Eigen::Block<Eigen::Array<float, Eigen::Dynamic,
    LIBFIVE_EVAL_ARRAY_SIZE, Eigen::RowMajor>, 1, Eigen::Dynamic> out)
{
    setUnderlyingArrayValues(out.cols());
    underlying->evalArray(out);
}

void TransformedOracle::checkAmbiguous(
    Eigen::Block<Eigen::Array<bool, 1, LIBFIVE_EVAL_ARRAY_SIZE>,
    1, Eigen::Dynamic> out)
{
    setUnderlyingArrayValues(out.cols());
    underlying->checkAmbiguous(out);
    out = out || xEvaluator.array.getAmbiguous(out.cols())
              || yEvaluator.array.getAmbiguous(out.cols())
              || zEvaluator.array.getAmbiguous(out.cols());
}

void TransformedOracle::evalDerivs(
    Eigen::Block<Eigen::Array<float, 3, Eigen::Dynamic>,
    3, 1, true> out, size_t index)
{
    Eigen::Matrix3f Jacobian;
    Jacobian << xEvaluator.deriv.deriv(points.col(index)).template head<3>(),
        yEvaluator.deriv.deriv(points.col(index)).template head<3>(),
        zEvaluator.deriv.deriv(points.col(index)).template head<3>();
    Eigen::Vector3f transformedPoint{
        xEvaluator.deriv.eval(points.col(index)),
        yEvaluator.deriv.eval(points.col(index)),
        zEvaluator.deriv.eval(points.col(index))};
    underlying->set(transformedPoint, index);
    underlying->evalDerivs(out, index);
    out = Jacobian * out.matrix();
}

void TransformedOracle::evalDerivArray(
    Eigen::Block<Eigen::Array<float, 3, LIBFIVE_EVAL_ARRAY_SIZE>,
    3, Eigen::Dynamic, true> out)
{
    TransformedOracle::setUnderlyingArrayValues(out.cols());

    auto xDerivs = xEvaluator.array.derivs(out.cols());
    auto yDerivs = yEvaluator.array.derivs(out.cols());
    auto zDerivs = zEvaluator.array.derivs(out.cols());

    underlying->evalDerivArray(out);

    for (auto i = 0; i < out.cols(); ++i)
    {
        Eigen::Matrix3f Jacobian;
        Jacobian << xDerivs.col(i).template head<3>(), 
                    yDerivs.col(i).template head<3>(),
                    zDerivs.col(i).template head<3>();
        out.col(i) = Jacobian * out.col(i).matrix();
    }
}

void TransformedOracle::evalFeatures(
    boost::container::small_vector<Feature, 4>& out)
{
    out.clear();
    auto pt = points.col(0);
    Eigen::Vector3f transformedPoint{
        xEvaluator.feature.eval(pt),
        yEvaluator.feature.eval(pt),
        zEvaluator.feature.eval(pt) };

    auto xFeatures = xEvaluator.feature.features_(pt);
    auto yFeatures = yEvaluator.feature.features_(pt);
    auto zFeatures = zEvaluator.feature.features_(pt);

    boost::container::small_vector<Feature, 4> underlyingOut;
    underlying->set(transformedPoint);
    underlying->evalFeatures(underlyingOut);

    /*  This O(n^4) loop should almost never have large values for all input
     *  sizes, barring intentionally pathological cases.
     */
    for (auto f1 : xFeatures)
    {
        for (auto f2 : yFeatures)
        {
            if (f1.check(f2))
            {
                Feature f12({ 0.f, 0.f, 0.f }, f1, f2);
                for (auto f3 : zFeatures)
                {
                    if (f3.check(f12))
                    {
                        Feature f123({ 0.f, 0.f, 0.f }, f12, f3);
                        Eigen::Matrix3f Jacobian;
                        Jacobian << f1.deriv, f2.deriv, f3.deriv;
                        for (auto f4 : underlyingOut)
                        {
                            Feature transformed(f4, Jacobian);
                            if (transformed.check(f123))
                            {
                                out.emplace_back(
                                    transformed.deriv, transformed, f123);
                            }
                        }
                    }
                }
            }
        }
    }
}

void TransformedOracle::evalAndPushInterval(Interval::I& out)
{
    auto xRange = xEvaluator.interval.evalAndPush(lower, upper);
    auto yRange = yEvaluator.interval.evalAndPush(lower, upper);
    auto zRange = zEvaluator.interval.evalAndPush(lower, upper);

    Eigen::Vector3f rangeLower{
        xRange.lower(), yRange.lower(), zRange.lower() };
    Eigen::Vector3f rangeUpper{
        xRange.upper(), yRange.lower(), zRange.lower() };

    underlying->set(rangeLower, rangeUpper);
    underlying->evalAndPushInterval(out);
}

void TransformedOracle::evalPoint(float& out, size_t index)
{
    Eigen::Vector3f transformedPoint{
        xEvaluator.feature.eval(points.col(index)),
        yEvaluator.feature.eval(points.col(index)),
        zEvaluator.feature.eval(points.col(index)) };

    underlying->set(transformedPoint, index);
    underlying->evalPoint(out, index);
}

void TransformedOracle::evalAndPushPoint(float& out, size_t index)
{
    Eigen::Vector3f transformedPoint{
        xEvaluator.feature.evalAndPush(points.col(index)),
        yEvaluator.feature.evalAndPush(points.col(index)),
        zEvaluator.feature.evalAndPush(points.col(index)) };

    underlying->set(transformedPoint, index);
    underlying->evalAndPushPoint(out, index);
}

void TransformedOracle::baseEvalPoint(float& out, size_t index)
{
    Eigen::Vector3f transformedPoint{
        xEvaluator.feature.baseEval(points.col(index)),
        yEvaluator.feature.baseEval(points.col(index)),
        zEvaluator.feature.baseEval(points.col(index)) };

    underlying->set(transformedPoint, index);
    underlying->baseEvalPoint(out, index);
}

void TransformedOracle::baseEvalDerivs(
    Eigen::Block<Eigen::Array<float, 3, Eigen::Dynamic>,
    3, 1, true> out, size_t index) {
    Eigen::Matrix3f Jacobian;
    Jacobian << xEvaluator.baseDeriv.baseDeriv(points.col(index)),
        yEvaluator.baseDeriv.baseDeriv(points.col(index)),
        zEvaluator.baseDeriv.baseDeriv(points.col(index));
    Eigen::Vector3f transformedPoint{
        xEvaluator.deriv.baseEval(points.col(index)),
        yEvaluator.deriv.baseEval(points.col(index)),
        zEvaluator.deriv.baseEval(points.col(index))};
    underlying->set(transformedPoint, index);
    underlying->baseEvalDerivs(out, index);
    out = Jacobian * out.matrix();
}

void TransformedOracle::pop() {
  xEvaluator.deriv.pop();
  yEvaluator.deriv.pop();
  zEvaluator.deriv.pop();
  underlying->pop();
}

void TransformedOracle::setUnderlyingArrayValues(int count)
{
    auto xPoints = xEvaluator.array.values(count);
    auto yPoints = yEvaluator.array.values(count);
    auto zPoints = zEvaluator.array.values(count);
    for (auto i = 0; i < count; ++i)
    {
        underlying->set({ xPoints(i), yPoints(i), zPoints(i) }, i);
    }
}


} //Namespace Kernel

