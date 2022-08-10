/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "catch.hpp"

#include "libfive/render/brep/mesh.hpp"
#include "libfive/render/brep/region.hpp"
#include "libfive/oracle/oracle_clause.hpp"
#include "libfive/oracle/oracle_storage.hpp"

#include "libfive/eval/deck.hpp"
#include "libfive/eval/evaluator.hpp"

#include "util/shapes.hpp"
#include "util/oracles.hpp"

using namespace libfive;

/*  In order to test the transformed primitives, we can't test that the meshing
 *  is identical when they're used, since numeric error can cause them not to
 *  be.  So instead, the evaluators are tested to ensure they return the same
 *  results, up to numeric error.
 */

void compareUnderTransformation(Tree oracleTree, Tree controlTree,
    std::function<Tree(Tree)> transformation,
    std::vector<Eigen::Vector3f> testPoints)
{
    auto transformedOracle = transformation(oracleTree);
    auto transformedControl = transformation(controlTree);
    auto oDeck = std::make_shared<Deck>(transformedOracle);
    auto cDeck = std::make_shared<Deck>(transformedControl);
    /* If a point is ambiguous in either tree, we can't compare
     * individual derivatives, since there is more than one valid
     * result and no guarantee both trees will give the same one.
     */
    Eigen::Array<bool, 1, 256> ambigPoints(testPoints.size());
    {
        DerivArrayEvaluator o(oDeck);
        DerivArrayEvaluator c(cDeck);
        for (unsigned i = 0; i < testPoints.size(); ++i)
        {
            o.set(testPoints[i], i);
            c.set(testPoints[i], i);
        }
        /*  getAmbiguous is numerically unstable, so unfortunately it
         *  cannot be tested.  We still calculate it, though, since
         *  when either is ambiguous that means the derivatives cannot be
         *  tested either except via features, since there is more than one
         *  possible correct answer.
         */
        auto oAmbig = o.getAmbiguous(testPoints.size());
        auto cAmbig = c.getAmbiguous(testPoints.size());
        ambigPoints.head(testPoints.size()) = oAmbig && cAmbig;
        auto oResults = o.derivs(testPoints.size());
        auto cResults = c.derivs(testPoints.size());
        for (unsigned i = 0; i < testPoints.size(); ++i)
        {
            CAPTURE(i);
            CAPTURE(oResults.col(i));
            CAPTURE(cResults.col(i));
            Eigen::Vector4f diff = oResults.col(i) - cResults.col(i);
            auto relativeClose = oResults.col(i).isApprox(cResults.col(i));
            if (!ambigPoints(i))
            {
                if (!relativeClose)
                {
                    REQUIRE(diff.norm() < 1e-4);
                }
                else
                {
                    REQUIRE(true); /*That's a successful test too.*/
                }
            }
        }
    }
    {
        FeatureEvaluator o(oDeck);
        FeatureEvaluator c(cDeck);

        for (auto point : testPoints)
        {
            auto oFeatures = o.features(point);
            auto cFeatures = c.features(point);
            CAPTURE(point);
            CAPTURE(oFeatures.size());
            CAPTURE(cFeatures.size());
            /*  We would like to capture which feature is different, but the
             *  fact that they may be permuted makes it more convenient
             *  to add any necessary statements only as needed.
             */

            auto valid = std::is_permutation(oFeatures.begin(), oFeatures.end(),
                cFeatures.begin(), [](
                    const Eigen::Vector3f& first,
                    const Eigen::Vector3f& second)
            {return (first - second).norm() < (1e-4) ||
                first.isApprox(second, 1e-4);});
            REQUIRE(valid);
        }
        /*  isInside is numerically unstable, likely in a way that does
         *  preclude testing it in any meaningful way.
         */
    }
    {
        DerivArrayEvaluator o(oDeck);
        DerivArrayEvaluator c(cDeck);
        for (unsigned i = 0; i < testPoints.size(); ++i)
        {
            auto oDeriv = o.deriv(testPoints[i]);
            auto cDeriv = c.deriv(testPoints[i]);
            CAPTURE(testPoints[i]);
            CAPTURE(cDeriv);
            CAPTURE(oDeriv);
            if (!ambigPoints(i))
            {
                if (oDeriv.isApprox(cDeriv, 1e-4))
                {
                    REQUIRE(true);
                }
                else
                {
                    REQUIRE((oDeriv - cDeriv).norm() < 1e-4);
                }
            }
        }
        /*  IntervalEvaluator cannot be comparison-tested; while both
         *  approaches should return results that contain the actual
         *  interval, neither is required to be exact.
         */
    }
}

TEST_CASE("Rotated Oracle: Render and compare (cube)")
{
    auto cube = max(max(
        max(-(Tree::X() + 1.5),
            Tree::X() - 1.5),
        max(-(Tree::Y() + 1.5),
            Tree::Y() - 1.5)),
        max(-(Tree::Z() + 1.5),
            Tree::Z() - 1.5));
    Tree cubeOracle = convertToOracleAxes(cube);
    compareUnderTransformation(cubeOracle, cube,
        [](Tree t) {return rotate2d(t, 10);},
        { {1.5, 1.5, 1.5},{ -1.5, -1.5, 1.5 },
        { -1.5, 1.5, -1.5 },{ 1.5, -1.5, -1.5 }, {0., 0., 0.} });
    /*  Numeric error means that with a rotation we can't really hit where
     *  interesting stuff is going to happen anyway, except for 0.
     */
}

TEST_CASE("Rotated Oracle: Render and compare (cube as oracle)")
{
    auto cube = max(max(
        max(-(Tree::X() + 1.5),
            Tree::X() - 1.5),
        max(-(Tree::Y() + 1.5),
            Tree::Y() - 1.5)),
        max(-(Tree::Z() + 1.5),
            Tree::Z() - 1.5));
    Tree cubeOracle(std::make_unique<CubeOracleClause>());

    compareUnderTransformation(cubeOracle, cube,
        [](Tree t) {return rotate2d(t, 10);},
        { { 1.5, 1.5, 1.5 },{ -1.5, -1.5, 1.5 },
        { -1.5, 1.5, -1.5 },{ 1.5, -1.5, -1.5 },
        { 0., 0., 0. } });
}

TEST_CASE("Abs and skew applied to Oracle: "
    "Render and compare (cube)")
{
    auto cube = max(max(
        max(-(Tree::X() + 1.5),
            Tree::X() - 1.5),
        max(-(Tree::Y() + 1.5),
            Tree::Y() - 1.5)),
        max(-(Tree::Z() + 1.5),
            Tree::Z() - 1.5));
    Tree cubeOracle = convertToOracleAxes(cube);
    compareUnderTransformation(cubeOracle, cube,
        [](Tree t) {
        return t.remap(Tree::Y(), Tree::X(),
            abs(Tree::Z() + Tree::X() * 0.2f));},
        { { 1.5f, 1.5f, 1.8f },{ -1.5f, -1.5f, 1.8f},
        { -1.5f, 1.5f, -1.8f },{ 1.5f, -1.5f, -1.8f }, { 0.f, 0.f, 0.f } });
        //  Some of these points are maps of corners, some are not.

}

TEST_CASE("Abs and skew applied to Oracle: "
          "Render and compare (cube as oracle)")
{
    auto cube = max(max(
        max(-(Tree::X() + 1.5),
            Tree::X() - 1.5),
        max(-(Tree::Y() + 1.5),
            Tree::Y() - 1.5)),
        max(-(Tree::Z() + 1.5),
            Tree::Z() - 1.5));
    Tree cubeOracle(std::make_unique<CubeOracleClause>());
    compareUnderTransformation(cubeOracle, cube,
        [](Tree t) {
        return t.remap(Tree::Y(), Tree::X(),
            abs(Tree::Z() + Tree::X() * 0.2f));},
        { { 1.5f, 1.5f, 1.8f },{ -1.5f, -1.5f, 1.8f },
        { -1.5f, 1.5f, -1.8f },{ 1.5f, -1.5f, -1.8f },{ 0.f, 0.f, 0.f } });
    //  Some of these points correspond to corners, some do not.
}

TEST_CASE("Jacobian-0 transform and abs applied to Oracle: "
    "Render and compare (cube)")
{
    auto cube = max(max(
        max(-(Tree::X() + 1.5),
            Tree::X() - 1.5),
        max(-(Tree::Y() + 1.5),
            Tree::Y() - 1.5)),
        max(-(Tree::Z() + 1.5),
            Tree::Z() - 1.5));
    Tree cubeOracle = convertToOracleAxes(cube);

    compareUnderTransformation(cubeOracle, cube,
        [](Tree t) {
        return t.remap(Tree::X() * Tree::Y(), Tree::Y() * abs(Tree::Z()),
            abs(Tree::Z()) * Tree::X());},
        { { 1.5, 1., 1. },{ -1., -1.5, 1. },
        { -1., 1.5, -1. },{ 1., -1., -1. },{ 0., 0., 0. } });
    //  None of these correspond to corners

   cube = max(max(
        max(-(Tree::X() + 1.),
            Tree::X() - 1.),
        max(-(Tree::Y() + 1.),
            Tree::Y() - 1.)),
        max(-(Tree::Z() + 1.),
            Tree::Z() - 1.));
    cubeOracle = convertToOracleAxes(cube);

    compareUnderTransformation(cubeOracle, cube,
        [](Tree t) {
        return t.remap(Tree::X() * Tree::Y(), Tree::Y() * abs(Tree::Z()),
            abs(Tree::Z()) * Tree::X());},
        { { 1., 1., 1. },{ -1., -1., 1. },
        { -1., 1., -1. },{ 1., -1., -1. },{ 0., 0., 0. } });
    //  But these do
}

TEST_CASE("Jacobian-0 transform and abs applied to Oracle: "
    "Render and compare (cube as oracle)")
{
    auto cube = max(max(
        max(-(Tree::X() + 1.5),
            Tree::X() - 1.5),
        max(-(Tree::Y() + 1.5),
            Tree::Y() - 1.5)),
        max(-(Tree::Z() + 1.5),
            Tree::Z() - 1.5));
    Tree cubeOracle(std::make_unique<CubeOracleClause>());

    compareUnderTransformation(cubeOracle, cube,
        [](Tree t) {
        return t.remap(Tree::X() * Tree::Y(), Tree::Y() * abs(Tree::Z()),
            abs(Tree::Z()) * Tree::X());},
        { { 1.5, 1., 1. },{ -1., -1.5, 1. },
        { -1., 1.5, -1. },{ 1., -1., -1. },{ 0., 0., 0. } });
    //  None of these correspond to corners
}
