/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once
#include <map>
#include <Eigen/Eigen>
#include <libfive/render/brep/simplex_dc/simplex_dc_tree.hpp>

namespace libfive {

/*  When performing simplex DC meshing, it is necessary to assign a vertex on
 *  each simplex edge that intersects the surface, and if two such edges meet
 *  at a very small angle and the surface is near (or at) where they meet, those
 *  two edges' intersections with the surface may be identical up to float
 *  precision, or so close as to cause issues.  Thus, this class serves to
 *  merge such intersections (as well as any that are further away but still
 *  identical).
 */

template <unsigned N>
class IntersectionMap
{
public:
    using Vector = Eigen::Matrix<float, N, 1>;

    struct CmpVectors {
        bool operator()(Vector a, Vector b) const
        {
            auto conv = [](Vector v) {
                if constexpr (N == 3) {
                    return std::array{ v.x(), v.y(), v.z() };
                }
                else {
                    return std::array{ v.x(), v.y() };
                }
            };
            return conv(a) < conv(b);
        }
    };

    using IntersectionVec = std::vector<SimplexDCIntersection<N>*>;

    /*  No aligned allocator is needed for our map, since the size of our 
     *  vector is 4 * N bytes, and thus not divisible by 16 (N is 2 or 3).
     */
    using BySub = std::map<Vector, IntersectionVec, CmpVectors>;

    /*  The maps with a key of nullptr represent those intersections not close 
     *  enough to any subspace to get special treatment.
     */
    std::map<const SimplexLeafSubspace<N>*, BySub> m;

    /*  This constructor is needed only for use as a thread-specific output type
     *  for the SimplexDCIntersecter.*/
    IntersectionMap(std::atomic<uint32_t>& c) {}
    /*  And this one for use as an all-threads output type.*/
    IntersectionMap() {}

    void insert(const SimplexLeafSubspace<N>* sub,
                Eigen::Matrix<float, N, 1> pt,
                SimplexDCIntersection<N>* intersection) {
        m[sub][pt].push_back(intersection);
    }

    /*
     *  Collects a set of IntersectionMaps into this one, similar to how
     *  BRep collects PerThreadBReps.  Unlike that collection, this one
     *  is single-threaded, because maps do not accomodate multithreading.
     *  (In the future, if necessary, we may use sorted vectors instead
     *  for the collected form to allow some degree of multithreading.)
     */
    void collect(const std::vector<IntersectionMap>& children)
    {
        m.clear();
        for (const auto& child : children) {
            for (const auto& bySub : child.m) {
                auto& target = m[bySub.first];
                for (const auto& byPt : bySub.second) {
                    target[byPt.first].insert(target[byPt.first].end(),
                        byPt.second.begin(), byPt.second.end());
                }
            }
        }
    }

    /*  Uses the stored map to consolidate intersections that have the same
     *  non-null subspace or a null subspace and the same vertex; they are
     *  then assigned indices (subject to the consolidation) and their vertices
     *  packed into a vector, which is returned.
     */
    std::vector<Eigen::Matrix<float, N, 1>,
        Eigen::aligned_allocator<Eigen::Matrix<float, N, 1>>> buildVecs() const
    {
        uint32_t index(1); // If this is made multithreaded, 
                           // this needs to become atomic.
        std::vector<Eigen::Matrix<float, N, 1>,
            Eigen::aligned_allocator<Eigen::Matrix<float, N, 1>>> out;
        out.push_back(Eigen::Matrix<float, N, 1>::Zero());
        for (const auto& bySub : m) {
            if (bySub.first != nullptr) {
                auto currentIndex = index++;
                auto sumPoint = Eigen::Matrix<float, N, 1>::Zero().eval();
                auto intersectionCount = 0;
                for (const auto& byPt : bySub.second) {
                    sumPoint += (byPt.first * byPt.second.size());
                    intersectionCount += byPt.second.size();
                    for (auto intersection : byPt.second) {
                        intersection->index = currentIndex;
                    }
                }
                assert(intersectionCount > 0);
                out.push_back((sumPoint / intersectionCount).eval());
            }
            else {
                for (const auto& byPt : bySub.second) {
                    auto currentIndex = index++;
                    for (auto intersection : byPt.second) {
                        intersection->index = currentIndex;
                    }
                    out.push_back(byPt.first);
                }
            }
        }
        return out;
    }
};

}   // namespace libfive
