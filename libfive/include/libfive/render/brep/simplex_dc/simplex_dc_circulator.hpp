/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"

namespace libfive {

/*  When performing simplex DC meshing, it is necessary at a couple of points
 *  to traverse all of the simplices around a given simplex edge.  This can be
 *  complicated when there are degenerate simplices involved (the degenerate
 *  simplices themselves can be ignored, but the resulting duplicate vertex can
 *  cause the traversal to follow complicated paths with respect to the cells
 *  of the tree), so a SimplexDCCirculator class is used to function as an
 *  iterator and handle the process that way.
 *  In 2 dimensions, the "circulator" actually only has two points, but using
 *  the same class as 3d makes the SimplexDCIntersecter easier to handle.
 *  Regardless of dimension, it should not be created around a pair of cell
 *  features (i.e. two of cell/face/edge/corner) that is not minimal-dimension
 *  for features with those particular (post-collapse) subspace vertices.*/

template <unsigned N, bool constTree>
class SimplexDCCirculator;

template <>
class SimplexDCCirculator<3, true>
{
public:

    /// Constructs the circulator by specifying a cell, and the inside
    /// and outside subspaces.  It is possible, due to cell merging,
    /// that the inside or outside is a cell but the cell passed must
    /// be a different one; in such case, the sub corresponding to the
    /// cell will be from 27 to 52, being 27 plus the neighbor position
    /// of the cell.
    SimplexDCCirculator(const SimplexDCTree<3>* cell,
        unsigned insideIdx, unsigned outsideIdx)
    {
        auto insideDim = NeighborIndex(insideIdx).dimension();
        auto outsideDim = NeighborIndex(insideIdx).dimension();
        assert(NeighborIndex(insideIdx).dimension() != 
               NeighborIndex(insideIdx).dimension());
        construct(cell, { insideIdx, outsideIdx });
        insideDimension = insideDim;
        outsideDimension = outsideDim;
        while (!isValid() && !hitError) {
            advance();
        }
        startSimplex = simplex();
    }

    SimplexDCCirculator& operator++() {
        do { advance(); } while (!isValid() && !hitError);
        if (simplex() == startSimplex) {
            didCycle = true;
        }
        return *this;
    }

    const DCSimplex<3>* simplex() const {
        return mutableSimplex();
    }

    bool isEnd() { return didCycle || hitError; }

protected:
    /// Helper method for constructor; does not set insideDimension or 
    /// outsideDimension, and thus the order of indexPr does not matter.
    void construct(const SimplexDCTree<3>* cell,
                   std::array<unsigned, 2> indexPr);

    void advance();

    bool isValid() const;

    DCSimplex<3>* mutableSimplex() const {
        if (hitError) {
            return nullptr;
        }
        const auto& edge = edgeCell->leaf->edgeFromReduced(
            edgeAxis, edgeRelativePos);
        assert(std::holds_alternative<SimplexDCMinEdge<3>*>(edge));
        return &std::get<SimplexDCMinEdge<3>*>(edge)->
            simplexWithEdgeReducedCell(
                edgeAxis, faceAxis, relativeCellPos, cornerPos);
    }

    /*  Gets the K-dimensional subspace associated with this simplex.*/
    template <unsigned K>
    SimplexLeafSubspace<3>* getSub() const;

    // We'll want to specialize our getSub's, so declare specializations
    // before using them.
    template <> SimplexLeafSubspace<3>* getSub<0>() const;
    template <> SimplexLeafSubspace<3>* getSub<1>() const;
    template <> SimplexLeafSubspace<3>* getSub<2>() const;
    template <> SimplexLeafSubspace<3>* getSub<3>() const;

    // getSub with runtime-determined dimension.
    SimplexLeafSubspace<3>* getSub(unsigned k) const {
        switch (k) {
        case 0:
            return getSub<0>();
        case 1:
            return getSub<1>();
        case 2:
            return getSub<2>();
        case 3:
            return getSub<3>();
        default:
            assert(false);
            return nullptr;
        }
    }


    unsigned insideDimension;
    unsigned outsideDimension;

    const SimplexDCTree<3>* edgeCell;

    /*  These could be packed a lot tighter, but circulators are short-lived, 
     *  so there's no need for efficient packing.*/
    Axis::Axis edgeAxis;
    int edgeRelativePos;
    Axis::Axis faceAxis;
    bool cornerPos;
    int relativeCellPos; // Relative to the edge axis.

    bool didCycle = false;

    bool hitError = false;

    const DCSimplex<3>* startSimplex;
};

template <>
class SimplexDCCirculator<2, true>
{
public:
    /// Constructs the circulator by specifying a cell, and the inside
    /// and outside subspaces.  It is possible, due to cell merging,
    /// that the inside or outside is a cell but the cell passed must
    /// be a different one; in such case, the sub corresponding to the
    /// cell will be from 9 to 16, being 9 plus the neighbor position
    /// of the cell.
    SimplexDCCirculator(const SimplexDCTree<2>* cell,
                        unsigned insideIdx, unsigned outsideIdx)
        : simplices{backwardSimplex(cell, insideIdx, outsideIdx),
                    forwardSimplex(cell, insideIdx, outsideIdx) }
    {}

    SimplexDCCirculator& operator++() {
        ++index;
        return *this;
    }

    const DCSimplex<2>* simplex() const {
        return mutableSimplex();
    }

    bool isEnd() { return index > 1; }

protected:
    DCSimplex<2>* mutableSimplex() const {
        return simplices[index % 2];
    }

    static DCSimplex<2>* forwardSimplex(
        const SimplexDCTree<2>* cell, unsigned insideIdx, unsigned outsideIdx);

    static DCSimplex<2>* backwardSimplex(
        const SimplexDCTree<2>* cell, unsigned insideIdx, unsigned outsideIdx) 
    { return forwardSimplex(cell, outsideIdx, insideIdx); }


    std::array<DCSimplex<2>*, 2> simplices;
    unsigned index;
};

template <unsigned N>
class SimplexDCCirculator<N, false> : SimplexDCCirculator<N, true>
{
public:

    template <typename ... Ts>
    SimplexDCCirculator(Ts... ts)
        : SimplexDCCirculator<N, true>(ts...) {}

    SimplexDCCirculator& operator++() {
        SimplexDCCirculator<N, true>::operator++();
        return *this;
    }

    DCSimplex<N>* simplex() const {
        return mutableSimplex();
    }
};

}   // namespace libfive
