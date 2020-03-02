/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once
#include <vector>
#include <optional>
#include <memory>
#include "libfive/render/brep/simplex/qef.hpp"

namespace libfive {

// Forward declarations
template <unsigned N, class Leaf> class SimplexTree;
template <unsigned N> struct SimplexDCLeaf;
template <unsigned N> struct SimplexLeafSubspace;
class Evaluator;
class Tape;

template <unsigned N>
using SimplexDCTree = SimplexTree<N, SimplexDCLeaf<N>>;

/*  During simplex DC meshing, it is important that any N+1 distinct simplex 
 *  vertices from a "chain" of subspaces (each including the next) not fit in an
 *  N-1-dimensional subspace; in practice, this is done by ensuring that each
 *  vertex is either in the interior of the subspace that generated it, or 
 *  collapsed into the vertex of the subspace that it would otherwise be in.
 *  However, in some cases this collapsing is impossible, either because the
 *  subspace in that direction is the result of merging smaller subspaces, or
 *  because doing so would go past a differently collapsing intermediate 
 *  subspace.  In cases like these, we need to find a new vertex position, and
 *  that is what this class does.  Minimizing the QEF over the entire subspace
 *  is not going to work, as that's what we did the first time, so instead we
 *  start at the center and work outward.  Because in this case we know that
 *  minimizing the QEF over the entire subspace gives a boundary point, we do
 *  not need to worry about missing a feature due to using a portion of the
 *  subspace.  We do, however, still need to worry about a situation in which
 *  the line between our vertex and one of its neighbors intersects the shape's
 *  surface at two points.  Normally, this would be impossible with small enough
 *  boxes because that would mean the line has an extremum of the shape's 
 *  function somewhere in its interior and that extremum should do a better
 *  job of minimizing the QEF than either endpoint, but in this case it may be
 *  that the extremum is outside our reduced box.  Thus, we use derivative 
 *  calculations and a bisection method to find if there is such an extremum
 *  (on the assumption our cell size is small enough that there is at most one 
 *  such), and if there is and it is outside our box, expand our box to include
 *  it and repeat.*/

template <unsigned N, unsigned K>
class SimplexDCRecalculator
{
public:
    /*  Constructor.  'tree' must be a leaf node, and targetSub should reference
     *  a subspace of dimension at least 2.*/
    SimplexDCRecalculator(
        Evaluator* eval, SimplexDCTree<N>* tree, 
        int targetSub, double boundCutoff);

    /*  Actually uses the recalculator to get a new vertex for the subspace.*/
    Eigen::Matrix<double, N, 1> getNewVertex();

protected:
    /*  Gets the extremum on the line between the current vertex and that of 
     *  borderSub.  If there is one, returns true and sets 'result'; otherwise,
     *  returns false and the value of 'result' is undefined.*/
    bool getLineExtremum(SimplexLeafSubspace<N>* borderSub, 
                         Eigen::Matrix<double, K, 1>& result) const;

    /*  Adds to border_subs all subspaces contained in sub number targetSub of 
     *  tree*/
    void addBorderSubs(SimplexDCTree<N>* tree, int targetSub);

protected:
    /*  Compares two subspaces that a higher-dimensional subspace may collapse
     *  to, in order to determine which is the better candidate.  This will be
     *  called only in the case where a high-dimension subspace has a vertex
     *  that should collapse to a lower-dimensional subspace, but there are
     *  multiple intermediate subspaces that do not collapse at all (so it 
     *  should be fairly rare).  Returns true if candidate A is better,
     *  and false if candidate B is better; these should be the neighbor
     *  indices (not relative to toCollapse) of the candidate subspaces to
     *  collapse to.  Result is undefined if the two candidates do not
     *  share at least a corner, or if they are the same.*/
    
    Evaluator* eval;
    Region<K> total_region;
    const std::shared_ptr<Tape>& tape;
    bool is_ambiguous;
    int missing_dim = 0; // Initialized so a bad constructor won't result in a
                         // vector going out-of-bounds.
    double missing_dim_value;
    std::vector<SimplexLeafSubspace<N>*> border_subs;
    QEF<K> summed_qef;
    Eigen::Matrix<double, K, 1> current_vertex;
    double bound_cutoff;
};

}   // namespace libfive
