/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2018  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include "libfive/render/brep/simplex_dc/simplex_dc_tree.hpp"
#include "../object_pool.inl"

namespace libfive {
template class ObjectPool<SimplexTree<2, SimplexDCLeaf<2>>, SimplexDCLeaf<2>, 
                          SimplexLeafSubspace<2>, SimplexDCMinEdge<2>, 
                          SimplexDCIntersection<2>>;
template class ObjectPool<SimplexDCLeaf<2>, SimplexLeafSubspace<2>, 
                          SimplexDCMinEdge<2>, SimplexDCIntersection<2>>;
template class ObjectPool<SimplexLeafSubspace<2>, SimplexDCMinEdge<2>,
                          SimplexDCIntersection<2>>;
template class ObjectPool<SimplexDCMinEdge<2>, SimplexDCIntersection<2>>;
template class ObjectPool<SimplexDCIntersection<2>>;

template SimplexTree<2, SimplexDCLeaf<2>>*
ObjectPool<SimplexTree<2, SimplexDCLeaf<2>>, SimplexDCLeaf<2>,
    SimplexLeafSubspace<2>, SimplexDCMinEdge<2>,
    SimplexDCIntersection<2>>::get(
        SimplexTree<2, SimplexDCLeaf<2>>*, unsigned, Region<2>);
template SimplexDCLeaf<2>* ObjectPool<SimplexDCLeaf<2>,
    SimplexLeafSubspace<2>, SimplexDCMinEdge<2>,
    SimplexDCIntersection<2>>::get();
template SimplexLeafSubspace<2>* ObjectPool<SimplexLeafSubspace<2>,
    SimplexDCMinEdge<2>, SimplexDCIntersection<2>>::get();
template SimplexDCMinEdge<2>* ObjectPool<SimplexDCMinEdge<2>,
    SimplexDCIntersection<2>>::get();
template SimplexDCIntersection<2>* ObjectPool<SimplexDCIntersection<2>>::get();
}   // namespace libfive
