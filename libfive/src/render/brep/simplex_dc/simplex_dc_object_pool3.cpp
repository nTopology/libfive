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
template class ObjectPool<SimplexTree<3, SimplexDCLeaf<3>>, SimplexDCLeaf<3>, 
                          SimplexLeafSubspace<3>, SimplexDCMinEdge<3>, 
                          SimplexDCIntersection<3>>;
template class ObjectPool<SimplexDCLeaf<3>, SimplexLeafSubspace<3>, 
                          SimplexDCMinEdge<3>, SimplexDCIntersection<3>>;
template class ObjectPool<SimplexLeafSubspace<3>, SimplexDCMinEdge<3>,
                          SimplexDCIntersection<3>>;
template class ObjectPool<SimplexDCMinEdge<3>, SimplexDCIntersection<3>>;
template class ObjectPool<SimplexDCIntersection<3>>;

template SimplexTree<3, SimplexDCLeaf<3>>*
ObjectPool<SimplexTree<3, SimplexDCLeaf<3>>, SimplexDCLeaf<3>,
           SimplexLeafSubspace<3>, SimplexDCMinEdge<3>, 
           SimplexDCIntersection<3>>::get(
        SimplexTree<3, SimplexDCLeaf<3>>*, unsigned, Region<3>);
template SimplexDCLeaf<3>* ObjectPool<SimplexDCLeaf<3>,
    SimplexLeafSubspace<3>, SimplexDCMinEdge<3>,
    SimplexDCIntersection<3>>::get();
template SimplexLeafSubspace<3>* ObjectPool<SimplexLeafSubspace<3>,
    SimplexDCMinEdge<3>, SimplexDCIntersection<3>>::get();
template SimplexDCMinEdge<3>* ObjectPool<SimplexDCMinEdge<3>, 
    SimplexDCIntersection<3>>::get();
template SimplexDCIntersection<3>* ObjectPool<SimplexDCIntersection<3>>::get();
}   // namespace libfive
