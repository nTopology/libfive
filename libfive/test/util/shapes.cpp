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
#include "libfive/tree/tree.hpp"

#include "shapes.hpp"

using namespace Kernel;

Tree rectangle(float xmin, float xmax, float ymin, float ymax,
               Eigen::Matrix4f M)
{
    auto x = M(0,0)*Tree::X() + M(0,1)*Tree::Y() + M(0,2)*Tree::Z() + M(0,3);
    auto y = M(1,0)*Tree::X() + M(1,1)*Tree::Y() + M(1,2)*Tree::Z() + M(1,3);

    return max(max(xmin - x, x - xmax), max(ymin - y, y - ymax));
}

Tree rotate2d(Tree t, float angle)
{
    return t.remap( cos(angle) * Tree::X() + sin(angle) * Tree::Y(),
                   -sin(angle) * Tree::X() + cos(angle) * Tree::Y(),
                   Tree::Z());
}

Tree move(Tree t, Eigen::Vector3f m)
{
    return t.remap(Tree::X() - m.x(), Tree::Y() - m.y(), Tree::Z() - m.z());
}

Kernel::Tree CSGUnion(Kernel::Tree tA, Kernel::Tree tB)
{
 return min(tA, tB);
}

Kernel::Tree CSGSubtract(Kernel::Tree tA, Kernel::Tree tB)
{
  return CSGIntersect(tA, -tB);
}

Kernel::Tree CSGIntersect(Kernel::Tree tA, Kernel::Tree tB)
{
  return max(tA, tB);
}

Kernel::Tree CSGUnionRound(Kernel::Tree tA, Kernel::Tree tB, float r)
{
  auto vc0 = r - tA;
  auto vc1 = r - tB;

  auto u0 = max(vc0, 0.f);
  auto u1 = max(vc1, 0.f);
   
  auto len = sqrt(square(u0) + square(u1));

  return max(r, min(tA, tB)) - len;
//   // The "Round" variant uses a quarter-circle to join the two objects smoothly:
//   float fOpUnionRound(float a, float b, float r)
//   {
//     vec2 u = max(vec2(r - a, r - b), vec2(0));
//     return max(r, min(a, b)) - length(u);
//   }
}

Kernel::Tree CSGUnionChamfer(Kernel::Tree tA, Kernel::Tree tB, float r)
{
  return min(min(tA, tB), (tA - r + tB)*sqrt(0.5f));
}

Kernel::Tree offset(Kernel::Tree t, float r)
{
  return t + r;
}

Kernel::Tree shell(Kernel::Tree t, float r)
{
  return clearence(t, t, r);
}

Kernel::Tree clearence(Kernel::Tree tA, Kernel::Tree tB, float r)
{
  return CSGSubtract(tA, offset(tB, r));
}

Kernel::Tree blend(Kernel::Tree tA, Kernel::Tree tB, float r)
{
  return CSGUnion(tA, 
                  CSGUnion(tB, (sqrt(abs(tA)) + sqrt(abs(tB))) - r));
}

Kernel::Tree loft(Kernel::Tree tA, Kernel::Tree tB, float zMin, float zMax)
{
  return max(max((Tree::Z() - zMax), 
                 (zMin - Tree::Z())),
             (((Tree::Z() - zMin)*tB) + 
             ((zMax - Tree::Z())*tA)) 
             / (zMax - zMin));
}

Tree recurse(float x, float y, float scale, Eigen::Matrix4f M, int i)
{
    auto base = rectangle(x - scale/2, x + scale/2,
                          y - scale/2, y + scale/2, M);

    if (i == 0)
    {
        return base;
    }
    else
    {
        auto j = i - 1;
        auto t = scale / 3;

        return min(base,
               min(recurse(x + scale, y, t, M, j),
               min(recurse(x - scale, y, t, M, j),
               min(recurse(x, y + scale, t, M, j),
               min(recurse(x, y - scale, t, M, j),
               min(recurse(x + scale, y + scale, t, M, j),
               min(recurse(x + scale, y - scale, t, M, j),
               min(recurse(x - scale, y + scale, t, M, j),
                   recurse(x - scale, y - scale, t, M, j)
               ))))))));
    }
}

Tree menger(int i)
{
    Eigen::Matrix3f m = Eigen::Matrix3f::Identity();
    Eigen::Matrix4f M = Eigen::Matrix4f::Zero();
    M.block<3,3>(0,0) = m;
    Tree a = recurse(0, 0, 1, M, i);

    m = Eigen::AngleAxisf(float(M_PI/2), Eigen::Vector3f::UnitX());
    M.block<3,3>(0,0) = m;
    Tree b = recurse(0, 0, 1, M, i);

    m = Eigen::AngleAxisf(float(M_PI/2), Eigen::Vector3f::UnitY());
    M.block<3,3>(0,0) = m;
    Tree c = recurse(0, 0, 1, M, i);

    auto cube = max(max(
                    max(-(Tree::X() + 1.5),
                          Tree::X() - 1.5),
                    max(-(Tree::Y() + 1.5),
                          Tree::Y() - 1.5)),
                    max(-(Tree::Z() + 1.5),
                          Tree::Z() - 1.5));

    auto cutout = -min(min(a, b), c);
    return max(cube, cutout);
}

Tree circle(float r)
{
    return sqrt(square(Tree::X()) + square(Tree::Y())) - r;
}

Tree sphere(float r, Eigen::Vector3f center)
{
    return sqrt(square(Tree::X() - center.x()) +
                square(Tree::Y() - center.y()) +
                square(Tree::Z() - center.z())) - r;
}

Tree box(const Eigen::Vector3f& lower, const Eigen::Vector3f& upper)
{
    return max(max(
               max(lower.x() - Tree::X(),
                   Tree::X() - upper.x()),
               max(lower.y() - Tree::Y(),
                   Tree::Y() - upper.y())),
               max(lower.z() - Tree::Z(),
                   Tree::Z() - upper.z()));
}

Tree CylinderYAxis(Eigen::Vector3f start, float r) {   
  return r*r 
    - square(Tree::X() - start.x()) 
    - square(Tree::Z() - start.z());
}
