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
#pragma once

#include <Eigen/Eigen>

#include "libfive/tree/tree.hpp"

//Ops:

Kernel::Tree OpMix(Kernel::Tree tA, Kernel::Tree tB, float a);

//Shapes:
Kernel::Tree rectangle(float xmin, float xmax, float ymin, float ymax,
                       Eigen::Matrix4f M=Eigen::Matrix4f::Identity());
Kernel::Tree menger(int i);
Kernel::Tree circle(float r);
Kernel::Tree sphere(float r, Eigen::Vector3f center=Eigen::Vector3f::Zero());
Kernel::Tree box(const Eigen::Vector3f& lower, const Eigen::Vector3f& upper);

Kernel::Tree CylinderYAxis(Eigen::Vector3f start, float r);

Kernel::Tree rotate2d(Kernel::Tree t, float angle);
Kernel::Tree move(Kernel::Tree t, Eigen::Vector3f d);


//CSG:
Kernel::Tree CSGUnion(Kernel::Tree tA, Kernel::Tree tB);
Kernel::Tree CSGSubtract(Kernel::Tree tA, Kernel::Tree tB);
Kernel::Tree CSGIntersect(Kernel::Tree tA, Kernel::Tree tB);

Kernel::Tree CSGUnionRound(Kernel::Tree tA, Kernel::Tree tB, float r);
Kernel::Tree CSGUnionChamfer(Kernel::Tree tA, Kernel::Tree tB, float r);

//Offset: expand or contract t by r
Kernel::Tree offset(Kernel::Tree t, float r);

//Shell: returns the shell of a given tree w/ a given offset
//shell t by r
Kernel::Tree shell(Kernel::Tree t, float r);

//clearence: make sure A clears B by amount R
//Expands Tb by r & then subtracts it from A 
Kernel::Tree clearence(Kernel::Tree tA, Kernel::Tree tB, float r);

//blend: blends tA with tB by r
Kernel::Tree blend(Kernel::Tree tA, Kernel::Tree tB, float r);

//loft: creates a surface between tA & tB in the z-Axis from zMin to zMax
Kernel::Tree loft(Kernel::Tree tA, Kernel::Tree tB, float zMin, float zMax);
Kernel::Tree loftBetween(Kernel::Tree tA, 
                         Kernel::Tree tB, 
                         const Eigen::Vector3f& lower, 
                         const Eigen::Vector3f& upper);


//Transforms
//Scale, Rotate, Move, Shear etc...
