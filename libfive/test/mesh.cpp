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
#include <chrono>

#include "catch.hpp"

#include "libfive/render/brep/mesh.hpp"
#include "libfive/render/brep/region.hpp"

#include "util/shapes.hpp"

using namespace Kernel;

TEST_CASE("Mesh::render (sphere normals)")
{
    Tree s = sphere(0.5);
    Region<3> r({-1, -1, -1}, {1, 1, 1});

    auto mesh = Mesh::render(s, r);

    float dot = 2;
    int pos = 0;
    int neg = 0;
    for (auto t : mesh->branes)
    {
        auto norm = (mesh->verts[t(1)] - mesh->verts[t(0)])
            .cross(mesh->verts[t(2)] - mesh->verts[t(0)])
            .normalized();
        auto center = ((mesh->verts[t(0)] +
                        mesh->verts[t(1)] +
                        mesh->verts[t(2)])).normalized();
        auto dot_ = norm.dot(center);
        neg += (dot_ < 0);
        pos += (dot_ > 0);
        dot = fmin(dot, dot_);
    }
    CAPTURE(neg);
    CAPTURE(pos);
    REQUIRE(dot > 0.9);
}

TEST_CASE("Mesh::render (cylinder)")
{
  auto cube = max(max(
    max(-(Tree::X() + 1.5),
      Tree::X() - 1.5),
    max(-(Tree::Y() + 4.5),
      Tree::Y() - 4.5)),
    max(-(Tree::Z() + 1.5),
      Tree::Z() - 1.5));
  auto c = CylinderYAxis({ 0.f,0.f,0.f }, 1.f);

  auto cyl = max(-c,cube);
  Region<3> r({ -5.5, -5.5, -5.5 }, { 5.5, 5.5, 5.5 });

  auto mesh = Mesh::render(cyl, r, .1);
}

TEST_CASE("Mesh::render (cube)")
{
    auto cube = max(max(
                    max(-(Tree::X() + 1.5),
                          Tree::X() - 1.5),
                    max(-(Tree::Y() + 1.5),
                          Tree::Y() - 1.5)),
                    max(-(Tree::Z() + 1.5),
                          Tree::Z() - 1.5));
    auto s = sphere(1.95f);

    cube = max(-s, cube);
    Region<3> r({-2.5, -2.5, -2.5}, {2.5, 2.5, 2.5});

    auto mesh = Mesh::render(cube, r,.1);
}

TEST_CASE("Mesh::render (performance)")
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed;

    auto spheres = sphere(1, { 1.5, 1.5, 1.5 });
    Tree sponge = max(menger(2), -spheres);

    sponge = max(sponge, -spheres);

    Region<3> r({-2.5, -2.5, -2.5}, {2.5, 2.5, 2.5});

    // Begin timekeeping
    start = std::chrono::system_clock::now();
    auto mesh = Mesh::render(sponge, r, 0.1);
    end = std::chrono::system_clock::now();

    elapsed = end - start;

    auto elapsed_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(elapsed);

    std::string log = "\nMade sponge mesh in " +
           std::to_string(elapsed.count()) + " sec";
    WARN(log);
}

TEST_CASE("Mesh::render (gyroid)")
{

  //   float scale = 5.f, radius = (10.f / 3.f) * scale, radiusSq = radius * radius, size = radius * 1.1f, minFeat = radius / 40.f;
  //   Kernel::Region<3> region({ -size, -size, -size }, { size, size, size });
  //   Kernel::Tree box = max(max(max(max(max(Kernel::Tree::X() - radius, -Kernel::Tree::X() - radius),
  //                          Kernel::Tree::Y() - radius), -Kernel::Tree::Y() - radius),
  //                          Kernel::Tree::Z() - radius), -Kernel::Tree::Z() - radius);
  //   Kernel::Tree schwarz = cos(Kernel::Tree::X() / scale) + cos(Kernel::Tree::Y() / scale) + cos(Kernel::Tree::Z() / scale);
  //   Kernel::Tree boxschwarz = max(box, schwarz);
  //   Kernel::Tree gyroid =
  //     sin(Kernel::Tree::X() / scale) * cos(Kernel::Tree::Y() / scale) +
  //     sin(Kernel::Tree::Y() / scale) * cos(Kernel::Tree::Z() / scale) +
  //     sin(Kernel::Tree::Z() / scale) * cos(Kernel::Tree::X() / scale);
  //   Kernel::Tree boxgyroid = max(box, gyroid);
  //   Kernel::Tree sphere = -1.21 * radiusSq +
  //     Kernel::Tree::X() * Kernel::Tree::X() +
  //     Kernel::Tree::Y() * Kernel::Tree::Y() +
  //     Kernel::Tree::Z() * Kernel::Tree::Z();
  //   Kernel::Tree boxsphere = max(box, -sphere);
  // 
  //    auto mesh = Kernel::Mesh::render(boxgyroid, region, minFeat, 0.000001);



  //Brad
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed;


  auto scale = .125f;
  auto radius = 1.5f;

  Kernel::Tree gyroid =
    sin(Kernel::Tree::X() / scale) * cos(Kernel::Tree::Y() / scale) +
    sin(Kernel::Tree::Y() / scale) * cos(Kernel::Tree::Z() / scale) +
    sin(Kernel::Tree::Z() / scale) * cos(Kernel::Tree::X() / scale);

  auto sphere1 = sphere(3.0f, { 0.f,0.f,0.f });


  Kernel::Tree boxschwarz = max(sphere1, gyroid);


  Region<3> r({ -5, -5, -5 }, { 5, 5, 5 });

  // Begin timekeeping
  start = std::chrono::system_clock::now();
  auto mesh = Mesh::render(boxschwarz, r, 0.05);
  end = std::chrono::system_clock::now();

  elapsed = end - start;

  auto elapsed_ms =
    std::chrono::duration_cast<std::chrono::milliseconds>(elapsed);

  std::string log = "\nMade gyroid in " +
    std::to_string(elapsed.count()) + " sec";
  WARN(log);

  mesh->saveSTL("gyroidBlnX.stl");

}

TEST_CASE("Mesh::render (schwartz)")
{

//   float scale = 5.f, radius = (10.f / 3.f) * scale, radiusSq = radius * radius, size = radius * 1.1f, minFeat = radius / 40.f;
//   Kernel::Region<3> region({ -size, -size, -size }, { size, size, size });
//   Kernel::Tree box = max(max(max(max(max(Kernel::Tree::X() - radius, -Kernel::Tree::X() - radius),
//                          Kernel::Tree::Y() - radius), -Kernel::Tree::Y() - radius),
//                          Kernel::Tree::Z() - radius), -Kernel::Tree::Z() - radius);
//   Kernel::Tree schwarz = cos(Kernel::Tree::X() / scale) + cos(Kernel::Tree::Y() / scale) + cos(Kernel::Tree::Z() / scale);
//   Kernel::Tree boxschwarz = max(box, schwarz);
//   Kernel::Tree gyroid =
//     sin(Kernel::Tree::X() / scale) * cos(Kernel::Tree::Y() / scale) +
//     sin(Kernel::Tree::Y() / scale) * cos(Kernel::Tree::Z() / scale) +
//     sin(Kernel::Tree::Z() / scale) * cos(Kernel::Tree::X() / scale);
//   Kernel::Tree boxgyroid = max(box, gyroid);
//   Kernel::Tree sphere = -1.21 * radiusSq +
//     Kernel::Tree::X() * Kernel::Tree::X() +
//     Kernel::Tree::Y() * Kernel::Tree::Y() +
//     Kernel::Tree::Z() * Kernel::Tree::Z();
//   Kernel::Tree boxsphere = max(box, -sphere);
// 
//    auto mesh = Kernel::Mesh::render(boxgyroid, region, minFeat, 0.000001);

  

  //Brad
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed;


 auto scale = .125f;
 auto radius = 1.5f;

 
 auto sphere1 = sphere(3.0f, { 0.f,0.f,0.f });

 Kernel::Tree box = max(max(max(max(max(Kernel::Tree::X() - radius, -Kernel::Tree::X() - radius),
                        Kernel::Tree::Y() - radius), -Kernel::Tree::Y() - radius),
                        Kernel::Tree::Z() - radius), -Kernel::Tree::Z() - radius);
 Kernel::Tree schwarz = cos(Kernel::Tree::X() / scale) + cos(Kernel::Tree::Y() / scale) + cos(Kernel::Tree::Z() / scale);
 Kernel::Tree boxschwarz = max(sphere1, schwarz);


 Region<3> r({ -5, -5, -5 }, { 5, 5, 5 });
 
 // Begin timekeeping
  start = std::chrono::system_clock::now();
 auto mesh = Mesh::render(boxschwarz, r, 0.05);
  end = std::chrono::system_clock::now();
 
  elapsed = end - start;
 
 auto elapsed_ms =
 std::chrono::duration_cast<std::chrono::milliseconds>(elapsed);
 
 std::string log = "\nMade schwartz in " +
 std::to_string(elapsed.count()) + " sec";
 WARN(log);
 
 mesh->saveSTL("schwartzBlnX.stl");

}

TEST_CASE("Mesh::render (face count in rectangular prism)")
{
    auto t = max(max(max(-Tree::X(), Tree::X() - 4),
                     max(-Tree::Y(), Tree::Y() - 1)),
                     max(-Tree::Z(), Tree::Z() - 0.25));
    auto m = Mesh::render(t, Region<3>({-1, -1, -1}, {5, 2, 1.25}), 0.125);
    REQUIRE(m->verts.size() == 9); // index 0 is unused
    REQUIRE(m->branes.size() == 12);
}

TEST_CASE("Mesh::render (sphere)")
{
    auto s = sphere(1);
    auto m = Mesh::render(s, Region<3>({-1.6, -1, -8}, {1.6, 1, 1}),
                          1/32.0f, pow(10, -3));
    REQUIRE(true);
}

TEST_CASE("Mesh::export (.stl)") {
  auto cube = max(max(
    max(-(Tree::X() + 1.5),
      Tree::X() - 1.5),
    max(-(Tree::Y() + 1.5),
      Tree::Y() - 1.5)),
    max(-(Tree::Z() + 1.5),
      Tree::Z() - 1.5));
  auto s = sphere(1.95f);

  cube = max(-s, cube);
  Region<3> r({ -2.5, -2.5, -2.5 }, { 2.5, 2.5, 2.5 });

  auto mesh = Mesh::render(cube, r, .5);

  Eigen::Matrix3f m;
  m = Eigen::AngleAxisf(float(M_PI / 4), Eigen::Vector3f::UnitY()) *
    Eigen::AngleAxisf(float(atan(1 / sqrt(2))), Eigen::Vector3f::UnitX());

  //     auto mesh_ = mesh->remap(
  //       m(0, 0)*Tree::X() + m(0, 1)*Tree::Y() + m(0, 2)*Tree::Z(),
  //       m(1, 0)*Tree::X() + m(1, 1)*Tree::Y() + m(1, 2)*Tree::Z(),
  //       m(2, 0)*Tree::X() + m(2, 1)*Tree::Y() + m(2, 2)*Tree::Z());

  mesh->saveSTL("cubeSphereXp5.stl");
  mesh->saveSTL("cubeSphereAp5.stl",false);

  auto cube2 = max(max(
    max(-(Tree::X() + 1.5),
      Tree::X() - 1.5),
    max(-(Tree::Y() + 4.5),
      Tree::Y() - 4.5)),
    max(-(Tree::Z() + 1.5),
      Tree::Z() - 1.5));
  auto c = CylinderYAxis({ 0.f,0.f,0.f }, 1.f);

  auto cyl = max(-c, cube2);
  Region<3> r2({ -5.5, -5.5, -5.5 }, { 5.5, 5.5, 5.5 });

  auto cylinderMesh = Mesh::render(cyl, r2, .1);

  cylinderMesh->saveSTL("cylinderMeshX.stl");
  cylinderMesh->saveSTL("cylinderMeshA.stl",false);
}
