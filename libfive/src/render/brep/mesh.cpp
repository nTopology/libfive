/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <numeric>
#include <fstream>
#include <boost/algorithm/string/predicate.hpp>

#include "libfive/eval/evaluator.hpp"

#include "libfive/render/brep/mesh.hpp"
#include "libfive/render/brep/dual.hpp"
#include "libfive/render/brep/region.hpp"
#include "libfive/render/brep/settings.hpp"

// Dual contouring
#include "libfive/render/brep/dc/dc_worker_pool.hpp"
#include "libfive/render/brep/dc/dc_mesher.hpp"

// Simplex meshing
#include "libfive/render/brep/simplex/simplex_worker_pool.hpp"
#include "libfive/render/brep/simplex/simplex_mesher.hpp"

// Hybrid meshing
#include "libfive/render/brep/hybrid/hybrid_worker_pool.hpp"
#include "libfive/render/brep/hybrid/hybrid_mesher.hpp"

// Simplex DC meshing
#include "libfive/render/brep/simplex_dc/simplex_dc_worker_pool.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_edges.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_sub_consolidator.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_collapser.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_intersecter.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_vertexer.hpp"
#include "libfive/render/brep/simplex_dc/simplex_dc_mesher.hpp"

namespace libfive {

std::unique_ptr<Mesh> Mesh::render(const Tree t, const Region<3>& r,
                                   const BRepSettings& settings)
{
    std::vector<Evaluator, Eigen::aligned_allocator<Evaluator>> es;
    es.reserve(settings.workers);
    for (unsigned i=0; i < settings.workers; ++i) {
        es.emplace_back(Evaluator(t));
    }

    return render(es.data(), r, settings);
}

std::unique_ptr<Mesh> Mesh::render(
        Evaluator* es,
        const Region<3>& r, const BRepSettings& settings)
{
    std::unique_ptr<Mesh> out;
    if (settings.alg == DUAL_CONTOURING)
    {
        if (settings.progress_handler) {
            // Pool::build, Dual::walk, t.reset
            settings.progress_handler->start({1, 1, 1});
        }
        auto t = DCWorkerPool<3>::build(es, r, settings);

        if (settings.cancel.load() || t.get() == nullptr) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            return nullptr;
        }

        // Perform marching squares
        out = Dual<3>::walk<DCMesher>(t, settings);

        // TODO: check for early return here again
        t.reset(settings);
    }
    else if (settings.alg == ISO_SIMPLEX)
    {
        if (settings.progress_handler) {
            // Pool::build, Dual::walk, t->assignIndices, t.reset
            settings.progress_handler->start({1, 1, 1});
        }
        auto t = SimplexWorkerPool<3>::build(es, r, settings);

        if (settings.cancel.load() || t.get() == nullptr) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            return nullptr;
        }

        t->assignIndices(settings);

        out = Dual<3>::walk_<SimplexMesher>(t, settings,
                [&](PerThreadBRep<3>& brep, int i) {
                    return SimplexMesher(brep, &es[i]);
                });
        t.reset(settings);
    }
    else if (settings.alg == HYBRID)
    {
        if (settings.progress_handler) {
            // Pool::build, Dual::walk, t->assignIndices, t.reset
            settings.progress_handler->start({1, 1, 1});
        }
        auto t = HybridWorkerPool<3>::build(es, r, settings);

        if (settings.cancel.load() || t.get() == nullptr) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            return nullptr;
        }

        t->assignIndices(settings);

        out = Dual<3>::walk_<HybridMesher>(t, settings,
                [&](PerThreadBRep<3>& brep, int i) {
                    return HybridMesher(brep, &es[i]);
                });
        t.reset(settings);
    }
    else if (settings.alg == SIMPLEX_DC) {
        if (settings.progress_handler) {
            // Pool::build, Edges, SubConsolidator Collapser 2 and 3, 
            // Intersecter, Vertexer, Mesher, t.reset
            settings.progress_handler->start({ 1, 1, 1, 1, 1, 1, 1, 1, 1 });
        }
        auto t = SimplexDCWorkerPool<3>::build(es, r, settings);
        if (settings.cancel.load() || t.get() == nullptr) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            return nullptr;
        }

        // Edges
        Dual<3>::walk<SimplexDCEdges<3>, decltype(t.pool())> (
            t, settings, t.pool());

        if (settings.cancel.load()) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            t.reset(settings);
            return nullptr;
        }

        Dual<3>::walk<SimplexDCSubConsolidator<3>>(t, settings);
        
        // Collapser for edges and faces
        Dual<3>::walk_<SimplexDCCollapser<3, 2>>(
            t, settings, 
            [&](SimplexDCCollapser<3, 2>::PerThreadOutput& out, int i) {
            return SimplexDCCollapser<3, 2>(out, &es[i]);
        });
        // Collapser for cells
        Dual<3>::walk_<SimplexDCCollapser<3, 3>>(
            t, settings,
            [&](SimplexDCCollapser<3, 3>::PerThreadOutput& out, int i) {
            return SimplexDCCollapser<3, 3>(out, &es[i]);
        });

        if (settings.cancel.load()) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            t.reset(settings);
            return nullptr;
        }

        // To do when it's time for it: Degeneracy handling (collapse or center)

        // Get the intersections
        auto intersectionsBRep = Dual<3>::walk_<SimplexDCIntersecter<3>>(
            t, settings, 
            [&](SimplexDCIntersecter<3>::PerThreadOutput& brep, int i) {
            return SimplexDCIntersecter<3>(
                brep, &es[i], t.pool(), {});
        });

        if (settings.cancel.load()) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            t.reset(settings);
            return nullptr;
        }

        // Get the vertices
        auto vertsBRep = Dual<3>::walk<SimplexDCVertexer<3>>(
            t, settings, intersectionsBRep->verts.size() - 1/*Vert offset*/);

        if (settings.cancel.load()) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            t.reset(settings);
            return nullptr;
        }

        // Connect the vertices
        out = Dual<3>::walk<SimplexDCMesher>(t, settings);

        if (settings.cancel.load()) {
            if (settings.progress_handler) {
                settings.progress_handler->finish();
            }
            t.reset(settings);
            return nullptr;
        }

        t.reset(settings);

        out->verts = std::move(intersectionsBRep->verts);
        assert(!vertsBRep->verts.empty());
        out->verts.insert(out->verts.end(), 
            vertsBRep->verts.begin() + 1, vertsBRep->verts.end());
    }

    if (settings.progress_handler) {
        settings.progress_handler->finish();
    }
    return out;
}

void Mesh::line(const Eigen::Vector3f& a, const Eigen::Vector3f& b)
{
    uint32_t a_ = verts.size();
    verts.push_back(a);
    uint32_t b_ = verts.size();
    verts.push_back(b);

    branes.push_back({a_, a_, b_});
}

////////////////////////////////////////////////////////////////////////////////

bool Mesh::saveSTL(const std::string& filename,
                   const std::list<const Mesh*>& meshes)
{
    if (!boost::algorithm::iends_with(filename, ".stl"))
    {
        std::cerr << "Mesh::saveSTL: filename \"" << filename
                  << "\" does not end in .stl" << std::endl;
    }
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::binary);
    if (!file.is_open())
    {
        std::cout << "Mesh::saveSTL: could not open " << filename
                  << std::endl;
        return false;
    }

    // File header (giving human-readable info about file type)
    std::string header = "This is a binary STL exported from libfive.";
    file.write(header.c_str(), header.length());

    // Pad the rest of the header to 80 bytes
    for (int i=header.length(); i < 80; ++i)
    {
        file.put(' ');
    }

    // Write the triangle count to the file
    uint32_t num = std::accumulate(meshes.begin(), meshes.end(), (uint32_t)0,
            [](uint32_t i, const Mesh* m){ return i + m->branes.size(); });
    file.write(reinterpret_cast<char*>(&num), sizeof(num));

    for (const auto& m : meshes)
    {
        for (const auto& t : m->branes)
        {
            // Write out the normal vector for this face (all zeros)
            float norm[3] = {0, 0, 0};
            file.write(reinterpret_cast<char*>(&norm), sizeof(norm));

            // Iterate over vertices (which are indices into the verts list)
            for (unsigned i=0; i < 3; ++i)
            {
                auto v = m->verts[t[i]];
                float vert[3] = {v.x(), v.y(), v.z()};
                file.write(reinterpret_cast<char*>(&vert), sizeof(vert));
            }

            // Write out this face's attribute short
            uint16_t attrib = 0;
            file.write(reinterpret_cast<char*>(&attrib), sizeof(attrib));
        }
    }

    return true;
}

bool Mesh::saveSTL(const std::string& filename) const
{
    return saveSTL(filename, {this});
}

}   // namespace libfive
