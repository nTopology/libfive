/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2019  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <atomic>
#include <memory>

namespace libfive {

// Forward declarations
class ProgressHandler;
class FreeThreadHandler;
class VolTree;

enum BRepAlgorithm {
    DUAL_CONTOURING,
    ISO_SIMPLEX,
    HYBRID,
    SIMPLEX_DC,
};

struct BRepSettings {
public:
    BRepSettings()
    {
        reset();
    }

    void reset() {
        min_feature = 0.1;
        max_err = 1e-8;
        simplex_dc_padding_rate = 1e-4;
        simplex_bounding_eigenvalue_cutoff = 0.;
        workers = 8;
        alg = DUAL_CONTOURING;
        free_thread_handler = nullptr;
        progress_handler = nullptr;
        cancel.store(false);
        vol = nullptr;
    }

    /*  The meshing region is subdivided until the smallest region edge
     *  is below min_feature in size.  Make this smaller to get a
     *  higher-resolution model. */
    double min_feature;

    /*  This value is used when deciding whether to collapse cells.  If it
     *  is very small, then only linear regions are merged.  Set as -1 to
     *  completely disable cell merging.  */
    double max_err;

    /*  This value is used in simplex DC meshing when deciding whether to 
     *  collapse nearby intersection vertices each other, as well as how close
     *  to the faces of a simplex the output vertex corresponding to that 
     *  simplex can be.  If it is very small, the resulting mesh will be more
     *  accurate, but will have more tiny triangles, as well as more triangles 
     *  overall.  It should never be more than 0.5.*/
    double simplex_dc_padding_rate;

    /*  This value is used in simplex and simplex DC meshing to set a custom
     *  (typically larger) relative eigenvalue cutoff when minimizing QEFs, but
     *  only to the extent required to ensure that the resulting point is in the
     *  appropriate region.  When assigned its default value of 0, this
     *  feature is not used, and the default cutoff is used regardless of whether
     *  the point will be in the appropriate region without it.*/
    double simplex_bounding_eigenvalue_cutoff;

    /*  Number of worker threads to use while meshing.  Set as 0 to use the
     *  platform-default number of threads. */
    unsigned workers;

    /*  This is the meshing algorithm */
    BRepAlgorithm alg;

    /*  Optional function called when a thread finds itself without anything
     *  to do.  This can be used to keep threads from spinning if libfive
     *  is embedded in a larger application with its own pooling system. */
    FreeThreadHandler* free_thread_handler;

    /*  Optional class that wraps a progress callback.  */
    ProgressHandler* progress_handler;

    /*  Optional acceleration structure */
    const VolTree* vol;

    mutable std::atomic_bool cancel;
};

}
