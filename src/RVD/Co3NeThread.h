//
// Created by hongbo on 09/04/25.
//

#ifndef SIMPLEXCVT_RECON_CO3NETHREAD_H
#define SIMPLEXCVT_RECON_CO3NETHREAD_H
#include <geogram/basic/process.h>
#include "RVD/Co3NeTools.h"
#include "RVD/Co3Ne.h"
#include <geogram/points/principal_axes.h>
namespace GEO_BASE{
    class Co3NeThread: public GEO::Thread {
    public:
        Co3NeThread(
                Co3Ne* master,
                GEO::index_t from, GEO::index_t to
        );

        void run() override;

        void run_tangent_and_reconstruct();
        void run_normals_and_reconstruct();

    private:
        Co3Ne* master_;
        GEO::index_t from_;
        GEO::index_t to_;
        Co3Ne::Co3NeMode mode_;
        GEO::PrincipalAxes3d least_squares_normal_;
        GEO::vector<GEO::index_t> triangles_;

    public:
        void set_mode(Co3Ne::Co3NeMode m) {
            mode_ = m;
        }

        GEO::vector<GEO::index_t>& triangles() {
            return triangles_;
        }

        GEO::index_t nb_triangles() const {
            return triangles_.size()/3;
        }

    };
}

#endif //SIMPLEXCVT_RECON_CO3NETHREAD_H
