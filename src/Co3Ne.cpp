//
// Created by hongbo on 09/04/25.
//
#include "RVD/Co3Ne.h"
#include "RVD/Co3NeThread.h"
#include "RVD/Co3NeTools.h"
#include "RVD/Co3NeManifoldExtraction.h"
#include <geogram/basic/command_line.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_io.h>
#include "Base/MeshAdaptor.h"
#include <Eigen/Dense>
#include <geogram/points/co3ne.h>
namespace GEO_BASE{
    GEO_BASE::Co3Ne::Co3Ne(GEO::Mesh& M):mesh_(M){
        RVD_.init(M);
        GEO::index_t nb = GEO::Process::maximum_concurrent_threads();
        thread_.clear();
        GEO::index_t batch_size = RVD_.nb_points() / nb;
        GEO::index_t cur = 0;
        GEO::index_t remaining = RVD_.nb_points();
        for(GEO::index_t i = 0; i < nb; i++) {
            GEO::index_t this_batch_size = batch_size;
            if(i == nb - 1) {
                this_batch_size = remaining;
            }
            thread_.push_back(
                    new Co3NeThread(
                            this, cur, cur + this_batch_size
                    )
            );
            cur += this_batch_size;
            remaining -= this_batch_size;
        }
        geo_assert(remaining == 0);

        // TODO: pass it as an argument and let Vorpaline's main.cpp
        // communicate with CmdLine.
        double alpha = GEO::CmdLine::get_arg_double("co3ne:max_N_angle");
        alpha = alpha * M_PI / 180.0;
        set_max_angle(alpha);
    }

    void Co3Ne::reconstructhd(double r) {
        GEO::ProgressTask progress("reconstruct",100);
        GEO::Stopwatch W("Co3Ne recons");
        GEO::Logger::out("Co3Ne")
                << "construct tangent plane"
                << std::endl;

        RVD_.set_nb_neighbors(40);
        RVD_.set_circles_radius(r);

        for(GEO::index_t t = 0; t < thread_.size(); t++) {
            thread_[t]->set_mode(CO3NE_TANGENT_AND_RECONSTRUCT);
            thread_[t]->triangles().clear();
        }
        progress.progress(1);
        run_threads();
        progress.progress(50);
    }

    void Co3Ne::reconstruct3d(double r) {
        GEO::mesh_repair(mesh_,GEO::MESH_REPAIR_COLOCATE, 1e-8*r);
        GEO::ProgressTask progress("reconstruct",100);
        GEO::Stopwatch W("Co3Ne recons");
        GEO::Logger::out("Co3Ne")
                << "construct tangent plane"
                << std::endl;

        RVD_.set_nb_neighbors(30);
        RVD_.set_circles_radius(r);

        for(GEO::index_t t = 0; t < thread_.size(); t++) {
            thread_[t]->set_mode(CO3NE_NORMALS_AND_RECONSTRUCT);
            thread_[t]->triangles().clear();
        }
        progress.progress(1);
        run_threads();
        progress.progress(50);
    }


    void Co3Ne::run_threads() {
        GEO::Process::run_threads(thread_);
    }

    void Co3Ne::extract_triangles() {
        GEO::ProgressTask progress("Manifold Extraction",100);
        GEO::Stopwatch W("Co3Ne manif.");
        RVD_.clear();  // reclaim memory used by ANN

        GEO::index_t nb_triangles = 0;
        for(GEO::index_t t = 0; t < thread_.size(); t++) {
            nb_triangles += thread_[t]->nb_triangles();
        }

        GEO::Logger::out("Co3Ne") << "Raw triangles: "
                             << nb_triangles
                             << std::endl;

        raw_triangles_.reserve(nb_triangles * 3);
        for(GEO::index_t th = 0; th < thread_.size(); th++) {
            GEO::vector<GEO::index_t>& triangles = thread_[th]->triangles();
            raw_triangles_.insert(
                    raw_triangles_.end(),
                    triangles.begin(), triangles.end()
            );
            thread_[th]->triangles().clear();
        }

        progress.progress(53);


        std::vector<GEO::index_t> not_so_good_triangles;
        Co3NeTools::co3ne_split_triangles_list(
                raw_triangles_, T3_triangles_, T12_triangles_
        );

//        GEO::Mesh M_t3;
//        M_t3.vertices.assign_points(
//                mesh_.vertices.point_ptr(0),
//                mesh_.vertices.dimension(),
//                mesh_.vertices.nb()
//        );
//        M_t3.facets.assign_triangle_mesh(T3_triangles_, false);
//        M_t3.vertices.set_dimension(3);
//        GEO::mesh_save(M_t3, "/Users/lihongbo/Desktop/code/SimplexCVT_recon/tmp/co3ne_T3.obj");
//
//        GEO::Mesh M_t12;
//        M_t12.vertices.assign_points(
//                mesh_.vertices.point_ptr(0),
//                mesh_.vertices.dimension(),
//                mesh_.vertices.nb()
//        );
//        M_t12.facets.assign_triangle_mesh(T12_triangles_, false);
//        M_t12.vertices.set_dimension(3);
//        GEO::mesh_save(M_t12, "/Users/lihongbo/Desktop/code/SimplexCVT_recon/tmp/co3ne_T12.obj");


        progress.progress(83);

        Co3NeManifoldExtraction manifold_extraction(
                mesh_, T3_triangles_, false
        );

        progress.progress(55);


        manifold_extraction.add_triangles(T12_triangles_);


        progress.progress(57);

        mesh_reorient(mesh_);

        progress.progress(100);

    }

    void Co3Ne::post_process() {
        GEO::mesh_repair(mesh_,GEO::MeshRepairMode(
                GEO::MESH_REPAIR_DEFAULT | GEO::MESH_REPAIR_RECONSTRUCT
        ));

//        GEO::mesh_repair(mesh_,GEO::MeshRepairMode(GEO::MESH_REPAIR_DEFAULT));

//        GEO::mesh_postprocess_RDT(mesh_, true);
    }

    void Co3Ne::save_raw_triangles(const std::string &filename) {
        GEO::Logger::out("Co3Ne") << ">> co3ne_raw"
                                  << std::endl;
        Base::MeshAdaptor::SaveMesh(mesh_, raw_triangles_, filename);
    }

    void Co3Ne::save_T3_triangles(const std::string &filename) {
        GEO::Logger::out("Co3Ne") << ">> output good triangles"
                             << std::endl;
        Base::MeshAdaptor::SaveMesh(mesh_, T3_triangles_, filename);
    }

    void Co3Ne::save_T12_triangles(const std::string &filename) {
        GEO::Logger::out("Co3Ne") << ">> output T12 triangles"
                             << std::endl;
        Base::MeshAdaptor::SaveMesh(mesh_, T12_triangles_, filename);
    }

    void Co3Ne::save_final_mesh(const std::string &filename) {
        GEO::Logger::out("Co3Ne") << ">> output Final triangles"
                                  << std::endl;
        Base::MeshAdaptor::SaveMesh(mesh_,filename);

    }
}