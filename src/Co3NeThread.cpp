//
// Created by hongbo on 10/04/25.
//
#include "RVD/Co3NeThread.h"
#include "RVD/Geometry.h"
#include <geogram/points/co3ne.h>
namespace GEO_BASE{
    Co3NeThread::Co3NeThread(GEO_BASE::Co3Ne *master,
                             GEO::index_t from,
                             GEO::index_t to):
            master_(master),
            from_(from),
            to_(to)  {

    }

    void Co3NeThread::run() {
        switch(mode_) {
//            case Co3Ne::CO3NE_NORMALS:
//            run_normals();
//            break;
//            case Co3Ne::CO3NE_SMOOTH:
//            run_smooth();
//            break;
//            case Co3Ne::CO3NE_RECONSTRUCT:
//                run_reconstruct();
//            break;
            case Co3Ne::CO3NE_NORMALS_AND_RECONSTRUCT:
                run_normals_and_reconstruct();
            break;
            case Co3Ne::CO3NE_TANGENT_AND_RECONSTRUCT:
                run_tangent_and_reconstruct();
            break;
            case Co3Ne::CO3NE_NONE:
            break;
        }
    }

    void Co3NeThread::run_tangent_and_reconstruct() {
        Co3NeRestrictedVoronoiDiagram& RVD = master_->RVD();
        GEO::index_t nb_neigh = RVD.nb_neighbors();
        GEO::vector<GEO::index_t> neigh(100);
        GEO::vector<double> sq_dist(100);
        Polygon_hd P(100);
        Polygon_hd Q(100);

        for(GEO::index_t i = from_; i < to_; i++) {
            RVD.get_neighbors(
                    i, neigh, sq_dist, nb_neigh
            );

            Eigen::MatrixXd nb_points(nb_neigh, master_->dim());
            for(GEO::index_t jj = 0; jj < neigh.size(); jj++) {
                for(GEO::index_t kk = 0; kk < 6; kk++) {
                    nb_points(jj, kk) = RVD.point_ptr(neigh[jj])[kk];
                }
            }

            Eigen::MatrixXd tangent_basis = Co3NeTools::PCA_fit(nb_points);

            RVD.get_RVC(i,tangent_basis,P,Q,neigh,sq_dist);

//
//            RVD.get_RVC(i, N, P, Q, neigh, sq_dist);
//            if(debug_RVD) {
//                for(index_t v = 0; v < P.nb_vertices(); ++v) {
//                    RVD_file << "v "
//                             << P.vertex(v).point().x
//                             << " "
//                             << P.vertex(v).point().y
//                             << " "
//                             << P.vertex(v).point().z
//                             << std::endl;
//                }
//                RVD_file << "f ";
//                for(index_t v = 0; v < P.nb_vertices(); ++v) {
//                    ++cur_v;
//                    RVD_file << cur_v << " ";
//                }
//                RVD_file << std::endl;
//                RVD_file << "#" << i << " ";
//                for(index_t v1 = 0; v1 < P.nb_vertices(); ++v1) {
//                    RVD_file << P.vertex(v1).adjacent_seed() << " ";
//                }
//                RVD_file << std::endl;
//            }
            for(GEO::index_t v1 = 0; v1 < P.nb_vertices(); v1++) {
                GEO::index_t v2 = P.next_vertex(v1);
                GEO::signed_index_t j = P.vertex(v1).adjacent_seed();
                GEO::signed_index_t k = P.vertex(v2).adjacent_seed();
                if(j >= 0 && k >= 0 && j != k) {
                    triangles_.push_back(i);
                    triangles_.push_back(GEO::index_t(j));
                    triangles_.push_back(GEO::index_t(k));
                }
            }
        }
        auto aaa = 0;
        &aaa;
    }

    void Co3NeThread::run_normals_and_reconstruct(){
        GEO::Attribute<double> normal;

        GEO::index_t cur_v = 0;

        Co3NeRestrictedVoronoiDiagram& RVD = master_->RVD();
        GEO::index_t nb_neigh = RVD.nb_neighbors();
        GEO::vector<GEO::index_t> neigh(100);
        GEO::vector<double> sq_dist(100);
        GEO_BASE::Polygon P(100);
        GEO_BASE::Polygon Q(100);

        for(GEO::index_t i = from_; i < to_; i++) {

            GEO::vec3 N;
            if(normal.is_bound()) {
                RVD.get_neighbors(
                        i, neigh, sq_dist, nb_neigh
                );
                N = GEO::vec3(normal[3*i], normal[3*i+1], normal[3*i+2]);
            } else {
                RVD.get_neighbors(
                        i, neigh, sq_dist, nb_neigh
                );
                least_squares_normal_.begin();
                for(GEO::index_t jj = 0; jj < neigh.size(); jj++) {
                    least_squares_normal_.add_point(RVD.point(neigh[jj]));
                }
                least_squares_normal_.end();
                N = least_squares_normal_.normal();
            }

            RVD.get_RVC(i, N, P, Q, neigh, sq_dist);

//            if(debug_RVD) {
//                for(index_t v = 0; v < P.nb_vertices(); ++v) {
//                    RVD_file << "v "
//                             << P.vertex(v).point().x
//                             << " "
//                             << P.vertex(v).point().y
//                             << " "
//                             << P.vertex(v).point().z
//                             << std::endl;
//                }
//                RVD_file << "f ";
//                for(index_t v = 0; v < P.nb_vertices(); ++v) {
//                    ++cur_v;
//                    RVD_file << cur_v << " ";
//                }
//                RVD_file << std::endl;
//                RVD_file << "#" << i << " ";
//                for(index_t v1 = 0; v1 < P.nb_vertices(); ++v1) {
//                    RVD_file << P.vertex(v1).adjacent_seed() << " ";
//                }
//                RVD_file << std::endl;
//            }


            for(GEO::index_t v1 = 0; v1 < P.nb_vertices(); v1++) {
                GEO::index_t v2 = P.next_vertex(v1);
                GEO::signed_index_t j = P.vertex(v1).adjacent_seed();
                GEO::signed_index_t k = P.vertex(v2).adjacent_seed();
                if(
                        j >= 0 && k >= 0 && j != k
                        ) {
                    triangles_.push_back(i);
                    triangles_.push_back(GEO::index_t(j));
                    triangles_.push_back(GEO::index_t(k));
                }
            }
        }

        if(normal.is_bound()) {
            normal.unbind();
        }
    }
}