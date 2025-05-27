//
// Created by hongbo on 10/04/25.
//
#include "RVD/Co3NeRestrictedVoronoiDiagram.h"
#include "RVD/Co3NeTools.h"

namespace GEO_BASE{
    Co3NeRestrictedVoronoiDiagram::Co3NeRestrictedVoronoiDiagram() :
            nb_points_(0),
            p_(nullptr),
            p_stride_(0),
            n_(nullptr),
            n_stride_(0),
            radius_(0.0),
            sqROS_(0.0),
            nb_neighbors_(0) {
    }

    void Co3NeRestrictedVoronoiDiagram::init(GEO::Mesh& M) {
        geo_assert(M.vertices.dimension() >= 3);
//        geo_assert(M.vertices.dimension() == 3 || NN_->stride_supported());

        NN_ = GEO::NearestNeighborSearch::create(M.vertices.dimension());
        NN3d_ = GEO::NearestNeighborSearch::create(3);

        double* normals_pointer = nullptr;
        {
            GEO::Attribute<double> normal;
            normal.bind_if_is_defined(M.vertices.attributes(), "normal");
            if(normal.is_bound() && normal.dimension() == 3) {
                normals_pointer = &normal[0];
            }
        }

        if(normals_pointer == nullptr) {
            init(
                    M.vertices.nb(),
                    M.vertices.point_ptr(0), M.vertices.dimension(),
                    nullptr, 0
            );
        } else {
            init(
                    M.vertices.nb(),
                    M.vertices.point_ptr(0), M.vertices.dimension(),
                    normals_pointer, 3
            );
        }
    }

    void Co3NeRestrictedVoronoiDiagram::init(
            GEO::index_t nb_points_in,
            double* p, GEO::index_t p_stride,
            double* n, GEO::index_t n_stride) {
        nb_points_ = nb_points_in;
        p_ = p;
        p_stride_ = p_stride;
        n_ = n;
        n_stride_ = n_stride;

        radius_ = 0.0;
        sqROS_ = 0.0;
        nb_neighbors_ = 0;

        NN_->set_points(nb_points_, p_, p_stride_);

        //create 3d nearest neighbor search
        p3d_ = Eigen::MatrixXd(nb_points_in, 3);
        for(auto i=0;i<nb_points_;i++){
            double* p_hd = p_ + i * p_stride_;
            p3d_.row(i) = Eigen::Vector3d(p_hd[0], p_hd[1], p_hd[2]);
        }

        NN3d_->set_points(nb_points_, p3d_.data(), 3);

    }

    void Co3NeRestrictedVoronoiDiagram::get_circle(GEO::index_t i,
                                                   GEO_BASE::Polygon &P,
                                                   const GEO::vec3 &N) const {
        P.clear();
        const GEO::vec3& pi = point(i);
        GEO::vec3 U = GEO::Geom::perpendicular(N);
        U = normalize(U);
        GEO::vec3 V = cross(N, U);
        V = normalize(V);

        for(GEO::index_t k = 0; k < Co3NeTools::sincos_nb; ++k) {
            double s = Co3NeTools::manifold_table(k, 0);
            double c = Co3NeTools::manifold_table(k, 1);
            GEO::vec3 p = pi + c * radius_ * U + s * radius_ * V;
            P.add_vertex(p);
        }
    }

    void Co3NeRestrictedVoronoiDiagram::get_RVC(GEO::index_t i, const GEO::vec3 &N, GEO_BASE::Polygon &P,
                                                GEO_BASE::Polygon &Q, GEO::vector<GEO::index_t> &neighbor,
                                                GEO::vector<double> &squared_dist) const {
        get_circle(i, P, N);

        GEO::index_t nb_neigh = std::min(GEO::index_t(nb_points() - 1), GEO::index_t(20));
        GEO::index_t jj = 0;

        // just in case, limit to 1000 neighbors.
        GEO::index_t max_neigh = std::min(GEO::index_t(1000), nb_points() - 1);

        while(nb_neigh < max_neigh) {
            if(P.nb_vertices() < 3) {
                return;
            }
            if(neighbor.size() < nb_neigh) {
                get_neighbors(i, neighbor, squared_dist, nb_neigh);
            }
            while(jj < nb_neigh && squared_dist[jj] < 1e-30) {
                jj++;
            }
            while(jj < nb_neigh) {
                if(squared_dist[jj] > sqROS_) {
                    return;
                }
                GEO::index_t j = neighbor[jj];
                double Rk = Co3NeTools::squared_radius(point(i), P);
                if(squared_dist[jj] > 4.0 * Rk) {
                    return;
                }
                Co3NeTools::clip_polygon_by_bisector(P, Q, point(i), point(j), j);
                jj++;
            }
            if(nb_neigh > 3) {
                nb_neigh += nb_neigh / 3;
            } else {
                nb_neigh++;
            }
            nb_neigh = std::min(nb_neigh, nb_points()-1);
        }
    }

    void Co3NeRestrictedVoronoiDiagram::get_RVC(GEO::index_t i, Eigen::MatrixXd tangent_basis,
                                                Polygon_hd& P,Polygon_hd& Q,
                                                GEO::vector<GEO::index_t>& neighbor,
                                                GEO::vector<double>& squared_dist){
        Eigen::MatrixXd temp = Co3NeTools::manifold_table *tangent_basis.transpose();

        Eigen::MatrixXd local_table = Co3NeTools::manifold_table*radius_;
        local_table = local_table*tangent_basis.transpose();
        Eigen::VectorXd center = point_hd(i);
        local_table = local_table.rowwise() + center.transpose();

        P.clear();
        for(GEO::index_t k = 0; k < local_table.rows(); k++) {
            Eigen::VectorXd p = local_table.row(k);
            Vertex_hd V(p);
            P.add_vertex(V);
        }

        GEO::index_t nb_neigh = std::min(GEO::index_t(nb_points() - 1), GEO::index_t(20));
        GEO::index_t jj = 0;

        // just in case, limit to 1000 neighbors.
        GEO::index_t max_neigh = std::min(GEO::index_t(1000), nb_points() - 1);

        while(nb_neigh < max_neigh) {
            if(P.nb_vertices() < 3) {
                return;
            }
            if(neighbor.size() < nb_neigh) {
                get_neighbors(i, neighbor, squared_dist, nb_neigh);
            }
            while(jj < nb_neigh && squared_dist[jj] < 1e-30) {
                jj++;
            }
            while(jj < nb_neigh) {

//                if(squared_dist[jj] > sqROS_) {
//                    return;
//                }
                GEO::index_t j = neighbor[jj];
                double Rk = Co3NeTools::squared_radius_hd(point_hd(i), P);
                if(squared_dist[jj] > 4.0 * Rk) {
                    return;
                }
                Co3NeTools::clip_polygon_by_bisector_hd(P, Q, point_hd(i), point_hd(j), j);
                jj++;
            }
            if(nb_neigh > 3) {
                nb_neigh += nb_neigh / 3;
            } else {
                nb_neigh++;
            }
            nb_neigh = std::min(nb_neigh, nb_points()-1);
        }
    }

    void Co3NeRestrictedVoronoiDiagram::get_neighbor_by_3d(GEO::index_t i,
                                                           GEO::vector<GEO::index_t>& neigh,
                                                           GEO::vector<double>& sq_dist, GEO::index_t nb) {
        GEO::vector<double> sq_dist_3d(nb);
        if(sq_dist.size() < nb) {
            sq_dist.resize(nb);
        }
        if(neigh.size() < nb) {
            neigh.resize(nb);
        }

        NN3d_->get_nearest_neighbors(nb,i,neigh.data(),sq_dist_3d.data());

        for(GEO::index_t j = 0; j < nb; ++j) {
            GEO::index_t k = neigh[j];
            double* p_hd = p_ + k * p_stride_;
            double* pi_hd = p_ + i * p_stride_;

            sq_dist[j] = Co3NeTools::L2_distance(p_hd, pi_hd, p_stride_);
        }

    }
}