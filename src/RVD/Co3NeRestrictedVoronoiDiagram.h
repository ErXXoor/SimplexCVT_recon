//
// Created by hongbo on 09/04/25.
//

#ifndef SIMPLEXCVT_RECON_CO3NERESTRICTEDVORONOIDIAGRAM_H
#define SIMPLEXCVT_RECON_CO3NERESTRICTEDVORONOIDIAGRAM_H
#include <geogram/basic/geometry.h>
#include <Eigen/Dense>
#include "RVD/Geometry.h"
#include <geogram/points/nn_search.h>
#include <geogram/mesh/mesh.h>
namespace GEO_BASE{
    class Co3NeRestrictedVoronoiDiagram{
    public:
        Co3NeRestrictedVoronoiDiagram();
        void init(GEO::Mesh& M);
        void init(GEO::index_t nb_points_in,
                  double* p, GEO::index_t p_stride,
                  double* n, GEO::index_t n_stride);

        void get_circle(GEO::index_t i, Polygon& P, const GEO::vec3& N) const;
        void get_RVC(
                GEO::index_t i, const GEO::vec3& N, Polygon& P,
                Polygon& Q,
                GEO::vector<GEO::index_t>& neighbor,
                GEO::vector<double>& squared_dist
        ) const;

        void get_RVC(GEO::index_t i, Eigen::MatrixXd tangent_basis,
                     Polygon_hd& P,Polygon_hd& Q,
                     GEO::vector<GEO::index_t>& neighbor,
                     GEO::vector<double>& squared_dist
        );

        void get_neighbor_by_3d(GEO::index_t i,
                                GEO::vector<GEO::index_t>& neigh,
                                GEO::vector<double>& sq_dist,
                                GEO::index_t nb);

        void get_neighbor_3d_dist(GEO::index_t i,
                                  GEO::vector<GEO::index_t>& neigh,
                                  GEO::vector<double>& sq_dist_3d,
                                  GEO::index_t nb);

    private:
        friend class Co3Ne;

        GEO::index_t nb_points_;
        double* p_;
        GEO::index_t p_stride_;
        double* n_;
        GEO::index_t n_stride_;
        double radius_;

        GEO::NearestNeighborSearch_var NN_;


        Eigen::MatrixXd p3d_;
        GEO::NearestNeighborSearch_var NN3d_;

        double sqROS_;
        GEO::index_t nb_neighbors_;

    public:
        void set_exact(bool x) {
            NN_->set_exact(x);
        }

        void set_circles_radius(double r) {
            radius_ = r;
            sqROS_ = 4.0 * radius_ * radius_*0.01*0.4;  // squared radius of security
            // when a neighbor is further away than ROS, then it cannot
            // clip a circle of radius r
        }



        GEO::index_t nb_points() const {
            return nb_points_;
        }

        void clear() {
            NN_.reset();
            nb_points_ = 0;
            p_ = nullptr;
            p_stride_ = 0;
            n_ = nullptr;
            n_stride_ = 0;
            nb_neighbors_ = 0;
        }

        void update(){
            init(nb_points_, p_, p_stride_, n_, n_stride_);
        }

        const GEO::vec3& point(GEO::index_t i) const {
            geo_debug_assert(i < nb_points());
            return *(GEO::vec3*) (p_ + i * p_stride_);
            // Yes I know, this is a bit ugly...
        }

        const double* point_ptr(GEO::index_t i) const {
            geo_debug_assert(i < nb_points());
            return p_ + i * p_stride_;
        }

        Eigen::VectorXd point_hd(GEO::index_t i) const {
            geo_debug_assert(i < nb_points());
            Eigen::VectorXd v(p_stride_);
            for(GEO::index_t j=0; j<p_stride_; ++j) {
                v(j) = p_[i * p_stride_ + j];
            }
            return v;
        }

        const GEO::vec3& normal(GEO::index_t i) const {
            geo_debug_assert(n_ != nullptr);
            geo_debug_assert(i < nb_points());
            return *(GEO::vec3*) (n_ + i * n_stride_);
            // Yes I know, this is a bit ugly...
        }

        void set_normal(GEO::index_t i, const GEO::vec3& N) const {
            geo_debug_assert(n_ != nullptr);
            geo_debug_assert(i < nb_points());
            double* n = n_ + i * n_stride_;
            n[0] = N.x;
            n[1] = N.y;
            n[2] = N.z;
        }
        void get_neighbors(
                GEO::index_t i,
                GEO::index_t* neigh,
                double* sq_dist,
                GEO::index_t nb
        ) const {
            return NN_->get_nearest_neighbors(
                    nb, i, neigh, sq_dist
            );
        }

        void get_neighbors(
                GEO::index_t i,
                GEO::vector<GEO::index_t>& neigh,
                GEO::vector<double>& sq_dist,
                GEO::index_t nb
        ) const {
            neigh.resize(nb);
            sq_dist.resize(nb);
            get_neighbors(i, neigh.data(), sq_dist.data(), nb);
        }

        void get_RVC(
                GEO::index_t i, Polygon& P,
                Polygon& Q,
                GEO::vector<GEO::index_t>& neighbor,
                GEO::vector<double>& squared_dist
        ) const {
            neighbor.resize(0);
            squared_dist.resize(0);
            get_RVC(i, normal(i), P, Q, neighbor, squared_dist);
        }

        GEO::index_t nb_neighbors() const {
            return std::min(nb_neighbors_,nb_points()-1);
        }

        void set_nb_neighbors(GEO::index_t x) {
            nb_neighbors_ = x;
        }

        ~Co3NeRestrictedVoronoiDiagram() {
            clear();
        }
    };
}


#endif //SIMPLEXCVT_RECON_CO3NERESTRICTEDVORONOIDIAGRAM_H
