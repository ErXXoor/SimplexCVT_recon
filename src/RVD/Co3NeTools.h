//
// Created by hongbo on 10/04/25.
//

#ifndef SIMPLEXCVT_RECON_CO3NETOOLS_H
#define SIMPLEXCVT_RECON_CO3NETOOLS_H
#include <geogram/basic/numeric.h>
#include <vector>
#include "RVD/Geometry.h"
#include <Eigen/Dense>
namespace GEO_BASE{
    class Co3NeTools{
    public:
        static void co3ne_split_triangles_list(
                GEO::vector<GEO::index_t>& triangles,
                GEO::vector<GEO::index_t>& good_triangles,
                GEO::vector<GEO::index_t>& not_so_good_triangles
        );

        static void clip_polygon_by_bisector(
                Polygon& Ping, Polygon& Pong,
                const GEO::vec3& pi, const GEO::vec3& pj, GEO::index_t j
        );

        static void clip_polygon_by_bisector_hd(
                Polygon_hd& Ping,Polygon_hd& Pong,
                const Eigen::VectorXd& pi, const Eigen::VectorXd& pj, GEO::index_t j);

        static double squared_radius(const GEO::vec3& p, const Polygon& P);
        static double squared_radius_hd(const Eigen::VectorXd& p, const Polygon_hd& P);

        static Eigen::MatrixXd PCA_fit(Eigen::MatrixXd &points);

        static double L2_distance(const double* p1, const double* p2, int dim);

        static Eigen::MatrixXd GenCircleTable(int numPoints);
        static constexpr GEO::index_t sincos_nb = 10;
        static Eigen::MatrixXd manifold_table;
    };
}
#endif //SIMPLEXCVT_RECON_CO3NETOOLS_H
