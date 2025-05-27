//
// Created by hongbo on 10/04/25.
//
#include "RVD/Co3NeTools.h"
#include "RVD/CompareTriangles.h"
namespace GEO_BASE{
    Eigen::MatrixXd Co3NeTools::GenCircleTable(int numPoints){
        Eigen::MatrixXd circle_table(numPoints, 2);
        for(int i=0; i<numPoints; i++){
            double angle = 2 * M_PI * i / (numPoints-1);
            circle_table(i, 0) = sin(angle);
            circle_table(i, 1) = cos(angle);
        }
        return circle_table;
    }

    Eigen::MatrixXd Co3NeTools::manifold_table = Co3NeTools::GenCircleTable(sincos_nb);

    void Co3NeTools::co3ne_split_triangles_list(GEO::vector<GEO::index_t> &triangles,
                                                GEO::vector<GEO::index_t> &good_triangles,
                                                GEO::vector<GEO::index_t> &not_so_good_triangles) {
        GEO::index_t nb_triangles = triangles.size()/3;

        // Step 1: normalize vertices order
        for(GEO::index_t i=0; i<triangles.size(); i+=3) {
            GEO::index_t* ptr = &triangles[i];
            std::sort(ptr, ptr+3);
        }

        // Step 2: sort the triangles in lexicographic order
        GEO::vector<GEO::index_t> t_sort(nb_triangles);
        for(GEO::index_t t=0; t<nb_triangles; ++t) {
            t_sort[t] = t ;
        }
        CompareTriangles compare_triangles(triangles);
        sort(t_sort.begin(), t_sort.end(), compare_triangles);


        // Step 3: select the triangles that appear exactly 3 times
        GEO::index_t if1 = 0;
        while(if1 < nb_triangles) {
            GEO::index_t if2 = if1 + 1;
            while(
                    if2 < nb_triangles &&
                    compare_triangles.is_same(t_sort[if1], t_sort[if2])
                    ) {
                if2++;
            }

            GEO::index_t t = t_sort[if1];
            if(if2 - if1 == 3) {
                good_triangles.push_back(triangles[3*t]);
                good_triangles.push_back(triangles[3*t+1]);
                good_triangles.push_back(triangles[3*t+2]);
            } else if(if2 - if1 <= 2) {
                not_so_good_triangles.push_back(triangles[3*t]);
                not_so_good_triangles.push_back(triangles[3*t+1]);
                not_so_good_triangles.push_back(triangles[3*t+2]);
            }
            if1 = if2;
        }
    }

    void Co3NeTools::clip_polygon_by_bisector(GEO_BASE::Polygon &Ping, GEO_BASE::Polygon &Pong, const GEO::vec3 &pi,
                                              const GEO::vec3 &pj, GEO::index_t j) {
        if(Ping.nb_vertices() == 0) {
            return;
        }
        Pong.clear();

        GEO::vec3 n(
                pi.x - pj.x,
                pi.y - pj.y,
                pi.z - pj.z
        );

        // Compute d = n . m, where n is the
        // normal vector of the bisector [pi,pj]
        // and m twice the middle point of the bisector.
        double d =
                n.x * (pi.x + pj.x) +
                n.y * (pi.y + pj.y) +
                n.z * (pi.z + pj.z);

        // The predecessor of the first vertex is the last vertex
        GEO::index_t prev_k = Ping.nb_vertices() - 1;
        const Vertex* prev_vk = &(Ping.vertex(prev_k));

        // We compute:
        //    prev_l = prev_vk . n
        double prev_l = dot(prev_vk->point(), n);

        // We compute:
        //    side1(pi,pj,q) = sign(2*q.n - n.m) = sign(2*l - d)
        GEO::Sign prev_status = GEO::geo_sgn(2.0 * prev_l - d);

        for(GEO::index_t k = 0; k < Ping.nb_vertices(); k++) {
            const Vertex* vk = &(Ping.vertex(k));

            // We compute: l = vk . n
            double l = dot(vk->point(), n);

            // We compute:
            //   side1(pi,pj,q) = sign(2*q.n - n.m) = sign(2*l - d)
            GEO::Sign status = GEO::geo_sgn(2.0 * l - d);

            // If status of edge extremities differ,
            // then there is an intersection.
            if(status != prev_status && (prev_status != 0)) {

                // Compute lambda1 and lambda2, the
                // barycentric coordinates of the intersection I
                // in the segment [prev_vk vk]
                // Note that d and l (used for the predicates)
                // are reused here.
                double denom = 2.0 * (prev_l - l);
                double lambda1, lambda2;

                // Shit happens ! [Forrest Gump]
                if(::fabs(denom) < 1e-20) {
                    lambda1 = 0.5;
                    lambda2 = 0.5;
                } else {
                    lambda1 = (d - 2.0 * l) / denom;
                    // Note: lambda2 is also given
                    // by (2.0*l2-d)/denom
                    // (but 1.0 - lambda1 is a bit
                    //  faster to compute...)
                    lambda2 = 1.0 - lambda1;
                }
                Vertex V;
                V.point().x =
                        lambda1 * prev_vk->point().x + lambda2 * vk->point().x;
                V.point().y =
                        lambda1 * prev_vk->point().y + lambda2 * vk->point().y;
                V.point().z =
                        lambda1 * prev_vk->point().z + lambda2 * vk->point().z;
                if(status > 0) {
                    V.set_adjacent_seed(prev_vk->adjacent_seed());
                } else {
                    V.set_adjacent_seed(GEO::signed_index_t(j));
                }
                Pong.add_vertex(V);
            }
            if(status > 0) {
                Pong.add_vertex(*vk);
            }
            prev_vk = vk;
            prev_status = status;
            prev_k = k;
            prev_l = l;
        }
        Ping.swap(Pong);
    }

    void Co3NeTools::clip_polygon_by_bisector_hd(GEO_BASE::Polygon_hd &Ping, GEO_BASE::Polygon_hd &Pong,
                                                 const Eigen::VectorXd &pi, const Eigen::VectorXd &pj, GEO::index_t j) {
        if (Ping.nb_vertices()==0){
            return;
        }
        Pong.clear();
        Eigen::VectorXd n = pi - pj;
        double d = n.dot((pi + pj));
        GEO::index_t prev_k = Ping.nb_vertices() - 1;
        const Vertex_hd* prev_vk = &(Ping.vertex(prev_k));
        double prev_l = prev_vk->point().dot(n);
        GEO::Sign prev_status = GEO::geo_sgn(2 * prev_l - d);

        for(GEO::index_t k=0;k<Ping.nb_vertices();k++) {
            const Vertex_hd* vk = &(Ping.vertex(k));
            double l = vk->point().dot(n);
            GEO::Sign status = GEO::geo_sgn(2.0 * l - d);

            if (status != prev_status && (prev_status != 0)) {
                double denom = 2.0 * (prev_l - l);
                double lambda1, lambda2;
                if (::fabs(denom) < 1e-20) {
                    lambda1 = 0.5;
                    lambda2 = 0.5;
                } else {
                    lambda1 = (d - 2.0 * l) / denom;
                    lambda2 = 1.0 - lambda1;
                }
                Vertex_hd V;
                V.point() = lambda1 * prev_vk->point() + lambda2 * vk->point();
                if (status > 0) {
                    V.set_adjacent_seed(prev_vk->adjacent_seed());
                } else {
                    V.set_adjacent_seed(GEO::signed_index_t(j));
                }
                Pong.add_vertex(V);
            }
            if (status > 0) {
                Pong.add_vertex(*vk);
            }
            prev_vk = vk;
            prev_status = status;
            prev_k = k;
            prev_l = l;
        }
        Ping.swap(Pong);
    }

    Eigen::MatrixXd Co3NeTools::PCA_fit(Eigen::MatrixXd &points) {
        int numPoints = points.rows();

        Eigen::VectorXd mean = points.colwise().mean();

        Eigen::MatrixXd centered = points.rowwise() - mean.transpose();

        Eigen::MatrixXd covariance = (centered.transpose() * centered) / numPoints;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(covariance);
        Eigen::MatrixXd full_eigen = solver.eigenvectors();

        return full_eigen.rightCols(2);
    }

    double Co3NeTools::L2_distance(const double *p1, const double *p2, int dim) {
        Eigen::Map<const Eigen::VectorXd> v1(p1, dim);
        Eigen::Map<const Eigen::VectorXd> v2(p2, dim);
        return (v1 - v2).norm();
    }

    double Co3NeTools::squared_radius(const GEO::vec3& p, const Polygon& P) {
        double result = 0.0;
        for(GEO::index_t i = 0; i < P.nb_vertices(); i++) {
            result = std::max(result, distance2(p, P.vertex(i).point()));
        }
        return result;
    }

    double Co3NeTools::squared_radius_hd(const Eigen::VectorXd& p, const Polygon_hd& P) {
        double result = 0.0;
        for(GEO::index_t i = 0; i < P.nb_vertices(); i++) {
            result = std::max(result, (p - P.vertex(i).point()).squaredNorm());
        }
        return result;
    }


}