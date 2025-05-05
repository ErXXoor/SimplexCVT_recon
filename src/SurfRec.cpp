//
// Created by hongbo on 09/04/25.
//
#include "RVD/SurfRec.h"
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/points/nn_search.h>
namespace RVD{
    void SurfRec::init(GEO::Mesh &point_in) {
//        //get avg distance of 60 nn
//        GEO::index_t k = 100;
//        GEO::NearestNeighborSearch_var nn_search;
//        nn_search = GEO::NearestNeighborSearch::create(point_in.vertices.dimension());
//        nn_search->set_points(point_in.vertices.nb(),
//                              point_in.vertices.point_ptr(0),
//                              point_in.vertices.dimension());
//        double total_distance = 0;
//        std::vector<double> distances;
//        std::vector<GEO::index_t> nn_points;
//        for ( auto i = 0; i < point_in.vertices.nb(); ++i) {
//            distances.clear();
//            distances.resize(k);
//
//            nn_points.clear();
//            nn_points.resize(k);
//
//            nn_search->get_nearest_neighbors(k,i,&nn_points[0],&distances[0]);
//            total_distance += distances[k-1];
//        }
//        m_radius = total_distance/point_in.vertices.nb();

        auto R = GEO::bbox_diagonal(point_in);
        m_radius = 0.05*R;

        m_co3ne = std::make_shared<GEO_BASE::Co3Ne>(point_in);
    }

    void SurfRec::Co3ne_rec3d(const std::string &output_file) {
       m_co3ne->reconstruct3d(m_radius);
    }

    void SurfRec::Co3ne_rec(const std::string& output_file,
                            int dim,
                            bool post_process) {
        if(dim==3){
            m_co3ne->reconstruct3d(m_radius);
        }
        else{
            m_co3ne->reconstructhd(m_radius);
        }

        m_co3ne->extract_triangles();

        if(!post_process){
            m_co3ne->save_raw_triangles(output_file);
        }
        else{
            m_co3ne->post_process();
            m_co3ne->save_final_mesh(output_file);
        }
    }
}