//
// Created by hongbo on 09/04/25.
//
#include "RVD/SurfRec.h"
#include <geogram/mesh/mesh_geometry.h>
namespace RVD{
    void SurfRec::init(GEO::Mesh &point_in) {
        auto R = GEO::bbox_diagonal(point_in);
        m_radius = 5*0.01*R;

        m_co3ne = std::make_shared<GEO_BASE::Co3Ne>(point_in);
    }

    void SurfRec::Co3ne_rec(std::string output_file) {
        m_co3ne->reconstructhd(m_radius);
        m_co3ne->extract_triangles();
        m_co3ne->save_raw_triangles(output_file);
    }
}