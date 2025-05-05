//
// Created by hongbo on 09/04/25.
//
#include "Base/MeshAdaptor.h"
#include <geogram/basic/line_stream.h>
#include <geogram/mesh/mesh_io.h>

namespace Base{
    void MeshAdaptor::LoadXYZ(const std::string &filepath, GEO::Mesh &mesh,int dim) {
        mesh.clear();
        mesh.vertices.set_double_precision();
        mesh.vertices.set_dimension(dim);
        GEO::LineInput in(filepath);
        if (!in.OK()) {
            std::cerr << "Error: cannot open file " << filepath << std::endl;
            return;
        }

        std::vector<double> P(dim);
        while (!in.eof() && in.get_line()) {
            in.get_fields();
            if (in.nb_fields() >= 1) {
                for (GEO::coord_index_t c = 0; c < dim; c++) {
                    if (GEO::index_t(c) < in.nb_fields()) {
                        P[c] = in.field_as_double(GEO::index_t(c));
                    } else {
                        P[c] = 0.0;
                    }

                }
                GEO::index_t v = mesh.vertices.create_vertex();
                double *p = mesh.vertices.point_ptr(v);
                for (GEO::index_t c = 0; c < dim; ++c) {
                    p[c] = P[c];
                }
            }
        }
    }

    void MeshAdaptor::SaveMesh(GEO::Mesh& mesh,
                               GEO::vector<GEO::index_t> triangles,
                               const std::string &filepath) {
        GEO::Mesh M;
        M.vertices.assign_points(
                mesh.vertices.point_ptr(0),
                mesh.vertices.dimension(),
                mesh.vertices.nb()
        );
        M.facets.assign_triangle_mesh(triangles, false);
        M.vertices.set_dimension(3);
        GEO::mesh_save(M, filepath);
    }

    void MeshAdaptor::SaveMesh(GEO::Mesh& mesh,
                               const std::string &filepath) {
        GEO::MeshIOFlags flags;
        flags.set_attribute(GEO::MESH_ALL_ATTRIBUTES);
        flags.set_dimension(3);
        GEO::mesh_save(mesh, filepath, flags);
    }
}