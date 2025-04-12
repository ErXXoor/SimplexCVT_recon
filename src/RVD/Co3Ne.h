//
// Created by hongbo on 09/04/25.
//

#ifndef SIMPLEXCVT_RECON_CO3NE_H
#define SIMPLEXCVT_RECON_CO3NE_H
#include <geogram/mesh/mesh.h>
#include <geogram/basic/process.h>
#include <memory>
#include "RVD/Co3NeRestrictedVoronoiDiagram.h"
namespace GEO_BASE{
    class Co3NeThread;
    class Co3Ne{
    public:
        Co3Ne(GEO::Mesh& M);
        void run_threads();
        void reconstructhd(double r);
        void extract_triangles();

        void save_raw_triangles(const std::string& filename);
        void save_T3_triangles(const std::string& filename);
        void save_T12_triangles(const std::string& filename);

    private:
        GEO::Mesh& mesh_;
        std::vector<double> new_vertices_;
        Co3NeRestrictedVoronoiDiagram RVD_;
        GEO::TypedThreadGroup<Co3NeThread> thread_;
        double min_cos_angle_;

        GEO::vector<GEO::index_t> raw_triangles_;
        GEO::vector<GEO::index_t> T3_triangles_;
        GEO::vector<GEO::index_t> T12_triangles_;

    public:
        enum Co3NeMode {
            CO3NE_NONE,    /**< uninitialized */
//            CO3NE_NORMALS, /**< estimate normals in pointset */
//            CO3NE_SMOOTH,  /**< smooth the pointset */
//            CO3NE_RECONSTRUCT, /**< reconstruct the triangles */
            CO3NE_NORMALS_AND_RECONSTRUCT,
            CO3NE_TANGENT_AND_RECONSTRUCT
            /**< combined normal estimation and reconstruction */
        };

        void set_max_angle(double angle){
            min_cos_angle_ = std::cos(angle);
        }

        Co3NeRestrictedVoronoiDiagram& RVD() {
            return RVD_;
        }
    };
}
#endif //SIMPLEXCVT_RECON_CO3NE_H
