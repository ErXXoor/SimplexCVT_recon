//
// Created by hongbo on 09/04/25.
//

#ifndef SIMPLEXCVT_RECON_MESHADAPTOR_H
#define SIMPLEXCVT_RECON_MESHADAPTOR_H
#include <geogram/mesh/mesh.h>

namespace Base{
    class MeshAdaptor{
    public:
        static void LoadHDXYZ(const std::string &filepath,
                              GEO::Mesh &mesh,
                              int dim);

        static void SaveMesh(GEO::Mesh &mesh,
                             GEO::vector<GEO::index_t> triangles,
                             const std::string &filepath);
    };

}

#endif //SIMPLEXCVT_RECON_MESHADAPTOR_H
