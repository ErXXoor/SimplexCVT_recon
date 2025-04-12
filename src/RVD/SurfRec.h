//
// Created by hongbo on 09/04/25.
//

#ifndef SIMPLEXCVT_RECON_RECONSTRUCTOR_H
#define SIMPLEXCVT_RECON_RECONSTRUCTOR_H
#include <geogram/mesh/mesh.h>
#include "RVD/Co3Ne.h"
namespace RVD{
    class SurfRec{
    public:
        SurfRec() = default;
        ~SurfRec() = default;

        void init(GEO::Mesh& point_in);
        void Co3ne_rec(std::string output_file);

    private:
        double m_radius;
        std::shared_ptr<GEO_BASE::Co3Ne> m_co3ne;
    };

}
#endif //SIMPLEXCVT_RECON_RECONSTRUCTOR_H
