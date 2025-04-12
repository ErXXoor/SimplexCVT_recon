//
// Created by hongbo on 10/04/25.
//

#ifndef SIMPLEXCVT_RECON_COMPARETRIANGLES_H
#define SIMPLEXCVT_RECON_COMPARETRIANGLES_H
#include <geogram/basic/numeric.h>
#include <geogram/basic/geometry.h>
namespace GEO_BASE
{
    class CompareTriangles {
    public:
        /**
         * \brief Constructs a new CompareFacets.
         * \param[in] triangles a const reference to a vector
         *  of indices triplets
         */
        explicit CompareTriangles(const GEO::vector<GEO::index_t>& triangles) :
                triangles_(triangles) {
        }

        /**
         * \brief Tests the lexicographic order of two facets by their indices.
         * \param[in] f1 index of the first facet
         * \param[in] f2 index of the second facet
         * \return true if facet \p f1 is before facet \p f2 according to
         *  the lexicographic order of its vertices, false otherwise.
         */
        bool is_before(GEO::index_t f1, GEO::index_t f2) const {
            for(GEO::index_t c=0; c<3; c++) {
                GEO::index_t v1 = triangles_[3*f1+c];
                GEO::index_t v2 = triangles_[3*f2+c];
                if(v1 > v2) {
                    return false;
                }
                if(v1 < v2) {
                    return true;
                }
            }
            return false;
        }

        /**
         * \brief Tests whether two facets are identical.
         * \param[in] f1 index of the first facet
         * \param[in] f2 index of the second facet
         * \return true if facets \p f1 and \p f2 have the same
         *  vertices, false otherwise
         */
        bool is_same(GEO::index_t f1, GEO::index_t f2) const {
            for(GEO::index_t c=0; c<3; c++) {
                GEO::index_t v1 = triangles_[3*f1+c];
                GEO::index_t v2 = triangles_[3*f2+c];
                if(v1 != v2) {
                    return false;
                }
            }
            return true;
        }

        /**
         * \brief Tests the lexicographic order of two facets by their indices.
         * \param[in] f1 index of the first facet
         * \param[in] f2 index of the second facet
         * \return true if facet \p f1 is before facet \p f2 according to
         *  the lexicographic order of its vertices, false otherwise.
         */
        bool operator() (GEO::index_t f1, GEO::index_t f2) const {
            return is_before(f1, f2);
        }

    private:
        const GEO::vector<GEO::index_t>& triangles_;
    };
}
#endif //SIMPLEXCVT_RECON_COMPARETRIANGLES_H
