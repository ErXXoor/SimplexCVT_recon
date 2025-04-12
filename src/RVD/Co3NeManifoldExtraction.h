//
// Created by hongbo on 11/04/25.
//

#ifndef SIMPLEXCVT_RECON_CO3NEMANIFOLDEXTRACTION_H
#define SIMPLEXCVT_RECON_CO3NEMANIFOLDEXTRACTION_H
#include <geogram/mesh/mesh.h>
namespace GEO_BASE {
    class Co3NeManifoldExtraction {
    public:
        static constexpr GEO::index_t NO_CORNER = GEO::index_t(-1);
        static constexpr GEO::index_t NO_FACET = GEO::index_t(-1);
        static constexpr GEO::index_t NO_CNX = GEO::index_t(-1);

        /**
         * \brief Initializes a new Co3NeManifoldExtraction with
         *  a list of triangles.
         * \param[in,out] target the target mesh. It needs to be already
         *  initialized with the vertices.
         * \param[in,out] good_triangles the good triangles reconstructed
         *  by Co3Ne. They are 'stealed' by the mesh (on exit, good_triangles
         *  is empty). If some non-manifold edges are detected, then all
         *  the triangles incident to any manifold edge are ignored.
         */
        Co3NeManifoldExtraction(
                GEO::Mesh &target,
                GEO::vector <GEO::index_t> &good_triangles,
                bool strict = false
        );

        /**
         * \brief Tentatively adds triangle from the specified list.
         * \details Some geometric and topological properties are
         *  verified by connect_and_validate_triangle() before accepting
         *  the triangle.
         * \see connect_and_validate_triangle()
         */
        void add_triangles(const GEO::vector <GEO::index_t> &not_so_good_triangles);

    protected:

        /**
         * \brief Initializes the combinatorial data
         *  structures and deletes all facets incident
         *  to a non-manifold edge.
         */
        void init_and_remove_non_manifold_edges();

        /**
         * \brief Tentatively connects a newly added triangle
         *  to the current mesh under construction. Accepted
         *  triangles satisfy the following criteria:
         *  - each new triangle should be either incident to at least
         *    two edges of existing triangles, or to one existing triangle
         *    and one isolated point.
         *  - the normals to the new triangle and its neighbor should
         *    not point to opposite directions.
         *  - inserting the new triangle should not generate 'by-excess'
         *    non-manifold vertices. A 'by-excess' non-manifold vertex
         *    has a closed loop of triangles in its neighbors plus
         *    additional triangles.
         *  - the orientation of the surface should be coherent (no Moebius
         *    strip).
         * \param[in] t index of the triangle
         * \param[out] classified true if the status of the triangle
         *  (accepted/rejected) could be completely determined,
         *  false if its status may still change during subsequent iterations
         * \retval true if all combinatorial and geometric tests succeeded
         * \retval false otherwise
         */
        bool connect_and_validate_triangle(GEO::index_t t, bool &classified);


        /**
         * \brief Tentatively enforces mesh orientation starting from a
         *  given triangle.
         * \details The triangle \p t is rejected if it is incident to
         *  the same connected component with two different orientations.
         * \param[in] t index of the triangle to start mesh orientation from
         * \retval true if the mesh could be coherently oriented
         * \retval false otherwise
         */
        bool enforce_orientation_from_triangle(GEO::index_t t);


        /**
         * \brief Adds a new triangle to the surface and to the
         *  combinatorial data structure.
         * \param[in] i first index of the triangle
         * \param[in] j second index of the triangle
         * \param[in] k third index of the triangle
         */
        GEO::index_t add_triangle(GEO::index_t i, GEO::index_t j, GEO::index_t k);

        /**
         * \brief Removes the latest triangle from both
         *  the mesh and the combinatorial data structure.
         */
        void rollback_triangle() {
            GEO::index_t t = M_.facets.nb() - 1;
            remove(t);
            M_.facets.pop();
        }


        /**
         * \brief Inverts the orientation of a triangle.
         * \param[in] t the index of the triangle to be flipped.
         */
        void flip_triangle(GEO::index_t t);

        /**
         * \brief Inserts a triangle of the mesh into the data structures
         *  used for topology checks.
         * \param[in] t index of the triangles to be inserted
         * \pre \p t is a valid triangle index in the mesh
         */
        void insert(GEO::index_t t);

        /**
         * \brief Removes a triangle of the mesh from the data structures
         *  used for topology/combinatorial checks.
         * \param[in] t index of the triangles to be removed
         * \param[in] disconnect if true, connections from the neighbors
         *  to t are set to -1 (facet_corners.adjacent_facet).
         * \pre \p t is a valid triangle index in the mesh
         */
        void remove(GEO::index_t t, bool disconnect = true);

        /**
         * \brief Gets the number of triangles incident
         *  to a vertex.
         * \param[in] v index of the vertex
         * \return the number of triangles incident to \p v
         */
        GEO::index_t nb_incident_triangles(GEO::index_t v) const {
            GEO::index_t result = 0;
            GEO::index_t c = v2c_[v];
            do {
                ++result;
                c = next_c_around_v_[c];
            } while (c != v2c_[v]);
            return result;
        }

        /**
         * \brief Tests whether a given vertex is non-manifold
         *  by excess.
         * \details A vertex is non-manifold by-excess if its
         *  set of incident triangles contains a closed loop
         *  of triangles and additional triangles.
         * \param[in] v index of the vertex to be tested
         * \retval true if \p v is non-manifold by excess
         * \retval false otherwise
         */
        bool vertex_is_non_manifold_by_excess(GEO::index_t v, bool &moebius);

        /**
         * \brief Gets the next corner around a vertex from a given
         *  corner.
         * \details This function works even for a mesh that has triangles
         *  that are not coherently oriented. In other words, for two
         *  corners c1, c2, if we have:
         *   - v1 = facet_corners.vertex(c1)
         *   - v2 = facet_corners.vertex(
         *       c1,facets.next_corner_around_facet(c2f(c1),c1)
         *   )
         *   - w1 = facet_corners.vertex(c2)
         *   - w2 = facet_corners.vertex(
         *         c2,facets.next_corner_around_facet(c2f(c2),c2)
         *   )
         *  then we can have:
         *   - v1=w2 and v2=w1 (as usual) or:
         *   - v1=v2 and w1=w2 ('inverted' configuration)
         * \param[in] v the vertex
         * \param[in] c1 a corner incident to \p v or pointing to \p v
         * \return another corner incident to the \p v
         */
        GEO::index_t next_around_vertex_unoriented(
                GEO::index_t v, GEO::index_t c1
        ) const;

        /**
         * \brief Gets the three corners adjacent to a triangle.
         * \details This function works even for a mesh that has triangles
         *  that are not coherently oriented. In other words, for two
         *  corners c1, c2, if we have:
         *   - v1 = facet_corners.vertex(c1)
         *   - v2 = facet_corners.vertex(
         *       c1,facets.next_corner_around_facet(c2f(c1),c1))
         *   - w1 = facet_corners.vertex(c2)
         *   - w2 = facet_corners.vertex(
         *       c2,facets.next_corner_around_facet(c2f(c2),c2))
         *  then c1 and c2 are adjacent if we have:
         *   - v1=w2 and v2=w1 (as usual) or:
         *   - v1=v2 and w1=w2 ('inverted' configuration)
         * \param[in] t1 index of the triangle
         * \param[out] adj_c index of the adjacent corners
         *  (array of 3 integers). Each entry contains a valid corner index
         *  or NO_CORNER if the corresponding edge is on the border.
         * \retval true if the three edges are manifold
         * \retval false otherwise (and then \p adj_c contains undefined
         *  values).
         */
        bool get_adjacent_corners(GEO::index_t t1, GEO::index_t *adj_c);

        /**
         * \brief Tentatively connect a triangle of the mesh with its
         *   neighbors.
         * \details This function is independent of triangles orientations,
         *  see get_adjacent_corners().
         * \param[in] t index of the triangle to be connected
         * \param[in] adj_c an array of three integers that indicate
         *  for each corner of the triangle the index of the adjacent
         *  corner or NO_CORNER if the corner is on the border.
         */
        void connect_adjacent_corners(GEO::index_t t, GEO::index_t *adj_c) {
            for (GEO::index_t i = 0; i < 3; ++i) {
                if (adj_c[i] != NO_CORNER) {
                    GEO::index_t c = M_.facets.corners_begin(t) + i;
                    M_.facet_corners.set_adjacent_facet(c, c2f(adj_c[i]));
                    M_.facet_corners.set_adjacent_facet(adj_c[i], t);
                }
            }
        }

        /**
         * \brief Tentatively connect a triangle of the mesh with its
         *   neighbors.
         * \details This function is independent of triangles orientations,
         *  see get_adjacent_corners().
         * \param[in] t index of the triangle to be connected
         * \retval false if the connection would have created non-manifold
         *  edges
         * \retval true otherwise
         */
        bool connect(GEO::index_t t) {
            GEO::index_t adj_c[3];
            if (!get_adjacent_corners(t, adj_c)) {
                return false;
            }
            connect_adjacent_corners(t, adj_c);
            return true;
        }

        /**
         * \brief Gets a facet index by corner index.
         * \details for a triangulated mesh, indexing is
         *  implicit, and we do not need to store a c2f array.
         * \param[in] c corner index
         * \return the index of the facet incident to c
         */
        GEO::index_t c2f(GEO::index_t c) const {
            geo_debug_assert(c != NO_CORNER);
            geo_debug_assert(c < M_.facet_corners.nb());
            return c / 3;
        }


        /**
         * \brief Tests whether two triangles have the
         *  same orientation.
         * \param[in] t1 first triangle
         * \param[in] t2 second triangle
         * \retval true if \p t1 and \p t2 have the same
         *  orientation
         * \retval false otherwise
         * \pre \p t1 and \p t2 share an edge
         */
        bool triangles_have_same_orientation(
                GEO::index_t t1,
                GEO::index_t t2
        );


        /**
         * \brief Tests whether the normals of two triangles that
         *  share an edge 'agree', i.e. whether they do not form
         *  a too sharp angle.
         * \param[in] t1 index of the first triangle
         * \param[in] t2 index of the second triangle
         * \retval true if the normals of both triangles do not
         *  point in opposite directions
         * \retval false otherwise
         * \pre the two triangles are incident to the same edge
         *  (they have two vertices in common)
         */
        bool triangles_normals_agree(
                GEO::index_t t1,
                GEO::index_t t2
        ) const;

        /**
         * \brief Merges two connected components.
         * \details The connected component incident to \p t
         *  is replaced with \p comp2.
         * \param [in] t index of a triangle incident
         *  to the first connected component
         * \param [in] comp2 index of the second connected
         *  component
         * \param [in] flip if true, flip the triangles
         * \pre At least one of the triangles adjacent to
         *  \p t (directly or not) is incident to
         *  component \p comp2
         */
        void merge_connected_component(
                GEO::index_t t,
                GEO::index_t comp2,
                bool flip
        );

        /**
         * \brief Initializes the date structures
         *  that represent the connected components.
         * \details This function computes cnx_ and
         *  cnx_size_. The array cnx_[f] gives for each
         *  facet f the index of the connected component
         *  that contains f, and the array cnx_size_[comp]
         *  gives for each connected component comp the
         *  number of facets in comp.
         */
        void init_connected_components();

    private:
        GEO::Mesh &M_;

        /**
         * \brief For each corner, next_c_around_v_[c]
         * chains the circular list of corners
         * incident to the same corner as c.
         */
        GEO::vector <GEO::index_t> next_c_around_v_;

        /**
         * \brief For each vertex v, v2c_[v] contains a
         *  corner incident to v, or NO_VERTEX if v is
         *  isolated.
         */
        GEO::vector <GEO::index_t> v2c_;


        /**
         * \brief For each triangle t, cnx_[t] contains
         *  the index of the connected component of the
         *  mesh incident to t.
         */
        GEO::vector <GEO::index_t> cnx_;

        /**
         * \brief For each connected component C,
         *  cnx_size_[C] contains the number of
         *  facets in C.
         */
        GEO::vector <GEO::index_t> cnx_size_;

        /**
         * \brief In strict mode, each inserted triangle
         *  is checked for non-manifold configuration.
         *  In non-strict mode, only T2 and T1 triangles are
         *  tested (those seen from only 2 or only 1 Voronoi
         *  cell), T3 triangles are inserted without test.
         */
        bool strict_;
    };
}
#endif //SIMPLEXCVT_RECON_CO3NEMANIFOLDEXTRACTION_H
