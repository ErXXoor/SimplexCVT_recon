//
// Created by hongbo on 11/04/25.
//
#include "RVD/Co3NeManifoldExtraction.h"
#include <geogram/basic/command_line.h>
#include <geogram/mesh/mesh_repair.h>
#include <stack>
namespace GEO_BASE{
    Co3NeManifoldExtraction::Co3NeManifoldExtraction(GEO::Mesh &target,
                                                     GEO::vector <GEO::index_t> &good_triangles,
                                                     bool strict): M_(target) {
        if (strict) {
            GEO::vector <GEO::index_t> first_triangle;
            for (GEO::index_t i = 0; i < 3; ++i) {
                first_triangle.push_back(*good_triangles.rbegin());
                good_triangles.pop_back();
            }

            M_.facets.assign_triangle_mesh(first_triangle, true);

            init_and_remove_non_manifold_edges();
            init_connected_components();
            add_triangles(good_triangles);
        } else {

            M_.facets.assign_triangle_mesh(good_triangles, true);

            init_and_remove_non_manifold_edges();
            init_connected_components();
        }
    }

    void Co3NeManifoldExtraction::add_triangles(const GEO::vector<GEO::index_t> &not_so_good_triangles) {

        GEO::index_t nb_triangles = not_so_good_triangles.size() / 3;
        GEO::Logger::out("Co3ne") << "Tentatively add "
                             << nb_triangles << " triangles" << std::endl;
        GEO::vector<bool> t_is_classified(nb_triangles, false);
        bool changed = true;
        GEO::index_t max_iter = strict_ ? 5000 : 50;
        GEO::index_t iter = 0;
        bool first = true;
        while (changed && iter < max_iter) {
            if (first) {
                GEO::CmdLine::ui_clear_line();
            } else {
                first = false;
            }

            GEO::Logger::out("Manifold Rec")
                        << "Iteration:" << iter << std::endl;

            changed = false;
            ++iter;
            for (GEO::index_t t = 0; t < nb_triangles; ++t) {
                if (!t_is_classified[t]) {
                    GEO::index_t i = not_so_good_triangles[3 * t];
                    GEO::index_t j = not_so_good_triangles[3 * t + 1];
                    GEO::index_t k = not_so_good_triangles[3 * t + 2];
                    GEO::index_t new_t = add_triangle(i, j, k);
                    bool classified = false;
                    if (connect_and_validate_triangle(new_t, classified)) {
                        changed = true;
                    } else {
                        rollback_triangle();
                    }
                    if (classified) {
                        t_is_classified[t] = true;
                    }
                }
            }
        }

        GEO::Logger::out("Manifold Rec")
                    << "Iteration:" << iter << std::endl;
    }

    void Co3NeManifoldExtraction::init_and_remove_non_manifold_edges() {
        next_c_around_v_.assign(
                M_.facet_corners.nb(), GEO::index_t(NO_CORNER)
        );
        v2c_.assign(
                M_.vertices.nb(), GEO::index_t(NO_CORNER)
        );
        for (GEO::index_t t = 0; t < M_.facets.nb(); ++t) {
            insert(t);
        }

        GEO::index_t nb_non_manifold = 0;
        GEO::vector <GEO::index_t> remove_t;
        for (GEO::index_t t = 0; t < M_.facets.nb(); ++t) {
            if (!connect(t)) {
                remove_t.resize(M_.facets.nb(), 0);
                remove_t[t] = 1;
                ++nb_non_manifold;
            }
        }

        GEO::mesh_reorient(M_, &remove_t);

        if (remove_t.size() == 0) {
            GEO::Logger::out("Co3Ne")
                    << "All edges are manifold and well oriented"
                    << std::endl;
        } else {
            GEO::index_t nb_remove_t = 0;
            for (GEO::index_t t = 0; t < M_.facets.nb(); ++t) {
                if (remove_t[t] != 0) {
                    ++nb_remove_t;
                }
            }
            GEO::index_t nb_moebius = nb_remove_t - nb_non_manifold;
            GEO::Logger::out("Co3Ne")
                    << "Removing " << nb_remove_t
                    << " triangles ("
                    << nb_non_manifold << " non_manifold, "
                    << nb_moebius
                    << " moebius)"
                    << std::endl;
            M_.facets.delete_elements(remove_t, false);
        }

        // We need to re-compute next_c_around_v_ and v2c_
        // since all the indices changed in the mesh
        // (even if remove_t is empty, because mesh_reorient() may
        // have changed triangles orientation).
        next_c_around_v_.assign(M_.facet_corners.nb(), GEO::index_t(NO_CORNER));
        v2c_.assign(M_.vertices.nb(), GEO::index_t(NO_CORNER));
        for (GEO::index_t t = 0; t < M_.facets.nb(); ++t) {
            insert(t);
        }
    }

    bool Co3NeManifoldExtraction::connect_and_validate_triangle(GEO::index_t t, bool &classified){
        GEO::index_t adj_c[3];
        classified = false;

        //   Combinatorial test (I): tests whether the three
        // candidate edges are manifold.
        if (!get_adjacent_corners(t, adj_c)) {
            classified = true;
            return false;
        }

        //   Geometric test: tests whether the angles formed with
        // the candidate neighbors do not indicate degenerate sharp
        // creases.
        for (GEO::index_t i = 0; i < 3; ++i) {
            if (adj_c[i] != NO_CORNER) {
                GEO::index_t t2 = c2f(adj_c[i]);
                if (!triangles_normals_agree(t, t2)) {
                    classified = true;
                    return false;
                }
            }
        }

        int nb_neighbors =
                (adj_c[0] != NO_CORNER) +
                (adj_c[1] != NO_CORNER) +
                (adj_c[2] != NO_CORNER);

        // Combinatorial test (II)
        switch (nb_neighbors) {
            // If the candidate triangle is adjacent to no other
            // triangle, reject it
            case 0: {
                return false;
            }
                // If the candidate triangle is adjacent to a single
                // triangle, reject it if the vertex opposite to
                // the common edge is not isolated.
            case 1: {
                // If not in strict mode, we reject the triangle.
                // Experimentally, it improves the result.
                if (!strict_) {
                    return false;
                }
                GEO::index_t other_vertex = GEO::index_t(-1);
                for (GEO::index_t i = 0; i < 3; ++i) {
                    if (adj_c[i] != NO_CORNER) {
                        other_vertex =
                                M_.facet_corners.vertex(
                                        M_.facets.corners_begin(t) + ((i + 2) % 3)
                                );
                    }
                }
                geo_debug_assert(other_vertex != GEO::index_t(-1));
                // Test whether other_vertex is isolated, reject
                // the triangle if other_vertex is NOT isolated.
                GEO::index_t nb_incident_T = nb_incident_triangles(other_vertex);
                geo_assert(nb_incident_T != 0); // There is at least THIS T.
                if (nb_incident_T > 1) {
                    return false;
                }
            }
        }

        connect_adjacent_corners(t, adj_c);

        // Combinatorial test (III): test non-manifold vertices
        for (
                GEO::index_t c = M_.facets.corners_begin(t);
                c < M_.facets.corners_end(t); ++c
                ) {
            GEO::index_t v = M_.facet_corners.vertex(c);
            bool moebius = false;
            if (vertex_is_non_manifold_by_excess(v, moebius)) {
                classified = true;
                return false;
            }
            // It should not occur since we remove all Moebius configs
            // from the T3s and forbid Moebius configs when inserting
            // the T12s. However, some transient moebius configurations
            // due to triangle t may appear (since the Moebius test is
            // right after the non-manifold test).
            if (moebius) {
                GEO::Logger::warn("Co3Ne")
                        << "Encountered Moebius configuration" << std::endl;
                classified = true;
                return false;
            }
        }

        // Combinatorial test (IV): orientability
        if (!enforce_orientation_from_triangle(t)) {
            return false;
        }

        classified = true;
        return true;
    }

    bool Co3NeManifoldExtraction::enforce_orientation_from_triangle(GEO::index_t t) {

        // Index of adjacent triangle
        // (or NO_FACET if no neighbor)
        GEO::index_t adj[3];

        // Index of adjacent connected component
        // (or NO_CNX if no neighbor)
        GEO::index_t adj_cnx[3];

        //   Orientation of adjacent triangle relative to
        // triangle t (or 0 if no neighbor)
        GEO::signed_index_t adj_ori[3];

        for (GEO::index_t i = 0; i < 3; ++i) {
            GEO::index_t c = M_.facets.corners_begin(t) + i;
            adj[i] = GEO::index_t(M_.facet_corners.adjacent_facet(c));
        }


        for (GEO::index_t i = 0; i < 3; ++i) {
            if (adj[i] == NO_FACET) {
                adj_ori[i] = 0;
                adj_cnx[i] = NO_CNX;
            } else {
                adj_ori[i] =
                        (triangles_have_same_orientation(t, adj[i])) ? 1 : -1;
                adj_cnx[i] = cnx_[adj[i]];
            }
        }

        //  If in the neighborhood the same connected component appears
        // with two opposite orientations, then connecting the triangle
        // would create a Moebius strip (the triangle is rejected)
        for (GEO::index_t i = 0; i < 3; ++i) {
            if (adj[i] != NO_FACET) {
                for (GEO::index_t j = i + 1; j < 3; ++j) {
                    if (
                            adj_cnx[j] == adj_cnx[i] &&
                            adj_ori[j] != adj_ori[i]
                            ) {
                        return false;
                    }
                }
            }
        }

        //  The triangle is accepted,
        // now reorient all the connected components and the
        // triangle coherently.

        // Find the largest component incident to t
        GEO::index_t largest_neigh_comp = NO_CNX;
        for (GEO::index_t i = 0; i < 3; ++i) {
            if (
                    adj_cnx[i] != NO_CNX && (
                            largest_neigh_comp == NO_CNX ||
                            cnx_size_[adj_cnx[i]] >
                            cnx_size_[adj_cnx[largest_neigh_comp]]
                    )
                    ) {

                largest_neigh_comp = i;
            }
        }
        geo_assert(largest_neigh_comp != NO_CNX);

        // Orient t like the largest incident component
        GEO::index_t comp = adj_cnx[largest_neigh_comp];

        cnx_.resize(std::max(t + 1, cnx_.size()));
        cnx_[t] = comp;
        ++cnx_size_[comp];
        if (adj_ori[largest_neigh_comp] == -1) {
            flip_triangle(t);
            for (GEO::index_t i = 0; i < 3; ++i) {
                adj_ori[i] = -adj_ori[i];
            }
        }

        // Merge (and reorient if need be) all the other incident
        // components
        for (GEO::index_t i = 0; i < 3; ++i) {
            if (
                    i != largest_neigh_comp &&
                    adj[i] != NO_FACET && cnx_[adj[i]] != comp
                    ) {
                merge_connected_component(
                        adj[i], comp, (adj_ori[i] == -1)
                );
            }
        }

        return true;
    }

    GEO::index_t Co3NeManifoldExtraction::add_triangle(GEO::index_t i, GEO::index_t j, GEO::index_t k) {
        GEO::index_t result = M_.facets.create_triangle(i, j, k);
        next_c_around_v_.push_back(GEO::index_t(NO_CORNER));
        next_c_around_v_.push_back(GEO::index_t(NO_CORNER));
        next_c_around_v_.push_back(GEO::index_t(NO_CORNER));
        insert(result);
        return result;
    }

    void Co3NeManifoldExtraction::flip_triangle(GEO::index_t t){

        // Remove t from the additional combinatorial data structure
        // (it is both simpler and more efficient to do that
        //  than updating it).
        remove(
                t,
                false // disconnect is set to false because
                // we will re-insert t right after.
        );

        GEO::index_t c1 = M_.facets.corners_begin(t);
        GEO::index_t c2 = c1 + 1;
        GEO::index_t c3 = c2 + 1;
        GEO::index_t v1 = M_.facet_corners.vertex(c1);
        GEO::index_t f1 = M_.facet_corners.adjacent_facet(c1);
        GEO::index_t f2 = M_.facet_corners.adjacent_facet(c2);
        GEO::index_t v3 = M_.facet_corners.vertex(c3);

        M_.facet_corners.set_vertex(c1, v3);
        M_.facet_corners.set_adjacent_facet(c1, f2);
        M_.facet_corners.set_adjacent_facet(c2, f1);
        M_.facet_corners.set_vertex(c3, v1);

        // Re-insert t into the additional combinatorial data structure.
        insert(t);
    }

    void Co3NeManifoldExtraction::insert(GEO::index_t t) {
        for (
                GEO::index_t c = M_.facets.corners_begin(t);
                c < M_.facets.corners_end(t); ++c
                ) {
            GEO::index_t v = M_.facet_corners.vertex(c);
            if (v2c_[v] == NO_CORNER) {
                v2c_[v] = c;
                next_c_around_v_[c] = c;
            } else {
                next_c_around_v_[c] = next_c_around_v_[v2c_[v]];
                next_c_around_v_[v2c_[v]] = c;
            }
        }
    }

    void Co3NeManifoldExtraction::remove(GEO::index_t t, bool disconnect) {
        if (disconnect) {
            for (
                    GEO::index_t c = M_.facets.corners_begin(t);
                    c < M_.facets.corners_end(t); ++c
                    ) {

                // Disconnect facet-facet link that point to t
                GEO::index_t t2 = M_.facet_corners.adjacent_facet(c);
                if (t2 != NO_FACET) {
                    for (
                            GEO::index_t c2 = M_.facets.corners_begin(GEO::index_t(t2));
                            c2 < M_.facets.corners_end(GEO::index_t(t2));
                            ++c2
                            ) {
                        if (
                                M_.facet_corners.adjacent_facet(c2) == t
                                ) {
                            M_.facet_corners.set_adjacent_facet(
                                    c2, NO_FACET
                            );
                        }
                    }
                }
            }
        }


        for (
                GEO::index_t c = M_.facets.corners_begin(t);
                c < M_.facets.corners_end(t); ++c
                ) {
            // Remove t from combinatorial data structures
            GEO::index_t v = M_.facet_corners.vertex(c);
            if (next_c_around_v_[c] == c) {
                v2c_[v] = NO_CORNER;
            } else {
                GEO::index_t c_pred = next_c_around_v_[c];
                while (next_c_around_v_[c_pred] != c) {
                    c_pred = next_c_around_v_[c_pred];
                }
                next_c_around_v_[c_pred] = next_c_around_v_[c];
                v2c_[v] = c_pred;
            }
        }
    }

    bool Co3NeManifoldExtraction::vertex_is_non_manifold_by_excess(GEO::index_t v, bool &moebius){
        GEO::index_t nb_v_neighbors = nb_incident_triangles(v);
        GEO::index_t c = v2c_[v];
        do {
            GEO::index_t loop_size = 0;
            GEO::index_t c_cur = c;
            do {
                ++loop_size;
                if (c_cur == NO_CORNER) {
                    break;
                }
                if (loop_size > 100) {
                    // Probably Moebious strip or something...
                    moebius = true;
                    break;
                }
                c_cur = next_around_vertex_unoriented(v, c_cur);
            } while (c_cur != c);

            if (c_cur == c && loop_size < nb_v_neighbors) {
                return true;
            }
            c = next_c_around_v_[c];
        } while (c != v2c_[v]);

        return false;
    }

    GEO::index_t Co3NeManifoldExtraction::next_around_vertex_unoriented(GEO::index_t v, GEO::index_t c1) const  {
        GEO::index_t f1 = c2f(c1);
        GEO::index_t v1 = M_.facet_corners.vertex(c1);
        GEO::index_t v2 = M_.facet_corners.vertex(
                M_.facets.next_corner_around_facet(f1, c1)
        );

        geo_debug_assert(v1 == v || v2 == v);

        GEO::index_t f2 = M_.facet_corners.adjacent_facet(c1);
        if (f2 != NO_FACET) {
            for (
                    GEO::index_t c2 = M_.facets.corners_begin(f2);
                    c2 < M_.facets.corners_end(f2);
                    ++c2
                    ) {
                GEO::index_t w1 = M_.facet_corners.vertex(c2);
                GEO::index_t w2 = M_.facet_corners.vertex(
                        M_.facets.next_corner_around_facet(f2, c2)
                );
                if (
                        (v1 == w1 && v2 == w2) ||
                        (v1 == w2 && v2 == w1)
                        ) {
                    if (w2 == v) {
                        return M_.facets.next_corner_around_facet(f2, c2);
                    } else {
                        geo_debug_assert(w1 == v);
                        return M_.facets.prev_corner_around_facet(f2, c2);
                    }
                }
            }
        }
        return NO_CORNER;
    }

    bool Co3NeManifoldExtraction::get_adjacent_corners(GEO::index_t t1, GEO::index_t *adj_c) {
        for (
                GEO::index_t c1 = M_.facets.corners_begin(t1);
                c1 < M_.facets.corners_end(t1); ++c1
                ) {
            GEO::index_t v2 = M_.facet_corners.vertex(
                    M_.facets.next_corner_around_facet(t1, c1)
            );

            *adj_c = NO_CORNER;

            // Traverse the circular incident edge list
            GEO::index_t c2 = next_c_around_v_[c1];
            while (c2 != c1) {
                GEO::index_t t2 = c2f(c2);
                GEO::index_t c3 = M_.facets.prev_corner_around_facet(t2, c2);
                GEO::index_t v3 = M_.facet_corners.vertex(c3);
                if (v3 == v2) {
                    // Found an adjacent edge
                    if (*adj_c == NO_CORNER) {
                        *adj_c = c3;
                        geo_debug_assert(c3 != c1);
                    } else {
                        // If there was already an adjacent edge,
                        // then this is a non-manifold configuration
                        return false;
                    }
                }

                // Check with the other (wrong) orientation
                c3 = M_.facets.next_corner_around_facet(t2, c2);
                v3 = M_.facet_corners.vertex(c3);
                if (v3 == v2) {
                    // Found an adjacent edge
                    if (*adj_c == NO_CORNER) {
                        *adj_c = c2;
                        geo_debug_assert(c2 != c1);
                    } else {
                        // If there was already an adjacent edge,
                        // then this is a non-manifold configuration
                        return false;
                    }
                }
                c2 = next_c_around_v_[c2];
            }
            ++adj_c;
        }
        return true;
    }

    bool Co3NeManifoldExtraction::triangles_have_same_orientation(GEO::index_t t1, GEO::index_t t2) {
        GEO::index_t c1 = M_.facets.corners_begin(t1);
        GEO::index_t i1 = M_.facet_corners.vertex(c1);
        GEO::index_t j1 = M_.facet_corners.vertex(c1 + 1);
        GEO::index_t k1 = M_.facet_corners.vertex(c1 + 2);

        GEO::index_t c2 = M_.facets.corners_begin(t2);
        GEO::index_t i2 = M_.facet_corners.vertex(c2);
        GEO::index_t j2 = M_.facet_corners.vertex(c2 + 1);
        GEO::index_t k2 = M_.facet_corners.vertex(c2 + 2);

        if (
                (i1 == i2 && j1 == j2) ||
                (i1 == k2 && j1 == i2) ||
                (i1 == j2 && j1 == k2) ||
                (k1 == k2 && i1 == i2) ||
                (k1 == j2 && i1 == k2) ||
                (k1 == i2 && i1 == j2) ||
                (j1 == j2 && k1 == k2) ||
                (j1 == i2 && k1 == j2) ||
                (j1 == k2 && k1 == i2)
                ) {
            return false;
        }

        return true;
    }

    bool Co3NeManifoldExtraction::triangles_normals_agree(GEO::index_t t1, GEO::index_t t2) const  {
        const GEO::vec3 *points =
                reinterpret_cast<const GEO::vec3 *>(M_.vertices.point_ptr(0));

        GEO::index_t c1 = M_.facets.corners_begin(t1);
        GEO::index_t i1 = M_.facet_corners.vertex(c1);
        GEO::index_t j1 = M_.facet_corners.vertex(c1 + 1);
        GEO::index_t k1 = M_.facet_corners.vertex(c1 + 2);

        GEO::index_t c2 = M_.facets.corners_begin(t2);
        GEO::index_t i2 = M_.facet_corners.vertex(c2);
        GEO::index_t j2 = M_.facet_corners.vertex(c2 + 1);
        GEO::index_t k2 = M_.facet_corners.vertex(c2 + 2);

        GEO::vec3 n1 = normalize(
                cross(
                        points[j1] - points[i1],
                        points[k1] - points[i1]
                )
        );

        GEO::vec3 n2 = normalize(
                cross(
                        points[j2] - points[i2],
                        points[k2] - points[i2]
                )
        );

        double d = dot(n1, n2);
        // Test for combinatorial orientation,
        // if t1 and t2 have opposite orientation,
        // then we flip one of the normals (i.e.,
        // we simply change the sign of the dot product).
        if (
                (i1 == i2 && j1 == j2) ||
                (i1 == k2 && j1 == i2) ||
                (i1 == j2 && j1 == k2) ||
                (k1 == k2 && i1 == i2) ||
                (k1 == j2 && i1 == k2) ||
                (k1 == i2 && i1 == j2) ||
                (j1 == j2 && k1 == k2) ||
                (j1 == i2 && k1 == j2) ||
                (j1 == k2 && k1 == i2)
                ) {
            d = -d;
        }
        return (d > -0.8);
    }

    void Co3NeManifoldExtraction::merge_connected_component(GEO::index_t t, GEO::index_t comp2, bool flip){
        geo_assert(comp2 != cnx_[t]);

        std::stack<GEO::index_t> S;
        GEO::index_t comp1 = cnx_[t];


        cnx_[t] = comp2;
        --cnx_size_[comp1];
        ++cnx_size_[comp2];
        if (flip) {
            flip_triangle(t);
        }
        S.push(t);
        while (!S.empty()) {
            GEO::index_t t1 = S.top();
            S.pop();
            for (
                    GEO::index_t c = M_.facets.corners_begin(t1);
                    c < M_.facets.corners_end(t1); ++c
                    ) {
                GEO::index_t t2 = M_.facet_corners.adjacent_facet(c);
                if (t2 != NO_FACET && cnx_[t2] == comp1) {
                    cnx_[t2] = comp2;
                    --cnx_size_[comp1];
                    ++cnx_size_[comp2];
                    if (flip) {
                        flip_triangle(t2);
                    }
                    S.push(t2);
                }
            }
        }
        geo_assert(cnx_size_[comp1] == 0);
    }

    void Co3NeManifoldExtraction::init_connected_components(){
        cnx_.assign(M_.facets.nb(), GEO::index_t(NO_CNX));
        cnx_size_.clear();
        for (GEO::index_t t = 0; t < M_.facets.nb(); ++t) {
            if (cnx_[t] == NO_CNX) {
                GEO::index_t cnx_id = cnx_size_.size();
                GEO::index_t nb = 0;
                std::stack<GEO::index_t> S;
                S.push(t);
                cnx_[t] = cnx_id;
                ++nb;
                while (!S.empty()) {
                    GEO::index_t t2 = S.top();
                    S.pop();
                    for (
                            GEO::index_t c = M_.facets.corners_begin(t2);
                            c < M_.facets.corners_end(t2); ++c
                            ) {
                        GEO::index_t t3 = M_.facet_corners.adjacent_facet(c);
                        if (t3 != NO_FACET && cnx_[t3] != cnx_id) {
                            geo_assert(cnx_[t3] == NO_CNX);
                            cnx_[t3] = cnx_id;
                            ++nb;
                            S.push(t3);
                        }
                    }
                }
                cnx_size_.push_back(nb);
            }
        }
        GEO::Logger::out("Co3Ne")
                << "Found " << cnx_size_.size() << " connected components"
                << std::endl;
    }
}