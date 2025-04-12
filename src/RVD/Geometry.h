//
// Created by hongbo on 10/04/25.
//

#ifndef SIMPLEXCVT_RECON_GEOMETRY_H
#define SIMPLEXCVT_RECON_GEOMETRY_H
#include <geogram/basic/geometry.h>
#include <Eigen/Dense>
namespace GEO_BASE{
    class Vertex {
    public:
        /**
         * \brief Constructs a new uninitialized Vertex.
         */
        Vertex() {
        }

        /**
         * \brief Constructs a Vertex from a 3d point.
         */
        Vertex(const GEO::vec3& v) :
                point_(v),
                adjacent_seed_(-1) {
        }

        /**
         * \brief Gets the 3d point associated with this vertex.
         * \return a const reference to the 3d point
         */
        const GEO::vec3& point() const {
            return point_;
        }

        /**
         * \brief Gets the 3d point associated with this vertex.
         * \return a const reference to the 3d point
         */
        GEO::vec3& point() {
            return point_;
        }

        /**
         * \brief Gets the index of the adjacent seed associated with
         *  this vertex.
         * \details Each vertex stores combinatorial information, i.e.
         *  the index of the adjacent Voronoi seed accros the edge
         *  starting from this vertex
         * \return the index of the adjacent Voronoi seed
         */
        GEO::signed_index_t adjacent_seed() const {
            return adjacent_seed_;
        }

        /**
         * \brief Sets the index of the adjacent seed associated with
         *  this vertex.
         * \details Each vertex stores combinatorial information, i.e.
         *  the index of the adjacent Voronoi seed accros the edge
         *  starting from this vertex
         * \param[in] x the index of the adjacent Voronoi seed
         */
        void set_adjacent_seed(GEO::signed_index_t x) {
            adjacent_seed_ = x;
        }

    private:
        GEO::vec3 point_;
        GEO::signed_index_t adjacent_seed_;
    };

    class Vertex_hd {
    public:
        /**
         * \brief Constructs a new uninitialized Vertex.
         */
        Vertex_hd() {
        }

        /**
         * \brief Constructs a Vertex from a hd point.
         */
        Vertex_hd(Eigen::VectorXd v) :
                point_(std::move(v)),
                adjacent_seed_(-1) {
        }

        /**
         * \brief Gets the 3d point associated with this vertex.
         * \return a const reference to the 3d point
         */
        const Eigen::VectorXd& point() const {
            return point_;
        }

        /**
         * \brief Gets the 3d point associated with this vertex.
         * \return a const reference to the 3d point
         */
        Eigen::VectorXd& point() {
            return point_;
        }

        /**
         * \brief Gets the index of the adjacent seed associated with
         *  this vertex.
         * \details Each vertex stores combinatorial information, i.e.
         *  the index of the adjacent Voronoi seed accros the edge
         *  starting from this vertex
         * \return the index of the adjacent Voronoi seed
         */
        GEO::signed_index_t adjacent_seed() const {
            return adjacent_seed_;
        }

        /**
         * \brief Sets the index of the adjacent seed associated with
         *  this vertex.
         * \details Each vertex stores combinatorial information, i.e.
         *  the index of the adjacent Voronoi seed accros the edge
         *  starting from this vertex
         * \param[in] x the index of the adjacent Voronoi seed
         */
        void set_adjacent_seed(GEO::signed_index_t x) {
            adjacent_seed_ = x;
        }

    private:
        Eigen::VectorXd point_;
        GEO::signed_index_t adjacent_seed_;
    };

    /**
     * \brief Internal representation of the polygons, that represent
     *  the intersection between the disks and the Voronoi cells.
     */
    class Polygon {
    public:
        /**
         * \brief Creates a new uninitialized polygon with a given
         *  number of vertices.
         * \param[in] size number of vertices
         */
        Polygon(GEO::index_t size) :
                vertices_(size) {
        }

        /**
         * \brief Gets the number of vertices.
         * \return the number of vertices of this Polygon
         */
        GEO::index_t nb_vertices() const {
            return vertices_.size();
        }

        /**
         * \brief Adds a new vertex to this Polygon.
         * \param[in] v the vertex to be added.
         */
        void add_vertex(const Vertex& v) {
            vertices_.push_back(v);
        }

        /**
         * \brief Gets a Vertex by its index.
         * \param[in] i the index of the Vertex
         * \return a reference to the Vertex
         */
        Vertex& vertex(GEO::index_t i) {
            return vertices_[i];
        }

        /**
         * \brief Gets a Vertex by its index.
         * \param[in] i the index of the Vertex
         * \return a const reference to the Vertex
         */
        const Vertex& vertex(GEO::index_t i) const {
            return vertices_[i];
        }

        /**
         * \brief Gets the index of the next vertex around
         *  the polygon.
         * \param[in] i index of the vertex
         * \return index of the next vertex (successor of \p i)
         *  around the Polygon.
         */
        GEO::index_t next_vertex(GEO::index_t i) const {
            return (i == nb_vertices() - 1) ? 0 : i + 1;
        }

        /**
         * \brief Removes all the vertices.
         */
        void clear() {
            vertices_.resize(0);
        }

        /**
         * \brief Swaps the vertices of this Polygon with
         *  the vertices of another polygon.
         * \param[in] P the other polygon
         */
        void swap(Polygon& P) {
            vertices_.swap(P.vertices_);
        }

    private:
        GEO::vector<Vertex> vertices_;
    };

    class Polygon_hd {
    public:
        /**
         * \brief Creates a new uninitialized polygon with a given
         *  number of vertices.
         * \param[in] size number of vertices
         */
        Polygon_hd(GEO::index_t size) :
                vertices_(size) {
        }

        /**
         * \brief Gets the number of vertices.
         * \return the number of vertices of this Polygon
         */
        GEO::index_t nb_vertices() const {
            return vertices_.size();
        }

        /**
         * \brief Adds a new vertex to this Polygon.
         * \param[in] v the vertex to be added.
         */
        void add_vertex(const Vertex_hd& v) {
            vertices_.push_back(v);
        }

        /**
         * \brief Gets a Vertex by its index.
         * \param[in] i the index of the Vertex
         * \return a reference to the Vertex
         */
        Vertex_hd& vertex(GEO::index_t i) {
            return vertices_[i];
        }

        /**
         * \brief Gets a Vertex by its index.
         * \param[in] i the index of the Vertex
         * \return a const reference to the Vertex
         */
        const Vertex_hd& vertex(GEO::index_t i) const {
            return vertices_[i];
        }

        /**
         * \brief Gets the index of the next vertex around
         *  the polygon.
         * \param[in] i index of the vertex
         * \return index of the next vertex (successor of \p i)
         *  around the Polygon.
         */
        GEO::index_t next_vertex(GEO::index_t i) const {
            return (i == nb_vertices() - 1) ? 0 : i + 1;
        }

        /**
         * \brief Removes all the vertices.
         */
        void clear() {
            vertices_.resize(0);
        }

        /**
         * \brief Swaps the vertices of this Polygon with
         *  the vertices of another polygon.
         * \param[in] P the other polygon
         */
        void swap(Polygon_hd& P) {
            vertices_.swap(P.vertices_);
        }

    private:
        GEO::vector<Vertex_hd> vertices_;
    };

}

#endif //SIMPLEXCVT_RECON_GEOMETRY_H
