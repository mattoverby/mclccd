// Copyright Matt Overby 2026.
// Distributed under the MIT License.

#ifndef MCL_CCD_SDF_HPP
#define MCL_CCD_SDF_HPP 1

#include <Eigen/Geometry>

#include <memory>

namespace mcl {
namespace ccd {

/// @brief Placeholder, only does surface projection for now
template<typename T>
class SDF
{
public:
    struct SDFData; // forward declare
    std::unique_ptr<SDFData> data;

    /// @brief Creates an empty signed distance field.
    SDF();

    /// @brief Creates the signed distance field.
    /// @param vertices nv x dim vertices
    /// @param triangles nt x 3 triangles
    template<typename DerivedV, typename DerivedP>
    SDF(const Eigen::MatrixBase<DerivedV>& vertices, const Eigen::MatrixBase<DerivedP>& triangles)
    {
        Eigen::Matrix<T,Eigen::Dynamic,3,Eigen::RowMajor> V(vertices.template cast<T>());
        Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> P(triangles.template cast<int>());
        create(V.data(), V.rows(), P.data(), P.rows());
    }

    /// @brief Creates the signed distance field.
    /// @param vertex_data nv x 3 vertex data
    /// @param num_vertices nv
    /// @param triangle_data nt x 3 triangle index data
    /// @param num_triangles nt
    void create(const T* vertex_data, int num_vertices, const int* triangle_data, int num_triangles);

    /// @brief Project to (rest) face: returns {nearest_face, barycoords, distance}
    /// nearest_face = -1 if there was an error.
    /// @param x vertex to project
    std::tuple<int, Eigen::Vector3<T>, T> project_to_surface(const Eigen::Vector3<T> &x);

};

} // end namespace ccd
} // end namespace mcl

#endif // MCL_CCD_SDF_HPP
