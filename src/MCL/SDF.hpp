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

    SDF();                                  ///< default constructor
    ~SDF();                                 ///< default destructor
    SDF<T>& operator=(SDF<T> const& other); ///< recalculates SDF
    SDF(const SDF<T>& other);               ///< recalculates SDF

    /// @brief Creates the signed distance field.
    /// @param vertices nv x 3 vertices
    /// @param triangles nf x 3 triangles
    void create(const Eigen::Matrix<T, Eigen::Dynamic, 3, Eigen::RowMajor>& vertices,
                const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& triangles);

    /// @brief Project to triangle surface: returns {nearest_face, barycoords, distance}
    /// nearest_face = -1 if there was an error.
    /// @param x vertex to project
    std::tuple<int, Eigen::Vector3<T>, T> project_to_surface(const Eigen::Vector3<T>& x);
};

} // end namespace ccd
} // end namespace mcl

#endif // MCL_CCD_SDF_HPP
