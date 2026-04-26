// Copyright Matt Overby 2026.
// Distributed under the MIT License.

#include "SDF.hpp"

#include "../../third-party/TriangleMeshDistance.hpp"

namespace mcl {
namespace ccd {

template<typename T>
struct SDF<T>::SDFData
{
    Eigen::Matrix<T, Eigen::Dynamic, 3, Eigen::RowMajor> vertices;
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> triangles;
    tmd::TriangleMeshDistance sdf;
};

template<typename T>
SDF<T>::SDF()
{
}

template<typename T>
SDF<T>::~SDF()
{
}

template<typename T>
SDF<T>&
SDF<T>::operator=(SDF<T> const& other)
{
    if (other.data != nullptr) {
        create(other.data->vertices, other.data->triangles);
    }
    return *this;
}

template<typename T>
SDF<T>::SDF(const SDF<T>& other)
{
    *this = other;
}

template<typename T>
void
SDF<T>::create(const Eigen::Matrix<T, Eigen::Dynamic, 3, Eigen::RowMajor>& vertices,
               const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& triangles)
{
    data = std::make_unique<SDFData>();
    data->vertices = vertices;
    data->triangles = triangles;
    data->sdf = tmd::TriangleMeshDistance(vertices.data(), vertices.rows(), triangles.data(), triangles.rows());
}

template<typename T>
std::tuple<int, Eigen::Vector3<T>, T>
SDF<T>::project_to_surface(const Eigen::Vector3<T>& x)
{
    if (!data) {
        return { -1, Eigen::Vector3<T>::Zero(), T(0) };
    }

    auto result = data->sdf.unsigned_distance({ x[0], x[1], x[2] });
    Eigen::Vector3<T> barys(result.barycentric[0], result.barycentric[1], result.barycentric[2]);
    return { result.triangle_id, std::move(barys), T(result.distance) };
}

} // end namespace ccd
} // end namespace mcl

template class mcl::ccd::SDF<double>;
template class mcl::ccd::SDF<float>;
