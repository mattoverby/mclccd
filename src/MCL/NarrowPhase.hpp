// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_CCD_NARROWPHASE_HPP
#define MCL_CCD_NARROWPHASE_HPP 1

#include <Eigen/Core>

namespace mcl {
namespace ccd {

template<typename T, int DIM>
class NarrowPhase
{
  public:
    using VecType = Eigen::Matrix<T, DIM, 1>;
    using Vec3 = Eigen::Vector3<T>;
    using Vec2 = Eigen::Vector2<T>;

    static bool hit_wrong_side_vf(const VecType* v0, const VecType* v1, T t);

    // Used in vertex-facet queries,
    // tests vertex as a ray against the box.
    static bool query_ray_box(const VecType& p_x0, const VecType& p_x1, const VecType& bmin, const VecType& bmax);

    // Discrete tests for face-face and edge-edge
    // only defined for DIM=3 and DIM=2 respectively

    // f0 and f1 are arrays of three vertices
    static bool discrete_tri_tri(const Vec3* f0, const Vec3* f1);

    // e0 and e1 are arrays of two vertices
    static bool discrete_edge_edge(const Vec2* e0, const Vec2* e1);

    /// @brief Returns true if the point is inside the tetrahedron.
    static bool point_in_tet(const Vec3 &p,
        const Vec3 &v0, const Vec3 &v1, const Vec3 &v2, const Vec3 &v3);
};

// Even though float is allowed as a template, everything is casted
// to double before evaluating. Kernels for floats are future work.
template<typename T, int DIM>
class NarrowPhaseCTCD
{
  public:
    typedef Eigen::Matrix<T, DIM, 1> VecType;

    // Returns:
    // 0 = not colliding
    // 1 = is colliding
    // -1 = colliding wrong side
    // Can only return -1 if test_wrong_side=true

    static int query_ccd_vf(const VecType* v0,
                            const VecType* v1,
                            const T& eta, // gap
                            bool test_wrong_side,
                            T& t_impact);

    static int query_ccd_ee(const VecType* v0,
                            const VecType* v1,
                            const T& eta, // gap
                            bool test_vv_and_ve,
                            T& t_impact);
};

// ACCD from https://doi.org/10.1145/3450626.3459767
// which is a variant of conservative advancement.
template<typename T, int DIM>
class NarrowPhaseACCD
{
  public:
    typedef Eigen::Matrix<T, DIM, 1> VecType;

    // Returns:
    // 0 = not colliding
    // 1 = is colliding

    static int query_ccd_vf(const VecType* v0,
                            const VecType* v1,
                            const T& eta, // gap
                            T& t_impact);

    static int query_ccd_ee(const VecType* v0,
                            const VecType* v1,
                            const T& eta, // gap
                            T& t_impact);

    // Alg 1. from doi: 10.1145/3450626.3459767
    static bool additive_ccd(const VecType* v0,
                             const VecType* v1,
                             const T& eta, // gap
                             bool is_vf,
                             T& t_impact);

    // Returns distance between pairs (VF or EE)
    // If error, returns -1
    static T pair_distance(const VecType* v, bool is_vf);
};

} // end namespace ccd
} // end namespace mcl

#endif
