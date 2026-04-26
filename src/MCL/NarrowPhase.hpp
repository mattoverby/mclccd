// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_CCD_NARROWPHASE_HPP
#define MCL_CCD_NARROWPHASE_HPP 1

#include <Eigen/Core>

namespace mcl {
namespace ccd {

template<typename T>
class NarrowPhase
{
  public:
    using Vec3 = Eigen::Vector3<T>;
    using Vec2 = Eigen::Vector2<T>;

    /// @brief Returns true if VF contact is on opposite side of normal.
    static bool hit_wrong_side_vf(const Vec2* v0, const Vec2* v1, T t);

    /// @brief Returns true if VF contact is on opposite side of normal.
    static bool hit_wrong_side_vf(const Vec3* v0, const Vec3* v1, T t);

    // Used in vertex-facet queries,
    // tests vertex as a ray against the box.
    static bool query_edge_box(const Vec3& p_x0, const Vec3& p_x1, const Vec3& bmin, const Vec3& bmax);

    /// @brief Returns true if two triangles intersect.
    static bool discrete_tri_tri(const Vec3& p0,
                                 const Vec3& p1,
                                 const Vec3& p2,
                                 const Vec3& q0,
                                 const Vec3& q1,
                                 const Vec3& q2);

    /// @brief Returns true if two edges are intersecting
    static bool discrete_edge_edge(const Vec2& p0, const Vec2& p1, const Vec2& q0, const Vec2& q1);

    /// @brief Returns true if the point is inside the tetrahedron.
    static bool point_in_tet(const Vec3& p, const Vec3& v0, const Vec3& v1, const Vec3& v2, const Vec3& v3);
};

// Even though float is allowed as a template, everything is casted
// to double before evaluating. Kernels for floats are future work.
template<typename T, int DIM>
class NarrowPhaseCTCD
{
  public:
    // Returns:
    // 0 = not colliding
    // 1 = is colliding
    // -1 = colliding wrong side
    // Can only return -1 if test_wrong_side=true

    static int query_ccd_vf(const Eigen::Vector<T, DIM>* v0,
                            const Eigen::Vector<T, DIM>* v1,
                            const T& eta, // gap
                            bool test_wrong_side,
                            T& t_impact);

    static int query_ccd_ee(const Eigen::Vector<T, DIM>* v0,
                            const Eigen::Vector<T, DIM>* v1,
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
    // Returns:
    // 0 = not colliding
    // 1 = is colliding

    static int query_ccd_vf(const Eigen::Vector<T, DIM>* v0,
                            const Eigen::Vector<T, DIM>* v1,
                            const T& eta, // gap
                            T& t_impact);

    static int query_ccd_ee(const Eigen::Vector<T, DIM>* v0,
                            const Eigen::Vector<T, DIM>* v1,
                            const T& eta, // gap
                            T& t_impact);

    // Alg 1. from doi: 10.1145/3450626.3459767
    static bool additive_ccd(const Eigen::Vector<T, DIM>* v0,
                             const Eigen::Vector<T, DIM>* v1,
                             const T& eta, // gap
                             bool is_vf,
                             T& t_impact);

    // Returns distance between pairs (VF or EE)
    // If error, returns -1
    static T pair_distance(const Eigen::Vector<T, DIM>* v, bool is_vf);
};

} // end namespace ccd
} // end namespace mcl

#endif
