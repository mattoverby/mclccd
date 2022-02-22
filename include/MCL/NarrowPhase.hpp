// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_NARROWPHASE_HPP
#define MCL_NARROWPHASE_HPP 1

#include <Eigen/Core>

namespace mcl
{

template <typename T, int DIM>
class NarrowPhase
{
public:
    typedef Eigen::Matrix<T,DIM,1> VecType;
    using Vec3 = Eigen::Matrix<T,3,1>;
    using Vec2 = Eigen::Matrix<T,2,1>;

    static bool hit_wrong_side_vf(const VecType *v0, const VecType *v1, T t);

    // Used in vertex-facet queries,
    // tests vertex as a ray against the box.
    static bool query_ray_box(
        const VecType &p_x0, const VecType &p_x1,
        const VecType &bmin, const VecType &bmax);

    // Discrete tests for face-face and edge-edge
    // only defined for DIM=3 and DIM=2 respectively

    // f0 and f1 are arrays of three vertices
    static bool discrete_tri_tri(
        const Vec3 *f0, const Vec3 *f1);

    // e0 and e1 are arrays of two vertices
    static bool discrete_edge_edge(
        const Vec2 *e0, const Vec2 *e1);
};

// Even though float is allowed as a template, everything is casted
// to double before evaluating. Kernels for floats are future work.
template <typename T, int DIM>
class NarrowPhaseCTCD
{
public:
    typedef Eigen::Matrix<T,DIM,1> VecType;
    using Vec3 = Eigen::Matrix<T,3,1>;
    using Vec2 = Eigen::Matrix<T,2,1>;

    // Returns:
    // 0 = not colliding
    // 1 = is colliding
    // -1 = colliding wrong side
    // Can only return -1 if test_wrong_side=true

    static int query_ccd_vf(
        const VecType *v0, const VecType *v1,
        const T &eta, // gap
        bool test_wrong_side,
        T &t_impact);

    static int query_ccd_ee(
        const VecType *v0, const VecType *v1,
        const T &eta, // gap
        bool test_vv_and_ve,
        T &t_impact);
};

// ACCD from https://doi.org/10.1145/3450626.3459767
// which is a variant of conservative advancement.
template <typename T, int DIM>
class NarrowPhaseACCD
{
public:
    typedef Eigen::Matrix<T,DIM,1> VecType;

    // Returns:
    // 0 = not colliding
    // 1 = is colliding
    // -1 = colliding wrong side
    // Can only return -1 if test_wrong_side=true

    static int query_ccd_vf(
        const VecType *v0, const VecType *v1,
        const T &eta, // gap
        bool test_wrong_side,
        T &t_impact);

    static int query_ccd_ee(
        const VecType *v0, const VecType *v1,
        const T &eta, // gap
        T &t_impact);
};

} // ns mcl

#endif
