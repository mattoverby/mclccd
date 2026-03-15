// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "NarrowPhase.hpp"
#include "Projection.hpp"

#include "../../third-party/CTCD.hpp"
#include "../../third-party/Distance.hpp"
#include "../../third-party/tt_isect.hpp"

#include <limits>

// I need to fix tabbing...

namespace mcl {
namespace ccd {

// ---------------------------------------------------------
//	Hit Wrong Side
// ---------------------------------------------------------

template<typename T>
bool
NarrowPhase<T>::hit_wrong_side_vf(const Vec2* verts0, const Vec2* verts1, T t)
{
    // Resting contact, relative velocity=0
    // and so there is no way to tell.
    if (std::abs(t) <= T(0)) {
        return false;
    }

    const Vec2& q0start = verts0[0];
    const Vec2& q1start = verts0[1];
    const Vec2& q2start = verts0[2];
    const Vec2& q0end = verts1[0];
    const Vec2& q1end = verts1[1];
    const Vec2& q2end = verts1[2];

    Vec2 p = q0start * (1.0 - t) + q0end * t;
    Vec2 e0 = q1start * (1.0 - t) + q1end * t;
    Vec2 e1 = q2start * (1.0 - t) + q2end * t;
    Vec2 n = Vec2(e1[1] - e0[1], -(e1[0] - e0[0])).stableNormalized();
    double denom = (e0 - e1).norm();
    if (denom < std::numeric_limits<T>::epsilon()) {
        return false; // can't check
    }

    Vec2 barys = Vec2((e1 - p).norm(), (e0 - p).norm()) / denom;
    Vec2 e0bary = barys[0] * q1start + barys[1] * q2start;
    Vec2 e1bary = barys[0] * q1end + barys[1] * q2end;
    Vec2 apex_vel = q0end - q0start;
    Vec2 edge_vel = e1bary - e0bary;
    Vec2 v = apex_vel - edge_vel;
    T vel_v_dot_e = n.dot(v);
    return vel_v_dot_e > T(0);
}

template<typename T>
bool
NarrowPhase<T>::hit_wrong_side_vf(const Vec3* verts0, const Vec3* verts1, T t)
{
    // Resting contact, relative velocity=0
    // and so there is no way to tell.
    if (std::abs(t) <= T(0)) {
        return false;
    }

    // Cast to Vector3d for CTCD kernel
    std::array<Eigen::Vector3d, 4> xt = { (verts0[0] * (T(1) - t) + verts1[0] * t).template cast<double>(),
                                          (verts0[1] * (T(1) - t) + verts1[1] * t).template cast<double>(),
                                          (verts0[2] * (T(1) - t) + verts1[2] * t).template cast<double>(),
                                          (verts0[3] * (T(1) - t) + verts1[3] * t).template cast<double>() };

    Eigen::Vector3d n = (xt[2] - xt[1]).cross(xt[3] - xt[1]);
    if (n.squaredNorm() <= 0) {
        return false;
    }
    n.stableNormalize();

    Eigen::Vector3d barys = Eigen::Vector3d::Zero();
    mcl::ctcd::vertexFaceDistance(xt[0], xt[1], xt[2], xt[3], barys[0], barys[1], barys[2]);

    Vec3 apex_vel = (verts1[0] - verts0[0]);
    Vec3 face_pt0 = (barys[0] * verts0[1] + barys[1] * verts0[2] + barys[2] * verts0[3]);
    Vec3 face_pt1 = (barys[0] * verts1[1] + barys[1] * verts1[2] + barys[2] * verts1[3]);
    Vec3 face_vel = face_pt1 - face_pt0;
    Vec3 v = apex_vel - face_vel;
    double vel_v_dot_f = n.dot(v.template cast<double>());
    bool hit_wrong_side = vel_v_dot_f > 0.0;
    return hit_wrong_side;
}

template<typename T>
bool
NarrowPhase<T>::query_edge_box(const Vec3& p_x0, const Vec3& p_x1, const Vec3& bmin, const Vec3& bmax)
{
    using namespace Eigen;
    Vec3 dir = p_x1 - p_x0;
    const Vec3& origin = p_x0;
    T t0 = 0;
    T t1 = dir.norm();

    // Starts inside box
    AlignedBox<T, 3> box;
    box.extend(bmin);
    box.extend(bmax);
    if (box.contains(p_x0) || box.contains(p_x1)) {
        return true;
    }

    dir.normalize();
    typedef Matrix<T, 1, 3> RowVector3S;
    const RowVector3S inv_dir(1. / dir(0), 1. / dir(1), 1. / dir(2));
    const std::array<bool, 3> sign = { inv_dir(0) < 0, inv_dir(1) < 0, inv_dir(2) < 0 };
    // http://people.csail.mit.edu/amy/papers/box-jgt.pdf
    // "An Efficient and Robust Ray–Box Intersection Algorithm"
    T tymin, tymax, tzmin, tzmax;
    std::array<RowVector3S, 2> bounds = { bmin, bmax };
    T tmin = (bounds[sign[0]](0) - origin(0)) * inv_dir(0);
    T tmax = (bounds[1 - sign[0]](0) - origin(0)) * inv_dir(0);
    tymin = (bounds[sign[1]](1) - origin(1)) * inv_dir(1);
    tymax = (bounds[1 - sign[1]](1) - origin(1)) * inv_dir(1);
    if ((tmin > tymax) || (tymin > tmax)) {
        return false;
    }
    if (tymin > tmin) {
        tmin = tymin;
    }
    if (tymax < tmax) {
        tmax = tymax;
    }
    tzmin = (bounds[sign[2]](2) - origin(2)) * inv_dir(2);
    tzmax = (bounds[1 - sign[2]](2) - origin(2)) * inv_dir(2);
    if ((tmin > tzmax) || (tzmin > tmax)) {
        return false;
    }
    if (tzmin > tmin) {
        tmin = tzmin;
    }
    if (tzmax < tmax) {
        tmax = tzmax;
    }
    if (!((tmin < t1) && (tmax > t0))) {
        return false;
    }
    return true;
}

// ---------------------------------------------------------
//	Discrete tests
// ---------------------------------------------------------

template<typename T>
bool
NarrowPhase<T>::discrete_tri_tri(const Vec3& p0,
                                 const Vec3& p1,
                                 const Vec3& p2,
                                 const Vec3& q0,
                                 const Vec3& q1,
                                 const Vec3& q2)
{
    // Cast to double for precision
    std::array<Eigen::Vector3d, 3> p;
    std::array<Eigen::Vector3d, 3> q;
    p[0] = p0.template cast<double>();
    p[1] = p1.template cast<double>();
    p[2] = p2.template cast<double>();
    q[0] = q0.template cast<double>();
    q[1] = q1.template cast<double>();
    q[2] = q2.template cast<double>();
    return tritri::tri_tri_overlap_test_3d((double*)p[0].data(),
                                           (double*)p[1].data(),
                                           (double*)p[2].data(),
                                           (double*)q[0].data(),
                                           (double*)q[1].data(),
                                           (double*)q[2].data());
}

template<typename T>
bool
NarrowPhase<T>::discrete_edge_edge(const Vec2& p0, const Vec2& p1, const Vec2& q0, const Vec2& q1)
{
    // From https://stackoverflow.com/a/565282
    using namespace Eigen;
    constexpr T eps = std::numeric_limits<T>::epsilon();
    Vec2 n(q0[0] - p0[0], q0[1] - p0[1]);
    Vec2 r(p1[0] - p0[0], p1[1] - p0[1]);
    Vec2 s(q1[0] - q0[0], q1[1] - q0[1]);
    T rxs = r[0] * s[1] - r[1] * s[0];
    if (std::abs(rxs) < eps) {
        return false;
    } // parallel
    T nxr = n[0] * r[1] - n[1] * r[0];
    if (std::abs(nxr) < eps) // collinear
    {
        return ((q0[0] - p0[0] < 0) != (q0[0] - p1[0] < 0)) || ((q0[1] - p0[1] < 0) != (q0[1] - p1[1] < 0));
    }
    T nxs = n[0] * s[1] - n[1] * s[0];
    T rxsr = T(1) / rxs;
    T t = nxs * rxsr;
    T u = nxr * rxsr;
    return (t >= 0) && (t <= 1) && (u >= 0) && (u <= 1);
}

template<typename T>
bool
NarrowPhase<T>::point_in_tet(const Vec3& p, const Vec3& v0, const Vec3& v1, const Vec3& v2, const Vec3& v3)
{
    auto scalar_triple_product = [](const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
        return (b - a).cross(c - a).dot(d - a);
    };

    T s0 = scalar_triple_product(v0, v1, v2, p);
    T s1 = scalar_triple_product(v0, v1, v3, p);
    T s2 = scalar_triple_product(v0, v2, v3, p);
    T s3 = scalar_triple_product(v1, v2, v3, p);
    bool pos = (s0 >= 0 && s1 >= 0 && s2 >= 0 && s3 >= 0);
    bool neg = (s0 <= 0 && s1 <= 0 && s2 <= 0 && s3 <= 0);
    return pos || neg;
}

// ---------------------------------------------------------
//	Query CCD with CTCD
// ---------------------------------------------------------

template<>
int
NarrowPhaseCTCD<double, 2>::query_ccd_vf(const Eigen::Vector2d* verts0,
                                         const Eigen::Vector2d* verts1,
                                         const double& eta,
                                         bool test_wrong_side,
                                         double& t_impact)
{
    // First test is to check AABBs
    Eigen::AlignedBox<double, 2> q_AABB;
    q_AABB.extend(verts0[1]);
    q_AABB.extend(verts1[1]);
    q_AABB.extend(verts0[2]);
    q_AABB.extend(verts1[2]);
    for (int i = 0; i < 2; ++i) {
        q_AABB.min()[i] -= eta;
        q_AABB.max()[i] += eta;
    }

    // TODO
    // if (!NarrowPhase<double>::query_edge_box(verts0[0], verts1[0], q_AABB.min(), q_AABB.max()))
    //    return 0;

    std::vector<double> all_toi;
    bool ve = mcl::ctcd::CTCD::vertexEdgeCTCD(
        verts0[0], verts0[1], verts0[2], verts1[0], verts1[1], verts1[2], eta, t_impact, &all_toi);

    // Wrong side collision?
    if (ve && test_wrong_side) {
        assert(all_toi.size() > 0);

        // Loop t. Find first t that is not wrong-side-collision
        std::sort(all_toi.begin(), all_toi.end());
        bool actually_hit = false;
        int nt = all_toi.size();
        for (int i = 0; i < nt; ++i) {
            bool wrongside = NarrowPhase<double>::hit_wrong_side_vf(verts0, verts1, all_toi[i]);
            if (!wrongside) {
                actually_hit = true;
                t_impact = all_toi[i];
                break;
            }
        }

        if (!actually_hit) {
            return -1;
        } else {
            return 1;
        }

    } // end test wrong side

    return ve;

} // end query ccd vf 2d

template<>
int
NarrowPhaseCTCD<double, 3>::query_ccd_vf(const Eigen::Vector3d* verts0,
                                         const Eigen::Vector3d* verts1,
                                         const double& eta,
                                         bool test_wrong_side,
                                         double& t_impact)
{
    using namespace Eigen;
    typedef AlignedBox<double, 3> AABB;

    // First test is to check AABBs
    {
        AABB q_AABB;
        q_AABB.extend(verts0[1]);
        q_AABB.extend(verts1[1]);
        q_AABB.extend(verts0[2]);
        q_AABB.extend(verts1[2]);
        q_AABB.extend(verts0[3]);
        q_AABB.extend(verts1[3]);
        for (int i = 0; i < 3; ++i) {
            q_AABB.min()[i] -= eta;
            q_AABB.max()[i] += eta;
        }
        if (!NarrowPhase<double>::query_edge_box(verts0[0], verts1[0], q_AABB.min(), q_AABB.max()))
            return 0;
    }

    std::vector<double> all_toi;
    auto ccd_vf = [&]() -> bool {
        if (mcl::ctcd::CTCD::vertexFaceCTCD(verts0[0],
                                            verts0[1],
                                            verts0[2],
                                            verts0[3],
                                            verts1[0],
                                            verts1[1],
                                            verts1[2],
                                            verts1[3],
                                            eta,
                                            t_impact,
                                            &all_toi)) {
            return true;
        }

        // Vertex-face edges
        for (int edge = 0; edge < 3; ++edge) {
            if (mcl::ctcd::CTCD::vertexEdgeCTCD(verts0[0],
                                                verts0[1 + (edge % 3)],
                                                verts0[1 + ((edge + 1) % 3)],
                                                verts1[0],
                                                verts1[1 + (edge % 3)],
                                                verts1[1 + ((edge + 1) % 3)],
                                                eta,
                                                t_impact,
                                                &all_toi)) {
                return true;
            }
        }

        // Vertex-face vertices
        for (int vert = 0; vert < 3; ++vert) {
            if (mcl::ctcd::CTCD::vertexVertexCTCD(
                    verts0[0], verts0[1 + vert], verts1[0], verts1[1 + vert], eta, t_impact, &all_toi)) {
                return true;
            }
        }

        return false;
    };

    int hit = ccd_vf();

    // Wrong side collision?
    if (hit && test_wrong_side) {
        assert((int)all_toi.size() > 0);

        // Loop t. Find first t that is not wrong-side-collision
        std::sort(all_toi.begin(), all_toi.end());
        bool actually_hit = false;
        int nt = all_toi.size();
        for (int i = 0; i < nt; ++i) {
            if (all_toi[i] < 0) {
                continue;
            }
            bool wrongside = NarrowPhase<double>::hit_wrong_side_vf(verts0, verts1, all_toi[i]);
            if (!wrongside) {
                actually_hit = true;
                t_impact = all_toi[i];
                break;
            }
        }

        if (!actually_hit) {
            return -1;
        } else {
            return 1;
        }

    } // end test wrong side

    assert(std::isfinite(t_impact));
    return hit;

} // end query ccd vf

template<>
int
NarrowPhaseCTCD<float, 2>::query_ccd_vf(const Eigen::Vector2f* verts0_,
                                        const Eigen::Vector2f* verts1_,
                                        const float& eta,
                                        bool test_wrong_side,
                                        float& t_impact)
{
    Eigen::Vector2d verts0[3], verts1[3];
    for (int i = 0; i < 3; ++i) {
        verts0[i] = verts0_[i].cast<double>();
        verts1[i] = verts1_[i].cast<double>();
    }
    double t = t_impact;
    int ret = NarrowPhaseCTCD<double, 2>::query_ccd_vf(verts0, verts1, eta, test_wrong_side, t);
    t_impact = t;
    return ret;
}

template<>
int
NarrowPhaseCTCD<float, 3>::query_ccd_vf(const Eigen::Vector3f* verts0_,
                                        const Eigen::Vector3f* verts1_,
                                        const float& eta,
                                        bool test_wrong_side,
                                        float& t_impact)
{
    Eigen::Vector3d verts0[4], verts1[4];
    for (int i = 0; i < 4; ++i) {
        verts0[i] = verts0_[i].cast<double>();
        verts1[i] = verts1_[i].cast<double>();
    }
    double t = t_impact;
    int ret = NarrowPhaseCTCD<double, 3>::query_ccd_vf(verts0, verts1, eta, test_wrong_side, t);
    t_impact = t;
    return ret;
}

template<>
int
NarrowPhaseCTCD<double, 3>::query_ccd_ee(const Eigen::Vector3d* verts0,
                                         const Eigen::Vector3d* verts1,
                                         const double& eta,
                                         bool test_vv_and_ve,
                                         double& t_impact)
{
    using namespace Eigen;
    typedef AlignedBox<double, 3> AABB;

    // First test is to check AABBs
    AABB p_AABB, q_AABB;
    p_AABB.extend(verts0[0]);
    p_AABB.extend(verts1[0]);
    p_AABB.extend(verts0[1]);
    p_AABB.extend(verts1[1]);
    q_AABB.extend(verts0[2]);
    q_AABB.extend(verts1[2]);
    q_AABB.extend(verts0[3]);
    q_AABB.extend(verts1[3]);
    for (int i = 0; i < 3; ++i) {
        p_AABB.min()[i] -= eta;
        p_AABB.max()[i] += eta;
        q_AABB.min()[i] -= eta;
        q_AABB.max()[i] += eta;
    }
    if (!p_AABB.intersects(q_AABB)) {
        return 0;
    }

    std::vector<double> all_toi;
    auto ccd_ee = [&]() -> bool {
        if (mcl::ctcd::CTCD::edgeEdgeCTCD(verts0[0],
                                          verts0[1],
                                          verts0[2],
                                          verts0[3],
                                          verts1[0],
                                          verts1[1],
                                          verts1[2],
                                          verts1[3],
                                          eta,
                                          t_impact,
                                          &all_toi)) {
            return true;
        }

        if (test_vv_and_ve) {
            if (mcl::ctcd::CTCD::vertexEdgeCTCD(
                    verts0[0], verts0[2], verts0[3], verts1[0], verts1[2], verts1[3], eta, t_impact, &all_toi))
                return true;

            if (mcl::ctcd::CTCD::vertexEdgeCTCD(
                    verts0[1], verts0[2], verts0[3], verts1[1], verts1[2], verts1[3], eta, t_impact, &all_toi))
                return true;

            if (mcl::ctcd::CTCD::vertexEdgeCTCD(
                    verts0[2], verts0[0], verts0[1], verts1[2], verts1[0], verts1[1], eta, t_impact, &all_toi))
                return true;

            if (mcl::ctcd::CTCD::vertexEdgeCTCD(
                    verts0[3], verts0[0], verts0[1], verts1[3], verts1[0], verts1[1], eta, t_impact, &all_toi))
                return true;

            if (mcl::ctcd::CTCD::vertexVertexCTCD(verts0[0], verts0[2], verts1[0], verts1[2], eta, t_impact, &all_toi))
                return true;

            if (mcl::ctcd::CTCD::vertexVertexCTCD(verts0[0], verts0[3], verts1[0], verts1[3], eta, t_impact, &all_toi))
                return true;

            if (mcl::ctcd::CTCD::vertexVertexCTCD(verts0[1], verts0[2], verts1[1], verts1[2], eta, t_impact, &all_toi))
                return true;

            if (mcl::ctcd::CTCD::vertexVertexCTCD(verts0[1], verts0[3], verts1[1], verts1[3], eta, t_impact, &all_toi))
                return true;
        }

        return false;
    };

    int hit = ccd_ee();

    // Collision result okay?
    if (hit) {
        assert((int)all_toi.size() > 0);
        std::sort(all_toi.begin(), all_toi.end());

        // Loop t. Find first t that is not parallel
        bool actually_hit = false;
        int nt = all_toi.size();
        for (int i = 0; i < nt; ++i) {
            if (all_toi[i] < 0) {
                continue;
            }

            // Check if the edges are parallel
            // The edgeEdgeCTCD sometimes returns toi for parallel edges!
            std::array<Vector3d, 4> vt = { (1.0 - all_toi[i]) * verts0[0] + all_toi[i] * verts1[0],
                                           (1.0 - all_toi[i]) * verts0[1] + all_toi[i] * verts1[1],
                                           (1.0 - all_toi[i]) * verts0[2] + all_toi[i] * verts1[2],
                                           (1.0 - all_toi[i]) * verts0[3] + all_toi[i] * verts1[3] };
            double cross_norm = ((vt[0] - vt[1]).cross(vt[2] - vt[3])).norm();
            if (cross_norm < std::numeric_limits<double>::epsilon())
                continue;

            // If not parallel and not wrong side, it's a hit
            actually_hit = true;
            t_impact = all_toi[i];
            break;
        }

        if (!actually_hit) {
            return -1;
        } else {
            return 1;
        }
    } // end test wrong side

    assert(std::isfinite(t_impact));
    return hit;
}

template<>
int
NarrowPhaseCTCD<double, 2>::query_ccd_ee(const Eigen::Vector2d*, const Eigen::Vector2d*, const double&, bool, double&)
{
    return 0;
}

template<>
int
NarrowPhaseCTCD<float, 3>::query_ccd_ee(const Eigen::Vector3f* verts0_,
                                        const Eigen::Vector3f* verts1_,
                                        const float& eta,
                                        bool test_vv_and_ve,
                                        float& t_impact)
{
    Eigen::Vector3d verts0[4], verts1[4];
    for (int i = 0; i < 4; ++i) {
        verts0[i] = verts0_[i].cast<double>();
        verts1[i] = verts1_[i].cast<double>();
    }
    double t = t_impact;
    int ret = NarrowPhaseCTCD<double, 3>::query_ccd_ee(verts0, verts1, eta, test_vv_and_ve, t);
    t_impact = t;
    return ret;
}

template<>
int
NarrowPhaseCTCD<float, 2>::query_ccd_ee(const Eigen::Vector2f*, const Eigen::Vector2f*, const float&, bool, float&)
{
    return 0;
}

// ---------------------------------------------------------
//	Query CCD with ACCD
// ---------------------------------------------------------

template<typename T, int DIM>
int
NarrowPhaseACCD<T, DIM>::query_ccd_vf(const Eigen::Vector<T, DIM>* verts0,
                                      const Eigen::Vector<T, DIM>* verts1,
                                      const T& eta,
                                      T& t_impact)
{
    constexpr bool is_vf = true;
    bool hit = NarrowPhaseACCD<T, DIM>::additive_ccd(verts0, verts1, eta, is_vf, t_impact);
    return int(hit);
}

template<typename T, int DIM>
int
NarrowPhaseACCD<T, DIM>::query_ccd_ee(const Eigen::Vector<T, DIM>* verts0,
                                      const Eigen::Vector<T, DIM>* verts1,
                                      const T& eta,
                                      T& t_impact)
{
    constexpr bool is_vf = false;
    bool hit = NarrowPhaseACCD<T, DIM>::additive_ccd(verts0, verts1, eta, is_vf, t_impact);
    return int(hit);
}

template<typename T, int DIM>
bool
NarrowPhaseACCD<T, DIM>::additive_ccd(const Eigen::Vector<T, DIM>* verts0,
                                      const Eigen::Vector<T, DIM>* verts1,
                                      const T& eta, // gap
                                      bool is_vf,
                                      T& t_impact)
{
    using namespace Eigen;
    constexpr int ns = DIM + 1; // size of stencil
    T s = 1;                    // scaling factor
    T xsi = eta;                // minimal sep
    T t_c = 1;                  // global min t (for line search)

    // Displacement vectors and current x
    Eigen::Vector<T, DIM> p[ns];
    Eigen::Vector<T, DIM> x[ns];
    Eigen::Vector<T, DIM> p_bar = Eigen::Vector<T, DIM>::Zero();
    for (int i = 0; i < ns; ++i) {
        p[i] = verts1[i] - verts0[i];
        x[i] = verts0[i];
        p_bar += T(1) / T(ns) * p[i];
    }

    for (int i = 0; i < ns; ++i) {
        p[i] -= p_bar;
    }

    T l_p_first = 0;
    T l_p_second = 0;
    for (int i = 0; i < ns; ++i) {
        bool is_first = is_vf ? i == 0 : i <= 1;
        if (is_first) {
            l_p_first = std::max(l_p_first, p[i].norm());
        } else {
            l_p_second = std::max(l_p_second, p[i].norm());
        }
    }
    T l_p = l_p_first + l_p_second;
    if (l_p <= 0) {
        return false;
    }

    // We use dist instead of squared dist (line 8 of Alg 1.)
    T d = pair_distance(x, is_vf);
    T g = s * (d * d - xsi * xsi) / (d - xsi);
    t_impact = 0;
    T t_l = (1 - s) * (d * d - xsi * xsi) / ((d + xsi) * l_p);

    // Before the loop check if we're already active
    if (d <= xsi) {
        return true;
    }

    int iter = 0;
    int max_iter = 1000;
    for (; iter < max_iter; ++iter) {
        for (int i = 0; i < ns; ++i) {
            x[i] = x[i] + t_l * p[i];
        }

        d = pair_distance(x, is_vf);
        T eps = (d * d - xsi * xsi) / (d + xsi);
        if (t_impact > 0 && iter > 0 && eps < g) {
            break;
        }

        t_impact += t_l;
        if (t_impact > t_c) {
            return false;
        }

        t_l = 0.9 * (d * d - xsi * xsi) / ((d + xsi) * l_p);
    }

    return true;
}

template<typename T, int DIM>
T
NarrowPhaseACCD<T, DIM>::pair_distance(const Eigen::Vector<T, DIM>* v, bool is_vf)
{
    using Vec3 = Eigen::Matrix<T, 3, 1>;
    using Vec2 = Eigen::Matrix<T, 2, 1>;

    if (DIM == 2) {
        Vec2 v0 = v[0].template head<2>();
        Vec2 v1 = v[1].template head<2>();
        Vec2 v2 = v[2].template head<2>();
        Vec2 pt = mcl::ccd::point_on_edge<T>(v0, v1, v2);
        return (v0 - pt).norm();
    } else {
        Vec3 v0 = v[0].template head<3>();
        Vec3 v1 = v[1].template head<3>();
        Vec3 v2 = v[2].template head<3>();
        Vec3 v3 = v[3].template head<3>();
        if (is_vf) {
            Vec3 pt = mcl::ccd::point_on_triangle<T>(v0, v1, v2, v3);
            return (v0 - pt).norm();
        } else {
            Eigen::Matrix<T, 4, 1> b = Eigen::Matrix<T, 4, 1>::Zero();
            Vec3 eed = mcl::ccd::edge_to_edge(v0, v1, v2, v3, b);
            return eed.norm();
        }
    }
    return -1;
}

// ---------------------------------------------------------
//	Defines
// ---------------------------------------------------------

template class mcl::ccd::NarrowPhase<double>;
template class mcl::ccd::NarrowPhase<float>;

template class mcl::ccd::NarrowPhaseCTCD<double, 2>;
template class mcl::ccd::NarrowPhaseCTCD<float, 2>;
template class mcl::ccd::NarrowPhaseCTCD<double, 3>;
template class mcl::ccd::NarrowPhaseCTCD<float, 3>;

template class mcl::ccd::NarrowPhaseACCD<double, 2>;
template class mcl::ccd::NarrowPhaseACCD<float, 2>;
template class mcl::ccd::NarrowPhaseACCD<double, 3>;
template class mcl::ccd::NarrowPhaseACCD<float, 3>;

} // end namespace ccd
} // end namespace mcl
