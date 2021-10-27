// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "NarrowPhase.hpp"
#include "ccd_internal/Assert.hpp"
#include "ccd_internal/CTCD.hpp"
#include "ccd_internal/Distance.hpp"
#include "ccd_internal/tt_isect.hpp"

#include <limits>

// I need to fix tabbing...

namespace mcl
{

// ---------------------------------------------------------
//	Hit Wrong Side
// ---------------------------------------------------------

template<>
bool NarrowPhase<double,2>::hit_wrong_side_vf(
    const Eigen::Vector2d *verts0,
    const Eigen::Vector2d *verts1,
    double t)
{
    using namespace Eigen;

	// Resting contact, relative velocity=0
	// and so there is no way to tell.
	if (std::abs(t) <= 0.0) { return false; }
	const Vector2d &q0start = verts0[0];
	const Vector2d &q1start = verts0[1];
	const Vector2d &q2start = verts0[2];
	const Vector2d &q0end = verts1[0];
	const Vector2d &q1end = verts1[1];
	const Vector2d &q2end = verts1[2];

    Vector2d p = q0start*(1.0-t) + q0end*t;
    Vector2d e0 = q1start*(1.0-t) + q1end*t;
    Vector2d e1 = q2start*(1.0-t) + q2end*t;
    Vector2d n = Vector2d( e1[1]-e0[1], -(e1[0]-e0[0]) ).normalized();
    Vector2d barys = Vector2d((e1-p).norm(), (e0-p).norm()) / (e0-e1).norm();

	Vector2d e0bary = barys[0]*q1start + barys[1]*q2start;
	Vector2d e1bary = barys[0]*q1end + barys[1]*q2end;
	Vector2d apex_vel = q0end-q0start;
	Vector2d edge_vel = e1bary-e0bary;
	Vector2d v = apex_vel - edge_vel;
	double vel_v_dot_e = n.dot(v);
	return vel_v_dot_e > 0.0;
}

template<>
bool NarrowPhase<double,3>::hit_wrong_side_vf(
	const Eigen::Vector3d *verts0,
	const Eigen::Vector3d *verts1,
	double t)
{
    using namespace Eigen;

	// Resting contact, relative velocity=0
	// and so there is no way to tell.
	if (std::abs(t) <= 0.0) { return false; }
	std::vector<Vector3d> xt = {
		verts0[0]*(1.0-t) + verts1[0]*t,
		verts0[1]*(1.0-t) + verts1[1]*t,
		verts0[2]*(1.0-t) + verts1[2]*t,
		verts0[3]*(1.0-t) + verts1[3]*t };

	Vector3d n = (xt[2]-xt[1]).cross(xt[3]-xt[1]);
    if (n.squaredNorm() <= 0) { return false; }
    n.stableNormalize();

	Vector3d barys = Vector3d::Zero();
	ctcd::vertexFaceDistance(xt[0],xt[1],xt[2],xt[3],barys[0],barys[1],barys[2]);

	Vector3d apex_vel = ( verts1[0]-verts0[0] );
	Vector3d face_pt0 = (barys[0]*verts0[1] + barys[1]*verts0[2] + barys[2]*verts0[3]);
	Vector3d face_pt1 = (barys[0]*verts1[1] + barys[1]*verts1[2] + barys[2]*verts1[3]);
	Vector3d face_vel = face_pt1 - face_pt0;
	Vector3d v = apex_vel - face_vel;
	double vel_v_dot_f = n.dot(v);
	bool hit_wrong_side = vel_v_dot_f > 0.0;
	return hit_wrong_side;
}

template<>
bool NarrowPhase<float,2>::hit_wrong_side_vf(
	const Eigen::Vector2f *verts0_,
	const Eigen::Vector2f *verts1_,
	float t)
{
    Eigen::Vector2d verts0[3], verts1[3];
    for (int i=0; i<3; ++i)
    {
        verts0[i] = verts0_[i].cast<double>();
        verts1[i] = verts1_[i].cast<double>();
    }
    return NarrowPhase<double,2>::hit_wrong_side_vf(verts0, verts1, t);
}

template<>
bool NarrowPhase<float,3>::hit_wrong_side_vf(
	const Eigen::Vector3f *verts0_,
	const Eigen::Vector3f *verts1_,
	float t)
{
    Eigen::Vector3d verts0[4], verts1[4];
    for (int i=0; i<4; ++i)
    {
        verts0[i] = verts0_[i].cast<double>();
        verts1[i] = verts1_[i].cast<double>();
    }
    return NarrowPhase<double,3>::hit_wrong_side_vf(verts0, verts1, t);
}


template <>
bool NarrowPhase<double,3>::query_ray_box(
	const Eigen::Vector3d &p_x0, const Eigen::Vector3d &p_x1,
	const Eigen::Vector3d &bmin, const Eigen::Vector3d &bmax)
{
    using namespace Eigen;

	Vector3d dir = p_x1 - p_x0;
	const Vector3d &origin = p_x0;
	double t0 = 0;
	double t1 = dir.norm();

    // Starts inside box
    AlignedBox<double,3> box;
    box.extend(bmin); box.extend(bmax);
    if (box.contains(p_x0) || box.contains(p_x1)) { return true; }

	dir.normalize();
	typedef Matrix<double,1,3>  RowVector3S;
	const RowVector3S inv_dir( 1./dir(0),1./dir(1),1./dir(2));
	const std::vector<bool> sign = { inv_dir(0)<0, inv_dir(1)<0, inv_dir(2)<0};
	// http://people.csail.mit.edu/amy/papers/box-jgt.pdf
	// "An Efficient and Robust Ray–Box Intersection Algorithm"
	double tymin, tymax, tzmin, tzmax;
	std::vector<RowVector3S> bounds = {bmin, bmax};
	double tmin = ( bounds[sign[0]](0) - origin(0)) * inv_dir(0);
	double tmax = ( bounds[1-sign[0]](0) - origin(0)) * inv_dir(0);
	tymin = (bounds[sign[1]](1)   - origin(1)) * inv_dir(1);
	tymax = (bounds[1-sign[1]](1) - origin(1)) * inv_dir(1);
	if ((tmin > tymax) || (tymin > tmax)) { return false; }
	if (tymin > tmin) { tmin = tymin; }
	if (tymax < tmax) { tmax = tymax; }
	tzmin = (bounds[sign[2]](2) - origin(2)) * inv_dir(2);
	tzmax = (bounds[1-sign[2]](2) - origin(2)) * inv_dir(2);
	if ((tmin > tzmax) || (tzmin > tmax)) { return false; }
	if (tzmin > tmin) { tmin = tzmin; }
	if (tzmax < tmax) { tmax = tzmax; }
	if(!((tmin < t1) && (tmax > t0))) { return false; }
	return true;
}

template <>
bool NarrowPhase<double,2>::query_ray_box(
	const Eigen::Vector2d &p_x0, const Eigen::Vector2d &p_x1,
	const Eigen::Vector2d &bmin, const Eigen::Vector2d &bmax )
{
	// TODO: this function
	(void)(p_x0); (void)(p_x1); (void)(bmin); (void)(bmax);
	return true;
}

template <>
bool NarrowPhase<float,2>::query_ray_box(
	const Eigen::Vector2f &p_x0, const Eigen::Vector2f &p_x1,
	const Eigen::Vector2f &bmin, const Eigen::Vector2f &bmax )
{
    return NarrowPhase<double,2>::query_ray_box(
        p_x0.cast<double>(), p_x1.cast<double>(),
        bmin.cast<double>(), bmax.cast<double>());
}

template <>
bool NarrowPhase<float,3>::query_ray_box(
	const Eigen::Vector3f &p_x0, const Eigen::Vector3f &p_x1,
	const Eigen::Vector3f &bmin, const Eigen::Vector3f &bmax )
{
    return NarrowPhase<double,3>::query_ray_box(
        p_x0.cast<double>(), p_x1.cast<double>(),
        bmin.cast<double>(), bmax.cast<double>());
}

// ---------------------------------------------------------
//	Query CCD
// ---------------------------------------------------------


template<>
int NarrowPhase<double,2>::query_ccd_vf(
    const Eigen::Vector2d *verts0,
    const Eigen::Vector2d *verts1,
    const double &eta,
    bool test_wrong_side,
    double &t_impact)
{
	// First test is to check AABBs
	Eigen::AlignedBox<double,2> q_AABB;
	q_AABB.extend(verts0[1]);  q_AABB.extend(verts1[1]);
	q_AABB.extend(verts0[2]);  q_AABB.extend(verts1[2]);
	for (int i=0; i<2; ++i)
	{
		q_AABB.min()[i] -= eta;
		q_AABB.max()[i] += eta;
	}

	if (!query_ray_box(verts0[0], verts1[0], q_AABB.min(), q_AABB.max()))
        return 0;

	std::vector<double> all_toi;
	bool ve = CTCD::vertexEdgeCTCD(
		verts0[0], verts0[1], verts0[2],
		verts1[0], verts1[1], verts1[2],
		eta, t_impact, &all_toi);

	// Wrong side collision?
	if (ve && test_wrong_side)
	{
		mclAssert(all_toi.size()>0);

		// Loop t. Find first t that is not wrong-side-collision
		std::sort(all_toi.begin(), all_toi.end());
		bool actually_hit = false;
		int nt = all_toi.size();
		for (int i=0; i<nt; ++i)
		{
			bool wrongside = hit_wrong_side_vf(verts0,verts1,all_toi[i]);
			if(!wrongside)
			{
				actually_hit = true;
				t_impact = all_toi[i];
				break;
			}
		}

		if (!actually_hit) { return -1; }
		else { return 1; }

	} // end test wrong side

	return ve;

} // end query ccd vf 2d


template<>
int NarrowPhase<double,3>::query_ccd_vf(
    const Eigen::Vector3d *verts0,
    const Eigen::Vector3d *verts1,
    const double &eta,
    bool test_wrong_side,
    double &t_impact)
{
    using namespace Eigen;
	typedef AlignedBox<double,3> AABB;

	// First test is to check AABBs
	{
		AABB q_AABB;
		q_AABB.extend(verts0[1]);  q_AABB.extend(verts1[1]);
		q_AABB.extend(verts0[2]);  q_AABB.extend(verts1[2]);
		q_AABB.extend(verts0[3]);  q_AABB.extend(verts1[3]);
		for (int i=0; i<3; ++i)
		{
			q_AABB.min()[i] -= eta;
			q_AABB.max()[i] += eta;
		}
		if (!query_ray_box(verts0[0], verts1[0], q_AABB.min(), q_AABB.max()))
            return 0;
	}

	std::vector<double> all_toi;
    auto ccd_vf = [&]()->bool
    {
        if (CTCD::vertexFaceCTCD(
            verts0[0], verts0[1], verts0[2], verts0[3],
            verts1[0], verts1[1], verts1[2], verts1[3],
            eta, t_impact, &all_toi)) { return true; }

        // Vertex-face edges
        for(int edge=0; edge<3; ++edge)
        {
            if (CTCD::vertexEdgeCTCD(
                verts0[0], verts0[1+(edge%3)], verts0[1+ ((edge+1)%3)],
                verts1[0], verts1[1+(edge%3)], verts1[1+ ((edge+1)%3)],
                eta, t_impact, &all_toi)) { return true; }
        }

        // Vertex-face vertices
        for (int vert=0; vert<3; ++vert)
        {
            if (CTCD::vertexVertexCTCD(
                verts0[0], verts0[1+vert],
                verts1[0], verts1[1+vert],
                eta, t_impact, &all_toi)) { return true; }
        }

        return false;
    };

	int hit = ccd_vf();

	// Wrong side collision?
	if (hit && test_wrong_side)
	{
		mclAssert((int)all_toi.size() > 0);

		// Loop t. Find first t that is not wrong-side-collision
		std::sort(all_toi.begin(), all_toi.end());
		bool actually_hit = false;
		int nt = all_toi.size();
		for (int i=0; i<nt; ++i)
		{
			if (all_toi[i] < 0) { continue; }
			bool wrongside = hit_wrong_side_vf(verts0,verts1,all_toi[i]);
			if(!wrongside)
			{
				actually_hit = true;
				t_impact = all_toi[i];
				break;
			}
		}

		if (!actually_hit) { return -1; }
		else { return 1; }

	} // end test wrong side

	mclAssert(std::isfinite(t_impact));
	return hit;

} // end query ccd vf

template<>
int NarrowPhase<float,2>::query_ccd_vf(
	const Eigen::Vector2f *verts0_,
	const Eigen::Vector2f *verts1_,
    const float &eta,
	bool test_wrong_side,
	float &t_impact)
{
    Eigen::Vector2d verts0[3], verts1[3];
    for (int i=0; i<3; ++i)
    {
        verts0[i] = verts0_[i].cast<double>();
        verts1[i] = verts1_[i].cast<double>();
    }
    double t = t_impact;
    int ret = NarrowPhase<double,2>::query_ccd_vf(verts0, verts1, eta, test_wrong_side, t);
    t_impact = t;
    return ret;
}

template<>
int NarrowPhase<float,3>::query_ccd_vf(
	const Eigen::Vector3f *verts0_,
	const Eigen::Vector3f *verts1_,
    const float &eta,
	bool test_wrong_side,
	float &t_impact)
{
    Eigen::Vector3d verts0[4], verts1[4];
    for (int i=0; i<4; ++i)
    {
        verts0[i] = verts0_[i].cast<double>();
        verts1[i] = verts1_[i].cast<double>();
    }
    double t = t_impact;
    int ret = NarrowPhase<double,3>::query_ccd_vf(verts0, verts1, eta, test_wrong_side, t);
    t_impact = t;
    return ret;
}

template<>
int NarrowPhase<double,3>::query_ccd_ee(
    const Eigen::Vector3d *verts0,
    const Eigen::Vector3d *verts1,
    const double &eta,
    bool test_vv_and_ve,
    double &t_impact)
{
    using namespace Eigen;
	typedef AlignedBox<double,3> AABB;

	// First test is to check AABBs
	AABB p_AABB, q_AABB;
	p_AABB.extend(verts0[0]);  p_AABB.extend(verts1[0]);
	p_AABB.extend(verts0[1]);  p_AABB.extend(verts1[1]);
	q_AABB.extend(verts0[2]);  q_AABB.extend(verts1[2]);
	q_AABB.extend(verts0[3]);  q_AABB.extend(verts1[3]);
	for (int i=0; i<3; ++i)
	{
		p_AABB.min()[i] -= eta;
		p_AABB.max()[i] += eta;
		q_AABB.min()[i] -= eta;
		q_AABB.max()[i] += eta;
	}
	if (!p_AABB.intersects(q_AABB))
		return 0;

	std::vector<double> all_toi;
    auto ccd_ee = [&]()->bool
    {
        if (CTCD::edgeEdgeCTCD(
            verts0[0], verts0[1], verts0[2], verts0[3],
            verts1[0], verts1[1], verts1[2], verts1[3],
            eta, t_impact, &all_toi)) { return true; }

        if (test_vv_and_ve)
        {
            if (CTCD::vertexEdgeCTCD(verts0[0], verts0[2], verts0[3], verts1[0], verts1[2], verts1[3], eta, t_impact, &all_toi))
                return true;

            if (CTCD::vertexEdgeCTCD(verts0[1], verts0[2], verts0[3], verts1[1], verts1[2], verts1[3], eta, t_impact, &all_toi))
                return true;

            if (CTCD::vertexEdgeCTCD(verts0[2], verts0[0], verts0[1], verts1[2], verts1[0], verts1[1], eta, t_impact, &all_toi))
                return true;

            if (CTCD::vertexEdgeCTCD(verts0[3], verts0[0], verts0[1], verts1[3], verts1[0], verts1[1], eta, t_impact, &all_toi))
                return true;

            if (CTCD::vertexVertexCTCD(verts0[0], verts0[2], verts1[0], verts1[2], eta, t_impact, &all_toi))
                return true;

            if (CTCD::vertexVertexCTCD(verts0[0], verts0[3], verts1[0], verts1[3], eta, t_impact, &all_toi))
                return true;

            if (CTCD::vertexVertexCTCD(verts0[1], verts0[2], verts1[1], verts1[2], eta, t_impact, &all_toi))
                return true;

            if (CTCD::vertexVertexCTCD(verts0[1], verts0[3], verts1[1], verts1[3], eta, t_impact, &all_toi))
                return true;
        }

        return false;
    };

    int hit = ccd_ee();

    // Collision result okay?
    if (hit)
    {
        mclAssert((int)all_toi.size() > 0);
		std::sort(all_toi.begin(), all_toi.end());

        // Loop t. Find first t that is not parallel
        bool actually_hit = false;
        int nt = all_toi.size();
        for(int i=0; i<nt; ++i)
        {
            if (all_toi[i] < 0) { continue; }

            // Check if the edges are parallel
            // The edgeEdgeCTCD sometimes returns toi for parallel edges!
            std::vector<Vector3d> vt = {
                (1.0-all_toi[i]) * verts0[0] + all_toi[i]*verts1[0],
                (1.0-all_toi[i]) * verts0[1] + all_toi[i]*verts1[1],
                (1.0-all_toi[i]) * verts0[2] + all_toi[i]*verts1[2],
                (1.0-all_toi[i]) * verts0[3] + all_toi[i]*verts1[3] };
            double cross_norm = ((vt[0]-vt[1]).cross(vt[2]-vt[3])).norm();
            if (cross_norm < std::numeric_limits<double>::epsilon())
                continue;

            // If not parallel and not wrong side, it's a hit
            actually_hit = true;
            t_impact = all_toi[i];
            break;
        }

        if (!actually_hit) { return -1; }
        else { return 1; }
    } // end test wrong side

	mclAssert(std::isfinite(t_impact));
	return hit;
}

template<>
int NarrowPhase<double,2>::query_ccd_ee(
    const Eigen::Vector2d*, const Eigen::Vector2d*,
    const double&, bool, double&) { return 0; }

template<>
int NarrowPhase<float,3>::query_ccd_ee(
    const Eigen::Vector3f *verts0_,
    const Eigen::Vector3f *verts1_,
    const float &eta,
    bool test_vv_and_ve,
    float &t_impact)
{
    Eigen::Vector3d verts0[4], verts1[4];
    for (int i=0; i<4; ++i)
    {
        verts0[i] = verts0_[i].cast<double>();
        verts1[i] = verts1_[i].cast<double>();
    }
    double t = t_impact;
    int ret = NarrowPhase<double,3>::query_ccd_ee(verts0, verts1, eta, test_vv_and_ve, t);
    t_impact = t;
    return ret;
}

template<>
int NarrowPhase<float,2>::query_ccd_ee(
    const Eigen::Vector2f*, const Eigen::Vector2f*,
    const float&, bool, float&) { return 0; }

// ---------------------------------------------------------
//	Discrete tests
// ---------------------------------------------------------

template<>
bool NarrowPhase<double,3>::discrete_tri_tri(
    const Eigen::Vector3d *f0, const Eigen::Vector3d *f1)
{
    return tritri::tri_tri_overlap_test_3d(
        (double*)f0[0].data(), (double*)f0[1].data(), (double*)f0[2].data(),
        (double*)f1[0].data(), (double*)f1[1].data(), (double*)f1[2].data());
}

template<>
bool NarrowPhase<float,3>::discrete_tri_tri(
    const Eigen::Vector3f *f0_, const Eigen::Vector3f *f1_)
{
    Eigen::Vector3d f0[3], f1[3];
    for (int i=0; i<3; ++i)
    {
        f0[i] = f0_[i].cast<double>();
        f1[i] = f1_[i].cast<double>();
    }
    return NarrowPhase<double,3>::discrete_tri_tri(f0,f1);
}

template<> bool NarrowPhase<double,2>::discrete_tri_tri(
    const Eigen::Vector3d*, const Eigen::Vector3d*) { return false; }


template<> bool NarrowPhase<float,2>::discrete_tri_tri(
    const Eigen::Vector3f*, const Eigen::Vector3f*) { return false; }


template<> bool NarrowPhase<double,2>::discrete_edge_edge(
    const Eigen::Vector2d *e0, const Eigen::Vector2d *e1)
{
	// From https://stackoverflow.com/a/565282
    using namespace Eigen;
	constexpr double eps = std::numeric_limits<double>::epsilon();

    const Vector2d& p0 = e0[0];
    const Vector2d& p1 = e0[1];
    const Vector2d& q0 = e1[0];
    const Vector2d& q1 = e1[1];

	Vector2d n(q0[0] - p0[0], q0[1] - p0[1]);
	Vector2d r(p1[0] - p0[0], p1[1] - p0[1]);
	Vector2d s(q1[0] - q0[0], q1[1] - q0[1]);
	double rxs = r[0] * s[1] - r[1] * s[0];
	if (std::abs(rxs) < eps) { return false; } // parallel
	double nxr = n[0] * r[1] - n[1] * r[0];
	if (std::abs(nxr) < eps) // collinear
	{
		return ((q0[0] - p0[0] < 0) != (q0[0] - p1[0] < 0)) ||
			((q0[1] - p0[1] < 0) != (q0[1] - p1[1] < 0));
	}
	double nxs = n[0] * s[1] - n[1] * s[0];
	double rxsr = 1.0 / rxs;
	double t = nxs * rxsr;
	double u = nxr * rxsr;
	return (t >= 0) && (t <= 1) && (u >= 0) && (u <= 1);
}

template<> bool NarrowPhase<float,2>::discrete_edge_edge(
    const Eigen::Vector2f *e0_, const Eigen::Vector2f *e1_)
{
    Eigen::Vector2d e0[3], e1[3];
    for (int i=0; i<2; ++i)
    {
        e0[i] = e0_[i].cast<double>();
        e1[i] = e1_[i].cast<double>();
    }
    return NarrowPhase<double,2>::discrete_edge_edge(e0,e1);
}

template<> bool NarrowPhase<double,3>::discrete_edge_edge(
    const Eigen::Vector2d*, const Eigen::Vector2d*) { return false; }

template<> bool NarrowPhase<float,3>::discrete_edge_edge(
    const Eigen::Vector2f*, const Eigen::Vector2f*) { return false; }

// ---------------------------------------------------------
//	Defines
// ---------------------------------------------------------

template class mcl::NarrowPhase<double,2>;
template class mcl::NarrowPhase<float,2>;
template class mcl::NarrowPhase<double,3>;
template class mcl::NarrowPhase<float,3>;

} // ns mcl
