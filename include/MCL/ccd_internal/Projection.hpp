// Copyright Matt Overby 2021.
// Distributed under the MIT License.

// Copied from https://github.com/mattoverby/mclgeom, March 2022.
// Added ccd_internal.

#ifndef MCL_PROJECTION_CCD_INTERNAL_HPP
#define MCL_PROJECTION_CCD_INTERNAL_HPP 1

#include <Eigen/Core>

namespace mcl
{
namespace ccd_internal
{

// Projection on Triangle
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_triangle(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &p1, const Eigen::Matrix<T,3,1> &p2, const Eigen::Matrix<T,3,1> &p3);

// Projection on Sphere
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_sphere(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &center, const T &rad);

// Projection on a Box
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_box(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &bmin, const Eigen::Matrix<T,3,1> &bmax);

// Project a point on to a plane
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_plane(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &plane_norm, const Eigen::Matrix<T,3,1> &plane_pt);

// Projection on an edge
template <typename T> 
static inline Eigen::Matrix<T,2,1> point_on_edge(const Eigen::Matrix<T,2,1> &point, const Eigen::Matrix<T,2,1> &p1, const Eigen::Matrix<T,2,1> &p2);

// Projection on an edge (3D)
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_edge(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &p1, const Eigen::Matrix<T,3,1> &p2);

// Computes the shortest vector from p toward q, as well as
// the barycentric coordinates of the closest points.
// The distance between the segments is the norm of this vector.
// Source: https://github.com/evouga/collisiondetection, license: public domain.
template <typename T>
static inline Eigen::Matrix<T,3,1> edge_to_edge(
	const Eigen::Matrix<T,3,1> &p0,
	const Eigen::Matrix<T,3,1> &p1,
	const Eigen::Matrix<T,3,1> &q0,
	const Eigen::Matrix<T,3,1> &q1,
	Eigen::Matrix<T,4,1> &bary);


//
//	Implementation
//

// I do not know the source of the function.
// If anyone knows, please tell me!
template <typename T>
Eigen::Matrix<T,3,1> point_on_triangle(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &p1, const Eigen::Matrix<T,3,1> &p2, const Eigen::Matrix<T,3,1> &p3)
{
	auto clamp_zero_one = [](const T &val){ return val < 0 ? 0 : (val > 1 ? 1 : val); };

	Eigen::Matrix<T,3,1> edge0 = p2 - p1;
	Eigen::Matrix<T,3,1> edge1 = p3 - p1;
	Eigen::Matrix<T,3,1> v0 = p1 - point;

	T a = edge0.dot(edge0);
	T b = edge0.dot(edge1);
	T c = edge1.dot(edge1);
	T d = edge0.dot(v0);
	T e = edge1.dot(v0);
	T det = a*c - b*b;
	T s = b*e - c*d;
	T t = b*d - a*e;

	const T zero(0);
	const T one(1);

	if ( s + t < det ) {
		if ( s < zero ) {
		    if ( t < zero ) {
			if ( d < zero ) {
			    s = clamp_zero_one( -d/a );
			    t = zero;
			}
			else {
			    s = zero;
			    t = clamp_zero_one( -e/c );
			}
		    }
		    else {
			s = zero;
			t = clamp_zero_one( -e/c );
		    }
		}
		else if ( t < zero ) {
		    s = clamp_zero_one( -d/a );
		    t = zero;
		}
		else {
		    T invDet = one / det;
		    s *= invDet;
		    t *= invDet;
		}
	}
	else {
		if ( s < zero ) {
		    T tmp0 = b+d;
		    T tmp1 = c+e;
		    if ( tmp1 > tmp0 ) {
			T numer = tmp1 - tmp0;
			T denom = a-T(2)*b+c;
			s = clamp_zero_one( numer/denom );
			t = one-s;
		    }
		    else {
			t = clamp_zero_one( -e/c );
			s = zero;
		    }
		}
		else if ( t < zero ) {
		    if ( a+d > b+e ) {
			T numer = c+e-b-d;
			T denom = a-T(2)*b+c;
			s = clamp_zero_one( numer/denom );
			t = one-s;
		    }
		    else {
			s = clamp_zero_one( -e/c );
			t = zero;
		    }
		}
		else {
		    T numer = c+e-b-d;
		    T denom = a-T(2)*b+c;
		    s = clamp_zero_one( numer/denom );
		    t = one - s;
		}
	}

	return (p1 + edge0*s + edge1*t);

} // end project triangle


template <typename T>
Eigen::Matrix<T,3,1> point_on_sphere(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &center, const T &rad)
{
	Eigen::Matrix<T,3,1> dir = point-center;
	dir.normalize();
	return (center + dir*rad);
} // end project sphere


template <typename T>
Eigen::Matrix<T,3,1> point_on_box(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &bmin, const Eigen::Matrix<T,3,1> &bmax)
{
	// Loops through axes and moves point to nearest surface
	Eigen::Matrix<T,3,1> x = point;
	T dx = std::numeric_limits<T>::max();
	for (int i=0; i<3; ++i)
	{
		T dx_max = std::abs(bmax[i]-point[i]);
		T dx_min = std::abs(bmin[i]-point[i]);
		if(dx_max < dx)
		{
			x = point;
			x[i] = bmax[i];
			dx = dx_max;
		}
		if (dx_min < dx)
		{
			x = point;
			x[i] = bmin[i];
			dx = dx_min;
		}
	}
	return x;
} // end project box


template <typename T>
Eigen::Matrix<T,3,1> point_on_plane(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &plane_norm, const Eigen::Matrix<T,3,1> &plane_pt)
{
    T d = -1 * plane_norm.dot(point-plane_pt);
    Eigen::Matrix<T,3,1> t_vec = plane_norm * d;
    return point + t_vec;
}


template <typename T>
static Eigen::Matrix<T,2,1> point_on_edge(const Eigen::Matrix<T,2,1> &p, const Eigen::Matrix<T,2,1> &e0, const Eigen::Matrix<T,2,1> &e1)
{
	Eigen::Matrix<T,2,1> e = (e1-e0);
	T e_len2 = e.dot(e);
	if(e_len2 <= 0.0) { return e0; } // zero length edge
	Eigen::Matrix<T,2,1> pe0 = (p-e0);
	T t = pe0.dot(e)/e_len2;
	if (t < 0.0){ return e0; }
	else if (t > 1.0) { return e1; }
	return e0 + t * e;
}


// Ericson, Real-Time Collision Detection
template <typename T>
static Eigen::Matrix<T,3,1> point_on_edge(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &p1, const Eigen::Matrix<T,3,1> &p2)
{
	Eigen::Matrix<T,3,1> ab = p2-p1;
	T t = ab.dot(point-p1);
	Eigen::Matrix<T,3,1> d = point; // result

	// c projects outside the [a,b] interval, on the a side; clamp to a
	if (t <= 0)
	{
		t = 0;
		d = p1;
	}
	else
	{
		T denom = ab.dot(ab); // Always nonnegative since denom = ||ab|| ∧ 2
		// c projects outside the [a,b] interval, on the b side; clamp to b
		if (t >= denom)
		{
			t = 1;
			d = p2;
		}
		// c projects inside the [a,b] interval; must do deferred divide now
		else
		{
			t = t/denom;
			d = p1+t*ab;
		}
	}
	return d;
}


template <typename T>
static Eigen::Matrix<T,3,1> edge_to_edge(
    const Eigen::Matrix<T,3,1> &p0,
    const Eigen::Matrix<T,3,1> &p1,
    const Eigen::Matrix<T,3,1> &q0,
    const Eigen::Matrix<T,3,1> &q1,
    Eigen::Matrix<T,4,1> &bary)
{
	typedef Eigen::Matrix<T,3,1> eteVec3;
	auto clamp01 = [](T x)
	{
		return std::min(T(1), std::max(x, T(0)));
	};

	eteVec3 d1 = p1-p0;
	eteVec3 d2 = q1-q0;
	eteVec3 r = p0-q0;
	T a = d1.squaredNorm();
	T e = d2.squaredNorm();
	T f = d2.dot(r);
	T s = 0;
	T t = 0;
	T c = d1.dot(r);
	T b = d1.dot(d2);
	T denom = a*e-b*b;
	if (denom != T(0)) {
		s = clamp01((b*f-c*e)/denom);
	}
	else { // Parallel/degenerate edges
		s = 0;
	}
	T tnom = b*s + f;
	if(tnom < 0 || e == 0)
	{
		t = 0;
		if(a == 0) {
			s = 0;
		} else {
			s = clamp01(-c/a);
		}
	}
	else if(tnom > e)
	{
		t = 1.0;
		if(a == 0) {
			s = 0;
		} else {
			s = clamp01((b-c)/a);
		}
	}
	else {
		t = tnom/e;
	}

	eteVec3 c1 = p0 + s*d1;
	eteVec3 c2 = q0 + t*d2;
	bary[0] = 1-s;
	bary[1] = s;
	bary[2] = 1-t;
	bary[3] = t;
	return c2-c1;
}

} // end namespace ccd_internal
} // end namespace mcl

#endif
