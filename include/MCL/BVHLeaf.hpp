// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_BVHLEAF_HPP
#define MCL_BVHLEAF_HPP 1

#include <Eigen/Geometry>

namespace mcl
{

template<typename T, int DIM>
class BoundingBox : public Eigen::AlignedBox<T,DIM>
{
public:
    bool active;

    // Box at t=0 and t=1
    Eigen::AlignedBox<T,DIM> t0, t1;
    
    BoundingBox() :
        Eigen::AlignedBox<T,DIM>(),
        active(true)
        {}

    virtual ~BoundingBox() {}

    BoundingBox merged(const BoundingBox& b) const;
};

template<typename T, int DIM>
class BVHLeaf
{
public:
    typedef BoundingBox<T,DIM> BoxType;

    // Index of the primitive
    int idx;

    // If rep-triangle, boolean indicator that v/e is represented
    Eigen::Vector3i v, e;

    BoxType box;

    BVHLeaf() :
        idx(-1),
        v(1,1,1),
        e(1,1,1)
        {}

    virtual ~BVHLeaf() {}
};

} // ns mcl

// Used for Eigen's BVH module
namespace Eigen
{
    template<typename T, int DIM>
    static inline typename mcl::BVHLeaf<T,DIM>::BoxType
    bounding_box(const mcl::BVHLeaf<T,DIM> &n)
    {
        return n.box;
    }
}

#endif
