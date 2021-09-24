// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "BVHLeaf.hpp"

namespace mcl
{

template<typename T, int DIM>
BoundingBox<T,DIM> BoundingBox<T,DIM>::merged(const BoundingBox& b) const
{
    BoundingBox<T,DIM> bb;
    bb.active = this->active || b.active;
    bb.extend(this->min());
    bb.extend(this->max());
    bb.extend(b.min());
    bb.extend(b.max());
    bb.t0.extend(this->t0);
    bb.t0.extend(b.t0);
    bb.t1.extend(this->t1);
    bb.t1.extend(b.t1);
    return bb;
}

} // ns mcl

template class mcl::BoundingBox<double,3>;
template class mcl::BoundingBox<double,2>;
template class mcl::BoundingBox<float,3>;
template class mcl::BoundingBox<float,2>;
template class mcl::BVHLeaf<double,3>;
template class mcl::BVHLeaf<double,2>;
template class mcl::BVHLeaf<float,3>;
template class mcl::BVHLeaf<float,2>;