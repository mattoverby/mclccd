// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_COLLISIONPAIR_HPP
#define MCL_COLLISIONPAIR_HPP 1

#include <Eigen/Core>
#include <type_traits> // is_same
#include <vector>

namespace mcl
{

// CollisionPair type
enum
{
    COLLISIONPAIR_INVALID,
    COLLISIONPAIR_VF,
    COLLISIONPAIR_EE,
    COLLISIONPAIR_NUM
};

template<typename T>
class CollisionPair
{
public:
    // vf = [v, f0, f1, f2], ee = [p0, p1, q0, q1]
    Eigen::Vector4i stencil; // unused set to -1 (i.e., DIM=2)
    int type; // vf/ee/ae/etc...
    T toi; // time of impact

    CollisionPair();
    virtual ~CollisionPair() {}

    virtual bool operator<(const CollisionPair &c) const; // sorts by toi then stencil
    virtual bool operator==(const CollisionPair &c) const; // only check stencil (does not check ee symmetry)
    virtual bool operator!=(const CollisionPair &c) const { return !(c==*this); }

    virtual std::string hash_str() const; // readable stencil
    virtual std::size_t hash() const;
};

} // ns mcl

// Hash function for using unordered_map/set
namespace std {
    template<typename T> struct hash<mcl::CollisionPair<T>> {
        std::size_t operator()(const mcl::CollisionPair<T>& c) const { return c.hash(); }
    };
}

#endif
