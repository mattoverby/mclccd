// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_COLLISIONPAIR_HPP
#define MCL_COLLISIONPAIR_HPP 1

#include <Eigen/Core>
#include <type_traits> // is_same

namespace mcl
{

// CollisionPair type
enum
{
    COLLISIONPAIR_INVALID,
    COLLISIONPAIR_VF,
    COLLISIONPAIR_EE,
    COLLISIONPAIR_AIR, // air-element
    COLLISIONPAIR_NUM
};

template<typename T, int DIM>
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

    // Compute signed barycoords for the pair at time t: vf=[1,-f0,-f1,-f2], ee=[p0,p1,-q0,-q1]
    // Does not deduce argument, e.g. if X0, X1 are MatrixXd:
    // Vector4d = pair.barys<MatrixXd>(X0, X1, 1);
    template <typename DerivedV>
    Eigen::Matrix<T,DIM+1,1> compute_barys(Eigen::Ref<const DerivedV> X0, Eigen::Ref<const DerivedV> X1, T t, bool clamp_barys=true) const
    {
        // Return zero if something wrong with input
        if (!bool(std::is_same<T,typename DerivedV::Scalar>::value) || !(X0.IsRowMajor == X1.IsRowMajor) ||
            (X0.rows()!=X1.rows()) || (X0.cols()!=DIM) || (X1.cols()!=DIM) )
            { return Eigen::Matrix<T,DIM+1,1>::Zero(); }
        return barys(X0.data(), X1.data(), X0.rows(), X0.IsRowMajor, t, clamp_barys);
    }

protected:

    // Kernel functions for barycoords
    Eigen::Matrix<T,DIM+1,1> compute_barys(const T *X0, const T *X1, int nx, bool rowmajor, T t, bool clamp_barys) const;

};

} // ns mcl

// Hash function for using unordered_map/set
namespace std {
    template<typename T, int DIM> struct hash<mcl::CollisionPair<T,DIM>> {
        std::size_t operator()(const mcl::CollisionPair<T,DIM>& c) const { return c.hash(); }
    };
}

#endif
