// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "CollisionPair.hpp"
#include "ccd_internal/Assert.hpp"
#include "ccd_internal/Distance.hpp"

#include <sstream>
#include <functional> // hash

namespace mcl
{

namespace collisionpair_internal
{
    template <typename T, int DIM> inline Eigen::Matrix<T,DIM,1> get_vec_cm(const T *X,  int rows, int idx)
    {
        Eigen::Matrix<T,DIM,1> r;
        for (int i=0; i<DIM; ++i) { r[i] = X[i*rows+idx]; }
        return r;
    }
    template <typename T, int DIM> inline Eigen::Matrix<T,DIM,1> get_vec_rm(const T *X, int idx)
    {
        Eigen::Matrix<T,DIM,1> r;
        for (int i=0; i<DIM; ++i) { r[i] = X[idx*DIM+i]; }
        return r;
    }
}

template<typename T, int DIM>
CollisionPair<T,DIM>::CollisionPair() :
    stencil(-1,-1,-1,-1),
    type(COLLISIONPAIR_INVALID),
    toi(-1)
    {}

template<typename T, int DIM>
bool CollisionPair<T,DIM>::operator<(const CollisionPair<T,DIM> &c) const
{
    // First, sort by TOI
    if (toi < c.toi) { return true; }
    if (toi > c.toi) { return false; }

    // If TOI the same, sort by stencil
    for (int i=0; i<4; ++i)
    {
        if (stencil[i] < c.stencil[i]) { return true; }
        if (stencil[i] > c.stencil[i]) { return false; }
    }

    // Finally, sort by hash
    return hash() < c.hash();
}

template<typename T, int DIM>
bool CollisionPair<T,DIM>::operator==(const CollisionPair<T,DIM> &c) const
{
    for (int i=0; i<4; ++i)
        if (stencil[i] != c.stencil[i])
            return false;

    mclAssert(hash() == c.hash(), "hash==hash, but sten!=sten");
    return true;
}

template<typename T, int DIM>
std::string CollisionPair<T,DIM>::hash_str() const
{
    std::stringstream ss;
    if (type == COLLISIONPAIR_VF && DIM == 2)
    {
        ss << "vf " << stencil[0] << ' ';
        if (stencil[1] < stencil[2]) { ss << stencil[1] << ' ' << stencil[2]; }
        else { ss << stencil[2] << ' ' << stencil[1]; }
    }
    else if (type == COLLISIONPAIR_VF && DIM == 3)
    {
        Eigen::Vector3i s(stencil[1], stencil[2], stencil[3]);
        if (s[0] > s[1]) { std::swap(s[0], s[1]); }
        if (s[1] > s[2]) { std::swap(s[1], s[2]); }
        if (s[0] > s[1]) { std::swap(s[0], s[1]); }
        ss << "vf " << stencil[0] << ' ' << s[0] << ' ' << s[1] << ' ' << s[2];

    }
    else if (type == COLLISIONPAIR_EE && DIM == 3)
    {
        Eigen::Vector4i s = stencil;
        if (s[1] < s[0]) { std::swap(s[0], s[1]); }
        if (s[3] < s[2]) { std::swap(s[2], s[3]); }
        if (s[2] < s[0]) { std::swap(s[0], s[2]); std::swap(s[1], s[3]); }
        ss << "ee " << s[0] << ' ' << s[1] << ' ' << s[2] << ' ' << s[3];
    }
    else if (type == COLLISIONPAIR_AIR)
    {
        Eigen::Vector4i s = stencil;
        std::sort(s.data(), s.data()+4);
        ss << "ae " << s[0] << ' ' << s[1] << ' ' << s[2];
        if (DIM==3) { ss << ' ' << s[3]; }
    }

    return ss.str();
}

template<typename T, int DIM>
std::size_t CollisionPair<T,DIM>::hash() const
{
    return std::hash<std::string>{}(hash_str());
}

template<typename T, int DIM>
Eigen::Matrix<T,DIM+1,1> CollisionPair<T,DIM>::compute_barys(const T *X0, const T *X1, int nx, bool rowmajor, T t, bool clamp_barys) const
{
    using namespace Eigen;
    using namespace collisionpair_internal;


    if (type == COLLISIONPAIR_AIR)
    {
        return Eigen::Matrix<T,DIM+1,1>::Ones();
    }

    Matrix<T,DIM,1> xt[DIM+1];
    for (int i=0; i<DIM+1; ++i)
    {
        Matrix<T,DIM,1> x0 = rowmajor ? get_vec_rm<T,DIM>(X0, stencil[i]) : get_vec_cm<T,DIM>(X0, nx, stencil[i]);
        Matrix<T,DIM,1> x1 = rowmajor ? get_vec_rm<T,DIM>(X1, stencil[i]) : get_vec_cm<T,DIM>(X1, nx, stencil[i]);
        xt[i] = (1.0-t) * x0 + t * x1;
    }

    if (DIM==2 && type == COLLISIONPAIR_VF)
    {
        Matrix<T,DIM,1> v0 = xt[2] - xt[1];
        Matrix<T,DIM,1> v2 = xt[0] - xt[1];
        T d00 = v0.dot(v0);
        T d20 = v2.dot(v0);
        T invDenom = 1.0 / d00;
        Matrix<T,DIM+1,1> r = Matrix<T,DIM+1,1>::Zero();
        r[0] = 1;
        r[1] = (d00 - d20)*invDenom;
        r[2] = 1.0 - r[0];
        if (clamp_barys)
        {
            r[1] = std::clamp(r[1], T(0), T(1));
            r[2] = std::clamp(r[2], T(0), T(1));
        }
        r[1] *= -1;
        r[2] *= -1;
        return r;
    }
    else if (DIM==3)
    {
        Matrix<double,DIM+1,1> r = Matrix<double,DIM+1,1>::Zero();
        if (type == COLLISIONPAIR_VF)
        {
            r[0] = 1;
            ctcd::vertexFaceDistance(
                xt[0].template cast<double>().head(3),
                xt[1].template cast<double>().head(3),
                xt[2].template cast<double>().head(3),
                xt[3].template cast<double>().head(3),
                r[1], r[2], r[3], clamp_barys);
                r[1] *= -1;
                r[2] *= -1;
                r[3] *= -1;
        }
        else if (type == COLLISIONPAIR_EE)
        {
            ctcd::edgeEdgeDistance(
                xt[0].template cast<double>().head(3),
                xt[1].template cast<double>().head(3),
                xt[2].template cast<double>().head(3),
                xt[3].template cast<double>().head(3),
                r[0], r[1], r[2], r[3], clamp_barys);
            r[2] *= -1;
            r[3] *= -1;
        }
        return r.template cast<T>();
    }
    return Matrix<T,DIM+1,1>::Zero();
}


} // ns mcl


template class mcl::CollisionPair<double,3>;
template class mcl::CollisionPair<double,2>;
template class mcl::CollisionPair<float,3>;
template class mcl::CollisionPair<float,2>;

/*
// This works but is super tedius

typedef mcl::CollisionPair<double,3> CollisionPair3d;
typedef mcl::CollisionPair<double,2> CollisionPair2d;
typedef mcl::CollisionPair<float,3> CollisionPair3f;
typedef mcl::CollisionPair<float,2> CollisionPair2f;

typedef Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> RowMatrixX3d;
typedef Eigen::Matrix<double,Eigen::Dynamic,2,Eigen::RowMajor> RowMatrixX2d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> RowMatrixXd;
typedef Eigen::Matrix<float,Eigen::Dynamic,3,Eigen::RowMajor> RowMatrixX3f;
typedef Eigen::Matrix<float,Eigen::Dynamic,2,Eigen::RowMajor> RowMatrixX2f;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> RowMatrixXf;

// 3D
template Eigen::Vector4d mcl::CollisionPair3d::barys(const Eigen::MatrixBase<Eigen::MatrixXd>&, const Eigen::MatrixBase<Eigen::MatrixXd>&, double) const;
template Eigen::Vector4d mcl::CollisionPair3d::barys(const Eigen::MatrixBase<RowMatrixXd>&, const Eigen::MatrixBase<RowMatrixXd>&, double) const;
template Eigen::Vector4d mcl::CollisionPair3d::barys(const Eigen::MatrixBase<RowMatrixX3d>&, const Eigen::MatrixBase<RowMatrixX3d>&, double) const;
template Eigen::Vector4f mcl::CollisionPair3d::barys(const Eigen::MatrixBase<Eigen::MatrixXf>&, const Eigen::MatrixBase<Eigen::MatrixXf>&, double) const;
template Eigen::Vector4f mcl::CollisionPair3d::barys(const Eigen::MatrixBase<RowMatrixXf>&, const Eigen::MatrixBase<RowMatrixXf>&, double) const;
template Eigen::Vector4f mcl::CollisionPair3d::barys(const Eigen::MatrixBase<RowMatrixX3f>&, const Eigen::MatrixBase<RowMatrixX3f>&, double) const;
*/
