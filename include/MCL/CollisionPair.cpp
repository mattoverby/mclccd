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

template<typename T>
CollisionPair<T>::CollisionPair() :
    stencil(-1,-1,-1,-1),
    type(COLLISIONPAIR_INVALID),
    toi(-1)
    {}

template<typename T>
bool CollisionPair<T>::operator<(const CollisionPair<T> &c) const
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

template<typename T>
bool CollisionPair<T>::operator==(const CollisionPair<T> &c) const
{
    for (int i=0; i<4; ++i)
        if (stencil[i] != c.stencil[i])
            return false;

    mclAssert(hash() == c.hash(), "hash==hash, but sten!=sten");
    return true;
}

template<typename T>
std::string CollisionPair<T>::hash_str() const
{
    bool is_2D = stencil[3] == -1;
    bool is_3D = stencil[3] > -1;
    std::stringstream ss;
    if (type == COLLISIONPAIR_VF && is_2D)
    {
        ss << "vf " << stencil[0] << ' ';
        if (stencil[1] < stencil[2]) { ss << stencil[1] << ' ' << stencil[2]; }
        else { ss << stencil[2] << ' ' << stencil[1]; }
    }
    else if (type == COLLISIONPAIR_VF && is_3D)
    {
        Eigen::Vector3i s(stencil[1], stencil[2], stencil[3]);
        if (s[0] > s[1]) { std::swap(s[0], s[1]); }
        if (s[1] > s[2]) { std::swap(s[1], s[2]); }
        if (s[0] > s[1]) { std::swap(s[0], s[1]); }
        ss << "vf " << stencil[0] << ' ' << s[0] << ' ' << s[1] << ' ' << s[2];

    }
    else if (type == COLLISIONPAIR_EE && is_3D)
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
        if (is_3D) { ss << ' ' << s[3]; }
    }

    return ss.str();
}

template<typename T>
std::size_t CollisionPair<T>::hash() const
{
    return std::hash<std::string>{}(hash_str());
}

template<typename T>
std::vector<T> CollisionPair<T>::compute_barys(const T *X0, const T *X1, int nx, int dim, bool rowmajor, T t, bool clamp_barys) const
{
    using namespace Eigen;
    using namespace collisionpair_internal;
    mclAssert(dim == 2 || dim == 3);

    if (type == COLLISIONPAIR_AIR)
    {
        return std::vector<T>(dim+1, 1);
    }

    if (dim==2 && type == COLLISIONPAIR_VF)
    {
        typedef Matrix<T,2,1> VecType;
        std::vector<VecType> xt(dim+1);
        for (int i=0; i<dim+1; ++i)
        {
            VecType x0 = rowmajor ? get_vec_rm<T,2>(X0, stencil[i]) : get_vec_cm<T,2>(X0, nx, stencil[i]);
            VecType x1 = rowmajor ? get_vec_rm<T,2>(X1, stencil[i]) : get_vec_cm<T,2>(X1, nx, stencil[i]);
            xt[i] = (1.0-t) * x0 + t * x1;
        }

        VecType v0 = xt[2] - xt[1];
        VecType v2 = xt[0] - xt[1];
        T d00 = v0.dot(v0);
        T d20 = v2.dot(v0);
        T invDenom = 1.0 / d00;
        std::vector<T> r(dim+1, 0);
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
    else if (dim==3)
    {
        typedef Matrix<T,3,1> VecType;
        std::vector<VecType> xt(dim+1);
        for (int i=0; i<dim+1; ++i)
        {
            VecType x0 = rowmajor ? get_vec_rm<T,3>(X0, stencil[i]) : get_vec_cm<T,3>(X0, nx, stencil[i]);
            VecType x1 = rowmajor ? get_vec_rm<T,3>(X1, stencil[i]) : get_vec_cm<T,3>(X1, nx, stencil[i]);
            xt[i] = (1.0-t) * x0 + t * x1;
        }
        Vector4d r = Vector4d::Zero();

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
        std::vector<T> rt(4);
        rt[0] = r[0];
        rt[1] = r[1];
        rt[2] = r[2];
        rt[3] = r[3];
        return rt;
    }
    return std::vector<T>();
}


} // ns mcl


template class mcl::CollisionPair<double>;
template class mcl::CollisionPair<float>;