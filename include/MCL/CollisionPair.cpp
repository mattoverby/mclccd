// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "CollisionPair.hpp"
#include "ccd_internal/Assert.hpp"
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

    return ss.str();
}

template<typename T>
std::size_t CollisionPair<T>::hash() const
{
    return std::hash<std::string>{}(hash_str());
}

} // ns mcl


template class mcl::CollisionPair<double>;
template class mcl::CollisionPair<float>;