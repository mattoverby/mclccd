// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_BVHTREE_HPP
#define MCL_BVHTREE_HPP 1

#include "BVHLeaf.hpp"
#include <atomic>
#include <memory>

namespace Eigen // Forward declare BVH type, then include in cpp
{ template <typename T,int DIM, class Object> class KdBVH; }

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

// Traversal struct for things like point-in-elem queries
template <typename T, int DIM>
class BVHTraverse
{
public:
    typedef BVHLeaf<T,DIM> ObjectType;
    typedef typename BVHLeaf<T,DIM>::BoxType VolumeType;
    virtual ~BVHTraverse() {}
    virtual bool intersectVolume(const VolumeType&) = 0; // true if move down branch
    virtual bool intersectObject(const ObjectType&) = 0; // true if stop traversing
};

template <typename T, int DIM>
class BVHTree
{
public:
    typedef BVHLeaf<T,DIM> LeafType;
    typedef Eigen::Matrix<T,DIM,1> VecType;
    typedef std::pair<Eigen::Vector4i,int> PairType;
    typedef std::pair<int, bool> NodeIndex; // [index, isleaf]

    std::vector<LeafType> leaves;
    std::shared_ptr<Eigen::KdBVH<T,DIM,LeafType>> tree;

    struct Options
    {
        T box_eta; // padding of tree boxes
        T ccd_eta; // gap for narrowphase <= box_eta
        bool reptri; // representative triangles
        bool threaded; // cpu-threaded traverse(...)
        bool discrete; // discrete check ff/ee at t=1
        bool continuous; // ccd check from t=0 to t=1
        bool vf_one_sided; // allow pass-through if norm dir
        Options() :
            box_eta(1e-6),
            ccd_eta(1e-6),
            reptri(true),
            threaded(true),
            discrete(true),
            continuous(true),
            vf_one_sided(false)
            {}
    } options;

    BVHTree();
    virtual ~BVHTree() {}

    // V is n x dim with V0 at t=0 and V1 at t=1.
    // P can be m x (whatever) depending on the prim (1=pt, 2=edge, 3=tri, etc)
    // Active is a per-vertex boolean if the vertex should be checked.
    // if active buffer empty, all vertices considered active.
    template <typename DerivedV, typename DerivedP>
    inline void update(
        const Eigen::MatrixBase<DerivedV> &V0,
        const Eigen::MatrixBase<DerivedV> &V1,
        const Eigen::MatrixBase<DerivedP> &P,
        const Eigen::VectorXi &active = Eigen::VectorXi());

    // Update BVH using raw data.
    void update(const T* V0, const T* V1, const int *P, int np, int pdim,
        const int *active = nullptr);

    // Traverses tree and checks from V0 to V1 with P = edges (DIM=2) or faces (DIM==3).
    // Calls narrow_phase, append_pair, and append_discrete.
    template <typename DerivedV, typename DerivedP>
    inline void traverse(
        const Eigen::MatrixBase<DerivedV> &V0,
        const Eigen::MatrixBase<DerivedV> &V1,
        const Eigen::MatrixBase<DerivedP> &P);
    
    // Traverse with raw data.
    void traverse(const T* V0, const T* V1, const int *P) const;
    
    // Traverse with iterator
    void traverse(BVHTraverse<T,DIM> *traverser) const;

    // This function is called from a thread during CCD if two primitives collide
    // (and options.continuous==true).
    std::function<void(const Eigen::Vector4i &sten, int type, const T& toi)> append_pair;

    // Return true if the candidate pair should be skipped before narrow phase.
    std::function<bool(const Eigen::Vector4i &sten, int type)> filter_pair;

    // This function is called from a thread if there is a discrete isect
    // (and options.discrete==true). Returns true to exit traveral immediately.
    std::function<bool(int p0, int p1)> append_discrete;

    // Peforms narrowphase to return toi (negative if no hit). If not set,
    // defaults to CTCD kernels.
    std::function<T(const Eigen::Vector4i &sten, int type)> narrow_phase;

protected:

    void get_children(const NodeIndex &idx, NodeIndex &l, NodeIndex &r) const;

    const typename LeafType::BoxType &get_box(const NodeIndex &idx) const;

    // Creates list of staring node pairs for traversal
    void make_frontlist(const NodeIndex &idx,
        std::vector<std::pair<NodeIndex,NodeIndex>> &queue) const;

    // Recursive traversal
    void collide(const T* V0, const T* V1, const int *P,
        const NodeIndex &left, const NodeIndex &right,
        std::atomic<int> &stop) const;

    // Creates list of broadphase pairs from rep-tris and no shared vertex
    // pairs = stencil, type
    void get_candidates(int p0, int p1, const int *P,
        std::vector<PairType> &pairs) const;

    bool boxes_intersect(const NodeIndex &left, const NodeIndex &right) const;

    T default_narrow_phase(const T* V0, const T* V1,
        const Eigen::Vector4i &sten, int type) const;

    bool default_discrete_test(const T* V, const int *p0, const int *p1) const;

}; // class BVHTree

// Convert matrix to needed type
template <typename T, int DIM>
template <typename DerivedV, typename DerivedP>
inline void BVHTree<T,DIM>::update(
    const Eigen::MatrixBase<DerivedV> &V0,
    const Eigen::MatrixBase<DerivedV> &V1,
    const Eigen::MatrixBase<DerivedP> &P,
    const Eigen::VectorXi &active)
{
    // TODO avoid copy. But MatrixBase doesn't have a data() function.
    if (V0.rows() != V1.rows()) { return; }
    const int *active_ptr = active.size() == V0.rows() ? active.data() : nullptr;
    typedef Eigen::Matrix<T,Eigen::Dynamic,DIM,Eigen::RowMajor> ScalarMat;
    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> IndexMat;
    ScalarMat V0cpy = V0.template cast<T>();
    ScalarMat V1cpy = V1.template cast<T>();
    IndexMat Pcpy = P.template cast<int>();
    update(V0cpy.data(), V1cpy.data(), Pcpy.data(), P.rows(), P.cols(), active_ptr);
}

// Convert matrix to needed type
template <typename T, int DIM>
template <typename DerivedV, typename DerivedP>
inline void BVHTree<T,DIM>::traverse(
    const Eigen::MatrixBase<DerivedV> &V0,
    const Eigen::MatrixBase<DerivedV> &V1,
    const Eigen::MatrixBase<DerivedP> &P)
{
    // TODO avoid copy. But MatrixBase doesn't have a data() function.
    if (V0.rows() != V1.rows()) { return; }
    typedef Eigen::Matrix<T,Eigen::Dynamic,DIM,Eigen::RowMajor> ScalarMat;
    typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> IndexMat;
    ScalarMat V0cpy = V0.template cast<T>();
    ScalarMat V1cpy = V1.template cast<T>();
    IndexMat Pcpy = P.template cast<int>();
    traverse(V0cpy.data(), V1cpy.data(), Pcpy.data());
}

} // end ns mcl

#endif
