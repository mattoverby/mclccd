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

// Wraps Eigen BVH and provides functions for CCD traversal
template <typename T, int DIM>
class BVHTree
{
public:
    typedef BVHLeaf<T,DIM> LeafType;
    typedef Eigen::Matrix<T,DIM,1> VecType;
    typedef std::pair<Eigen::Vector4i,bool> PairType; // [sten, is_vf]
    typedef std::pair<int, bool> NodeIndex; // [index, isleaf]

    std::vector<LeafType> leaves;
    std::shared_ptr<Eigen::KdBVH<T,DIM,LeafType> > tree;

    struct Options
    {
        T box_eta; // padding of tree boxes
        T vf_ccd_eta; // gap for vf narrowphase <= box_eta
        T ee_ccd_eta; // gap for ee narrowphase <= box_eta
        bool reptri; // representative triangles
        bool threaded; // cpu-threaded traverse(...)
        bool discrete; // discrete check ff/ee at V1
        bool continuous; // ccd check from V0 to V1
        bool vf_one_sided; // allow pass-through if norm dir
        bool ee_robust; // check VV and VE
        Options() :
            box_eta(1e-6),
            vf_ccd_eta(1e-6),
            ee_ccd_eta(1e-6),
            reptri(true),
            threaded(true),
            discrete(true),
            continuous(true),
            vf_one_sided(false),
            ee_robust(true)
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
    void traverse(const T* V0, const T* V1, const int *P) const;

    // Traverse but with Eigen vectors (makes RowMajor copies)
    template <typename DerivedV, typename DerivedP>
    inline void traverse(
        const Eigen::MatrixBase<DerivedV> &V0,
        const Eigen::MatrixBase<DerivedV> &V1,
        const Eigen::MatrixBase<DerivedP> &P);
    
    // Traverse with iterator
    void traverse(BVHTraverse<T,DIM> *traverser) const;

    // This function is called from a thread during CCD if two primitives collide
    // (and options.continuous==true). It's how you retrieve continuous collisions.
    // It is called from parallel threads (if options.threaded==true).
    // If is_vf==false, it is an edge-edge collision.
    std::function<void(const Eigen::Vector4i &sten, bool is_vf, const T& toi)> append_pair;

    // This function is called from a thread if there is a discrete isect
    // (and options.discrete==true). It's how you retrieve discrete collisions.
    // It is called from parallel threads (if options.threaded==true).
    // Return true to exit traversal immediately.
    std::function<bool(int p0, int p1)> append_discrete;

    // Optional:
    // Return true if the candidate pair should be skipped before narrow phase.
    // If is_vf==false, it is an edge-edge query.
    std::function<bool(const Eigen::Vector4i &sten, bool is_vf)> filter_pair;

    // Optional:
    // Performs narrowphase to return time of impact (negative if no hit).
    // If not set, defaults to CTCD kernels.
    // If is_vf==false, it is an edge-edge query.
    std::function<T(const Eigen::Vector4i &sten, bool is_vf)> narrow_phase;

protected:

    void get_children(const NodeIndex &idx, NodeIndex &l, NodeIndex &r) const;

    const typename LeafType::BoxType &get_box(const NodeIndex &idx) const;

    // Creates list of staring node pairs for traversal
    void make_frontlist(const NodeIndex &idx,
        std::vector<std::pair<NodeIndex,NodeIndex> > &queue) const;

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
        const Eigen::Vector4i &sten, bool is_vf) const;

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
