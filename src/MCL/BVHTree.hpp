// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_CCD_BVHTREE_HPP
#define MCL_CCD_BVHTREE_HPP 1

#include "BVHLeaf.hpp"

#include <atomic>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Eigen // Forward declare BVH type, then include in cpp
{
template<typename T, int DIM, class Object>
class KdBVH;
}

namespace mcl {
namespace ccd {

/// @brief Traversal struct for things like point-in-element queries
template<typename T, int DIM, int PDIM = 3>
class BVHTraverse
{
  public:
    typedef BVHLeaf<T, DIM> ObjectType;
    typedef typename BVHLeaf<T, DIM>::BoxType VolumeType;
    virtual ~BVHTraverse() {}
    virtual bool intersectVolume(const VolumeType&) = 0; // true if move down branch
    virtual bool intersectObject(const ObjectType&) = 0; // true if stop traversing
};

/// @brief Wraps Eigen BVH and provides functions for CCD traversal
template<typename T, int DIM, int PDIM = 3>
class BVHTree
{
  public:
    typedef BVHLeaf<T, DIM> LeafType;
    typedef std::pair<Eigen::Vector4i, int> PairType;      // [sten, -1=invalid, 0=ee, 1=vf]
    typedef std::pair<int, bool> NodeIndex;                // [index, isleaf]
    static const size_t NumCandidates = DIM == 2 ? 6 : 15; // narrowphase candidates
    std::vector<LeafType> leaves;
    std::unique_ptr<Eigen::KdBVH<T, DIM, LeafType>> tree;

    struct Options
    {
        T vf_ccd_eta;    ///< gap for vf narrowphase
        T ee_ccd_eta;    ///< gap for ee narrowphase
        bool parallel;   ///< cpu-threaded traverse(...)
        bool discrete;   ///< discrete check at V1
        bool continuous; ///< ccd check from V0 to V1 (3D triangles only)
        Options()
            : vf_ccd_eta(1e-6)
            , ee_ccd_eta(1e-6)
            , parallel(true)
            , discrete(true)
            , continuous(true)
        {
        }
    } options;

    /// @brief Constructor
    BVHTree();

    /// @brief Destructor
    virtual ~BVHTree();

    /// @brief Updates the BVH.
    /// @param V0 nv x dim vertices at t=0
    /// @param V1 nv x dim vertices at t=1
    /// @param P np x pdim primitives (pdim: 3=triangles, 4=tets)
    /// @param active per-vertex boolean if the vertex should be checked, ignored if empty.
    template<typename DerivedV, typename DerivedP>
    inline void update(const Eigen::MatrixBase<DerivedV>& V0,
                       const Eigen::MatrixBase<DerivedV>& V1,
                       const Eigen::MatrixBase<DerivedP>& P,
                       const Eigen::VectorXi& active = Eigen::VectorXi());

    /// @brief Update BVH as above
    void update(const T* V0, const T* V1, const int* P, int np, const int* active = nullptr);

    /// @brief Traverses tree and checks from V0 to V1 with P = edges (DIM=2) or faces (DIM==3).
    /// Calls narrow_phase, append_pair, and append_discrete.
    void traverse(const T* V0, const T* V1, const int* P) const;

    /// @brief Traverse but with Eigen vectors (makes RowMajor copies)
    template<typename DerivedV, typename DerivedP>
    inline void traverse(const Eigen::MatrixBase<DerivedV>& V0,
                         const Eigen::MatrixBase<DerivedV>& V1,
                         const Eigen::MatrixBase<DerivedP>& P);

    /// @brief Traverse with iterator
    void traverse(BVHTraverse<T, DIM, PDIM>* traverser) const;

    /// @brief This function is called from a thread during CCD if two primitives collide
    /// (and options.continuous==true). It's how you retrieve continuous collisions.
    /// It is called from parallel threads (if options.parallel==true).
    /// If is_vf==false, it is an edge-edge collision.
    std::function<void(const Eigen::Vector4i& sten, bool is_vf, const T& toi)> append_pair;

    /// @brief This function is called from a thread if there is a discrete isect
    /// (and options.discrete==true). It's how you retrieve discrete collisions.
    /// It is called from parallel threads (if options.parallel==true).
    /// Return true to exit traversal immediately.
    std::function<bool(int p0, int p1)> append_discrete;

    /// @brief Optional: return true if the candidate pair should be skipped
    /// before narrow phase. If is_vf==false, it is an edge-edge query.
    std::function<bool(const Eigen::Vector4i& sten, bool is_vf)> filter_pair;

    /// @brief Optional: performs narrowphase to return time of impact
    /// (negative if no hit). If not set, defaults to ACCD kernels.
    /// If is_vf==false, it is an edge-edge query.
    std::function<T(const Eigen::Vector4i& sten, bool is_vf)> narrow_phase;

    /// @brief Optional: performs primitive-to-primitive collision. If not
    /// defined, uses default edge-edge, tri-tri, or point-in-tet (not tet-tet),
    /// based on PDIM. Returns true if primitives are intersecting.
    std::function<bool(int p0, int p1)> discrete_test;

  protected:
    /// @brief Helper function to get the children of the idx'th node.
    void get_children(const NodeIndex& idx, NodeIndex& l, NodeIndex& r) const;

    /// @brief Helper function to get the box of the idx'th node
    const typename LeafType::BoxType& get_box(const NodeIndex& idx) const;

    /// @brief Creates list of staring node pairs for traversal
    void make_frontlist(const NodeIndex& idx, std::vector<std::pair<NodeIndex, NodeIndex>>& queue) const;

    /// @brief Recursive traversal
    void collide(const T* V0,
                 const T* V1,
                 const int* P,
                 const NodeIndex& left,
                 const NodeIndex& right,
                 std::atomic<int>& stop) const;

    /// @brief Creates list of broadphase pairs from rep-tris and no shared vertex. pairs = stencil, type
    void get_candidates(int p0, int p1, const int* P, std::array<PairType, NumCandidates>& pairs) const;

    /// @brief Returns true if boxes intersect
    bool boxes_intersect(const NodeIndex& left, const NodeIndex& right) const;

    /// @brief If narrow_phase function pointer is not defined, this function is used (ACCD)
    T default_narrow_phase(const T* V0, const T* V1, const Eigen::Vector4i& sten, bool is_vf) const;

    /// @brief If discrete_test function pointer is not defined, this function is used.
    bool default_discrete_test(const T* V, const int* p0, const int* p1) const;

}; // class BVHTree

// Convert matrix to needed type
template<typename T, int DIM, int PDIM>
template<typename DerivedV, typename DerivedP>
inline void
BVHTree<T, DIM, PDIM>::update(const Eigen::MatrixBase<DerivedV>& V0,
                              const Eigen::MatrixBase<DerivedV>& V1,
                              const Eigen::MatrixBase<DerivedP>& P,
                              const Eigen::VectorXi& active)
{
    // TODO update this to avoid copy when Derived is row major.
    if (V0.rows() != V1.rows()) {
        return;
    }
    if (V0.cols() != DIM || V1.cols() != DIM) {
        throw std::runtime_error("V.cols != DIM");
    }
    if (P.cols() != PDIM) {
        throw std::runtime_error("P.cols != PDIM");
    }
    const int* active_ptr = active.size() == V0.rows() ? active.data() : nullptr;
    typedef Eigen::Matrix<T, Eigen::Dynamic, DIM, Eigen::RowMajor> ScalarMat;
    typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> IndexMat;
    ScalarMat V0cpy = V0.template cast<T>();
    ScalarMat V1cpy = V1.template cast<T>();
    IndexMat Pcpy = P.template cast<int>();
    update(V0cpy.data(), V1cpy.data(), Pcpy.data(), P.rows(), active_ptr);
}

// Convert matrix to needed type
template<typename T, int DIM, int PDIM>
template<typename DerivedV, typename DerivedP>
inline void
BVHTree<T, DIM, PDIM>::traverse(const Eigen::MatrixBase<DerivedV>& V0,
                                const Eigen::MatrixBase<DerivedV>& V1,
                                const Eigen::MatrixBase<DerivedP>& P)
{
    // TODO update this to avoid copy when Derived is row major.
    if (V0.rows() != V1.rows()) {
        return;
    }
    if (V0.cols() != DIM || V1.cols() != DIM) {
        throw std::runtime_error("V.cols != DIM");
    }
    if (P.cols() != PDIM) {
        throw std::runtime_error("P.cols != PDIM");
    }
    typedef Eigen::Matrix<T, Eigen::Dynamic, DIM, Eigen::RowMajor> ScalarMat;
    typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> IndexMat;
    ScalarMat V0cpy = V0.template cast<T>();
    ScalarMat V1cpy = V1.template cast<T>();
    IndexMat Pcpy = P.template cast<int>();
    traverse(V0cpy.data(), V1cpy.data(), Pcpy.data());
}

} // end namespace ccd
} // end namespace mcl

#endif
