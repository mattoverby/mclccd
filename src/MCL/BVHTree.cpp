// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "BVHTree.hpp"
#include "NarrowPhase.hpp"

// We use our own edited version of KdBVH
// so it must be included before Eigen/BVH.
#include "KdBVH.hpp"
#include <unsupported/Eigen/BVH>

#include <chrono>
#include <string>
#include <tbb/parallel_for.h>

namespace mcl {
namespace ccd {

class Timer
{
  public:
    std::chrono::steady_clock::time_point start_time;
    Timer()
        : start_time(std::chrono::steady_clock::now())
    {
    }
    double elapsed_ms() const
    {
        auto now = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time);
        return double(duration.count());
    }
};

/// @brief Returns a primtive
template<int PDIM>
void
get_primitive(int prim_index, const int* P, int prim[PDIM])
{
    prim[0] = P[prim_index * PDIM + 0];
    if constexpr (PDIM > 1) {
        prim[1] = P[prim_index * PDIM + 1];
    }
    if constexpr (PDIM > 2) {
        prim[2] = P[prim_index * PDIM + 2];
    }
    if constexpr (PDIM > 3) {
        prim[3] = P[prim_index * PDIM + 3];
    }
}

/// @brief Returns a vertex
template<typename T, typename DerivedV>
void
get_vertex(int vert_index, const T* V, Eigen::MatrixBase<DerivedV>& vert)
{
    constexpr int DIM = Eigen::MatrixBase<DerivedV>::SizeAtCompileTime;
    if (vert_index < 0) {
        return;
    }
    vert[0] = V[vert_index * DIM + 0];
    vert[1] = V[vert_index * DIM + 1];
    if constexpr (DIM > 2) {
        vert[2] = V[vert_index * DIM + 2];
    }
}

/// @brief Returns all vertices of a primitive (up to PDIM=4)
template<typename T, int DIM, int PDIM>
std::array<Eigen::Vector<T, DIM>, PDIM>
get_vertices(const T* V, const int* prim)
{
    std::array<Eigen::Vector<T, DIM>, PDIM> verts;
    get_vertex(prim[0], V, verts[0]);
    if constexpr (PDIM > 1) {
        get_vertex(prim[1], V, verts[1]);
    }
    if constexpr (PDIM > 2) {
        get_vertex(prim[2], V, verts[2]);
    }
    if constexpr (PDIM > 3) {
        get_vertex(prim[3], V, verts[3]);
    }
    return verts;
}

/// @brief Returns true if two primitives share a vertex (up to PDIM=4)
template<int PDIM>
bool
prims_share_vertex(const int* p0, const int* p1)
{
    if (p0[0] == p1[0]) {
        return true;
    }
    if constexpr (PDIM > 1) {
        if (p0[0] == p1[1] || p0[1] == p1[1] || p0[1] == p1[0]) {
            return true;
        }
    }
    if constexpr (PDIM > 2) {
        if (p0[0] == p1[2] || p0[1] == p1[2] || p0[2] == p1[0] || p0[2] == p1[1] || p0[2] == p1[2]) {
            return true;
        }
    }
    if constexpr (PDIM > 3) {
        if (p0[0] == p1[3] || p0[1] == p1[3] || p0[2] == p1[3] || //
            p0[3] == p1[3] || p0[3] == p1[0] || p0[3] == p1[1] || //
            p0[3] == p1[2]) {
            return true;
        }
    }
    return false;
}

/// @brief Returns true if the vertex index is in the primitive (up to PDIM=4)
template<int PDIM>
bool
vertex_prim_share_vertex(int vi, const int* p)
{
    if (vi == p[0]) {
        return true;
    }
    if constexpr (PDIM > 1) {
        if (vi == p[1]) {
            return true;
        }
    }
    if constexpr (PDIM > 2) {
        if (vi == p[2]) {
            return true;
        }
    }
    if constexpr (PDIM > 3) {
        if (vi == p[3]) {
            return true;
        }
    }
    return false;
}

template<typename T, int DIM, int PDIM>
BVHTree<T, DIM, PDIM>::BVHTree()
{
    tree = std::make_unique<Eigen::KdBVH<T, DIM, LeafType>>();
}

template<typename T, int DIM, int PDIM>
BVHTree<T, DIM, PDIM>::~BVHTree() = default;

template<typename T, int DIM, int PDIM>
void
BVHTree<T, DIM, PDIM>::update(const T* V0, const T* V1, const int* P, int np, const int* active)
{
    if (np == 0) {
        leaves.clear();
        tree = std::make_unique<Eigen::KdBVH<T, DIM, LeafType>>();
        return;
    }

    assert(options.vf_ccd_eta >= 0);
    assert(options.ee_ccd_eta >= 0);
    T box_eta = options.vf_ccd_eta;
    if (DIM == 3) {
        box_eta = std::max(box_eta, options.ee_ccd_eta);
    }

    bool update_reptri = false;
    if ((int)leaves.size() != np) {
        update_reptri = true;
        leaves.resize(np);
    }

    // Update representative triangles
    if (update_reptri) {

        int n_verts_guess = np; // reallocate as needed
        std::vector<int> seen_verts(n_verts_guess, 0);
        std::vector<std::vector<int>> seen_edges(n_verts_guess);

        for (int i = 0; i < np; ++i) {
            BVHLeaf<T, DIM>& leaf = leaves[i];
            leaf.v.setZero();
            leaf.e.setZero();
            int prim[PDIM];
            get_primitive<PDIM>(i, P, prim);
            int maxInd = *std::max_element(prim, prim + PDIM);
            while (maxInd >= n_verts_guess) {
                n_verts_guess *= 2;
                seen_verts.resize(n_verts_guess, 0);
                seen_edges.resize(n_verts_guess);
            }

            // TODO: Update for PDIM != 3
            for (int j = 0; j < 3; ++j) {
                int vi = prim[j];
                if (seen_verts[vi] == 0) {
                    seen_verts[vi] = 1;
                    leaf.v[j] = 1;
                }

                int e0 = vi;
                int e1 = prim[(j + 1) % PDIM];
                if (e1 < e0) {
                    std::swap(e0, e1);
                }

                std::vector<int>& neighbors = seen_edges[e0];
                if (std::find(neighbors.begin(), neighbors.end(), e1) == neighbors.end()) {
                    neighbors.emplace_back(e1);
                    leaf.e[j] = 1;
                }
            }
        }
    }

    // Update bounding boxes
    tbb::parallel_for(tbb::blocked_range<int>(0, np), [&](const tbb::blocked_range<int>& r) {
        for (int i = r.begin(); i < r.end(); ++i) {
            BVHLeaf<T, DIM>& leaf = leaves[i];

            // Is active if all active, or any vertex in p active
            leaf.box.active = active == nullptr;
            leaf.idx = i;
            leaf.box.setEmpty();
            leaf.box.t0.setEmpty();
            leaf.box.t1.setEmpty();
            for (int j = 0; j < PDIM; ++j) {
                int vi = P[i * PDIM + j];
                Eigen::Vector<T, DIM> xi_t0 = Eigen::Vector<T, DIM>::Zero();
                Eigen::Vector<T, DIM> xi_t1 = Eigen::Vector<T, DIM>::Zero();
                get_vertex(vi, V0, xi_t0);
                get_vertex(vi, V1, xi_t1);
                leaf.box.t0.extend(xi_t0);
                leaf.box.t1.extend(xi_t1);
                leaf.box.extend(xi_t0);
                leaf.box.extend(xi_t1);
                if (active != nullptr) {
                    if (active[vi]) {
                        leaf.box.active = true;
                    }
                }
            }
            leaf.box.t0.min().array() -= box_eta;
            leaf.box.t1.min().array() -= box_eta;
            leaf.box.min().array() -= box_eta;
            leaf.box.t0.max().array() += box_eta;
            leaf.box.t1.max().array() += box_eta;
            leaf.box.max().array() += box_eta;
        }
    });

    // Initialize BVH
    tree->init(leaves.begin(), leaves.end());
}

template<typename T, int DIM, int PDIM>
void
BVHTree<T, DIM, PDIM>::traverse(const T* V0, const T* V1, const int* P) const
{
    if (leaves.size() <= 1) {
        return;
    }

    // Create a list of initial BVH intersection lists
    std::vector<std::pair<NodeIndex, NodeIndex>> queue;
    queue.reserve(leaves.size());
    make_frontlist(std::make_pair(tree->getRootIndex(), false), queue);
    if (queue.empty()) {
        return;
    }

    std::atomic<int> stop(0);
    int nq = queue.size();
    if (options.parallel) {
        tbb::parallel_for(tbb::blocked_range<int>(0, nq), [&](const tbb::blocked_range<int>& r) {
            for (int i = r.begin(); i < r.end(); ++i) {
                collide(V0, V1, P, queue[i].first, queue[i].second, stop);
            }
        });
    } else {
        for (int i = 0; i < nq; ++i) {
            collide(V0, V1, P, queue[i].first, queue[i].second, stop);
        }
    }
}

template<typename T, int DIM, int PDIM>
void
BVHTree<T, DIM, PDIM>::traverse(BVHTraverse<T, DIM, PDIM>* traverser) const
{
    const Eigen::KdBVH<T, DIM, LeafType>& tree_ref = *tree.get();
    Eigen::BVIntersect(tree_ref, *traverser);
}

template<typename T, int DIM, int PDIM>
void
BVHTree<T, DIM, PDIM>::make_frontlist(const NodeIndex& idx, std::vector<std::pair<NodeIndex, NodeIndex>>& queue) const
{
    if (idx.second) { // is leaf node
        return;
    }

    NodeIndex l, r;
    get_children(idx, l, r);

    if (l.first >= 0) {
        make_frontlist(l, queue);
    }

    if (r.first >= 0) {
        make_frontlist(r, queue);
    }

    // if both leaf, should not have same idx
    assert((l.second == r.second) ? (l.first != r.first) : true);

    queue.emplace_back(l, r);
}

template<typename T, int DIM, int PDIM>
void
BVHTree<T, DIM, PDIM>::get_children(const NodeIndex& idx, NodeIndex& l, NodeIndex& r) const
{
    typedef typename Eigen::KdBVH<T, DIM, BVHLeaf<T, DIM>>::VolumeIterator VolumeIter;
    typedef typename Eigen::KdBVH<T, DIM, BVHLeaf<T, DIM>>::ObjectIterator ObjectIter;
    l.first = -1;
    l.second = false;
    r.first = -1;
    r.second = false;

    if (idx.second) // is a leaf node
    {
        return;
    }

    VolumeIter vBegin = nullptr, vEnd = nullptr;
    ObjectIter oBegin = nullptr, oEnd = nullptr;
    tree->getChildren(idx.first, vBegin, vEnd, oBegin, oEnd);
    // int num_children = 0;
    bool is_l = true;
    for (; vBegin != vEnd; ++vBegin) {
        if (is_l) {
            l = NodeIndex(*vBegin, false);
            is_l = false;
        } else {
            r = NodeIndex(*vBegin, false);
        }
        // num_children++;
    }
    for (; oBegin != oEnd; ++oBegin) {
        if (is_l) {
            l = NodeIndex(oBegin->idx, true);
            is_l = false;
        } else {
            r = NodeIndex(oBegin->idx, true);
        }
        // num_children++;
    }
    // assert(num_children == 2); // assumes binary tree
}

template<typename T, int DIM, int PDIM>
const typename BVHLeaf<T, DIM>::BoxType&
BVHTree<T, DIM, PDIM>::get_box(const NodeIndex& idx) const
{
    if (idx.second) {
        return leaves[idx.first].box;
    }
    return tree->getVolume(idx.first);
}

template<typename T, int DIM, int PDIM>
void
BVHTree<T, DIM, PDIM>::collide(const T* V0,
                               const T* V1,
                               const int* P,
                               const NodeIndex& left,
                               const NodeIndex& right,
                               std::atomic<int>& stop) const
{
    // Check if we should stop traversing, only possible if discrete test
    if (options.discrete && append_discrete != nullptr) {
        if (stop.load() > 0) {
            return;
        }
    }

    // if both leaf, not eq
    // assert(left.first >= 0 && right.first >= 0);
    // assert((left.second == right.second) ? (left.first != right.first) : true);

    // Check if nodes intersect or inactive
    if (!boxes_intersect(left, right)) {
        return;
    }

    // Both are leaf nodes
    if (left.second && right.second) {
        if (options.continuous) {
            std::array<PairType, NumCandidates> pairs;
            get_candidates(left.first, right.first, P, pairs);
            for (size_t i = 0; i < NumCandidates; ++i) {
                if (pairs[i].second < 0) {
                    continue;
                }

                T toi = -1;
                if (narrow_phase != nullptr) {
                    toi = narrow_phase(pairs[i].first, pairs[i].second);
                } else {
                    toi = default_narrow_phase(V0, V1, pairs[i].first, pairs[i].second);
                }
                if (toi >= 0 && append_pair != nullptr) {
                    append_pair(pairs[i].first, pairs[i].second, toi);
                }
            }
        }

        if (options.discrete) {
            bool prims_intersected = false;
            if (discrete_test != nullptr) {
                prims_intersected = discrete_test(left.first, right.first);
            } else {
                int p0[PDIM], p1[PDIM];
                get_primitive<PDIM>(left.first, P, p0);
                get_primitive<PDIM>(right.first, P, p1);
                prims_intersected = default_discrete_test(V1, p0, p1);
            }
            if (prims_intersected && append_discrete != nullptr) {
                bool stop_traverse = append_discrete(left.first, right.first);
                if (stop_traverse) {
                    stop++;
                }
            }
        }

        return;
    }

    if (left.second) // left is leaf
    {
        NodeIndex l, r;
        get_children(right, l, r);
        collide(V0, V1, P, left, l, stop);
        collide(V0, V1, P, left, r, stop);
    } else {
        NodeIndex l, r;
        get_children(left, l, r);
        collide(V0, V1, P, l, right, stop);
        collide(V0, V1, P, r, right, stop);
    }
}

template<typename T, int DIM, int PDIM>
void
BVHTree<T, DIM, PDIM>::get_candidates(int p0, int p1, const int* P, std::array<PairType, NumCandidates>& pairs) const
{
    if constexpr (PDIM != 3)
    {
        printf("TODO: continuous collision for PDIM != 3 (rep tris)");
        return;
    }

    const LeafType& l0 = leaves[p0];
    const LeafType& l1 = leaves[p1];

    for (size_t i = 0; i < NumCandidates; ++i) {
        pairs[i].first.setZero();
        pairs[i].second = -1;
    }

    int f0[PDIM], f1[PDIM];
    get_primitive<PDIM>(p0, P, f0);
    get_primitive<PDIM>(p1, P, f1);

    int pair_index = 0;

    // VF, f0 -> f1
    for (int i = 0; i < PDIM; ++i) {
        if (!l0.v[i]) // v not represented
        {
            pairs[pair_index++].second = -1;
            continue;
        }

        if (vertex_prim_share_vertex<PDIM>(f0[i], f1)) {
            pairs[pair_index++].second = -1;
            continue;
        }

        PairType& pair = pairs[pair_index++];
        pair.second = 1; // is_vf

        pair.first = Eigen::Vector4i(-1, -1, -1, -1);
        pair.first[0] = f0[i];
        for (int j = 0; j < PDIM; ++j)
            pair.first[j + 1] = f1[j];

        if (filter_pair != nullptr)
            if (filter_pair(pair.first, pair.second))
                pair.second = -1;
    }

    // VF, f1 -> f0
    for (int i = 0; i < PDIM; ++i) {
        if (!l1.v[i]) // v not represented
        {
            pairs[pair_index++].second = -1;
            continue;
        }

        if (vertex_prim_share_vertex<PDIM>(f1[i], f0)) {
            pairs[pair_index++].second = -1;
            continue;
        }

        PairType& pair = pairs[pair_index++];
        pair.second = 1; // is_vf

        pair.first = Eigen::Vector4i(-1, -1, -1, -1);
        pair.first[0] = f1[i];
        for (int j = 0; j < PDIM; ++j)
            pair.first[j + 1] = f0[j];

        if (filter_pair != nullptr) {
            if (filter_pair(pair.first, pair.second)) {
                pair.second = -1;
            }
        }
    }

    // Skip edge-edge if not in 3D
    if (DIM != 3) {
        return;
    }

    // EE
    for (int i = 0; i < 3; ++i) {
        if (!l0.e[i]) {
            continue;
        }
        for (int j = 0; j < 3; ++j) {
            if (!l1.e[j]) {
                pairs[pair_index++].second = -1;
                continue;
            }

            Eigen::Vector4i sten(i, (i + 1) % 3, j, (j + 1) % 3);

            // shares vertex?
            if (f0[sten[0]] == f1[sten[2]] || f0[sten[0]] == f1[sten[3]] || f0[sten[1]] == f1[sten[2]] ||
                f0[sten[1]] == f1[sten[3]]) {
                pairs[pair_index++].second = -1;
                continue;
            }

            PairType& pair = pairs[pair_index++];
            pair.second = 0; // is not vf
            pair.first = Eigen::Vector4i(f0[sten[0]], f0[sten[1]], f1[sten[2]], f1[sten[3]]);

            if (filter_pair != nullptr)
                if (filter_pair(pair.first, pair.second))
                    pair.second = -1;
        }
    }
}

template<typename T, int DIM, int PDIM>
bool
BVHTree<T, DIM, PDIM>::boxes_intersect(const NodeIndex& left, const NodeIndex& right) const
{
    const typename BVHLeaf<T, DIM>::BoxType& lbox = get_box(left);
    const typename BVHLeaf<T, DIM>::BoxType& rbox = get_box(right);

    // Check if both branches inactive
    if (!lbox.active && !rbox.active) {
        return false;
    }

    // Check nodes intersect
    if (!lbox.intersects(rbox)) {
        return false;
    }

    // Check time varying boxes if not empty
    if (lbox.t0.isEmpty() || rbox.t0.isEmpty()) {
        return true;
    }
    if (lbox.t1.isEmpty() || rbox.t1.isEmpty()) {
        return true;
    }

    for (int i = 0; i < DIM; ++i) {
        if ((lbox.t0.min()[i] > rbox.t0.max()[i] && lbox.t1.min()[i] > rbox.t1.max()[i]) ||
            (lbox.t0.max()[i] < rbox.t0.min()[i] && lbox.t1.max()[i] < rbox.t1.min()[i]))
            return false;
    }

    return true;
}

template<typename T, int DIM, int PDIM>
T
BVHTree<T, DIM, PDIM>::default_narrow_phase(const T* V0, const T* V1, const Eigen::Vector4i& sten, bool is_vf) const
{
    auto verts_t0 = get_vertices<T, DIM, 4>(V0, sten.data());
    auto verts_t1 = get_vertices<T, DIM, 4>(V1, sten.data());
    T toi = -1;
    int hit = 0;
    if (is_vf) {
        hit = NarrowPhaseACCD<T, DIM>::query_ccd_vf(verts_t0.data(), verts_t1.data(), options.vf_ccd_eta, toi);
    } else {
        hit = NarrowPhaseACCD<T, DIM>::query_ccd_ee(verts_t0.data(), verts_t1.data(), options.ee_ccd_eta, toi);
    }

    if (hit == 1) {
        return toi;
    }

    return -1;
}

template<typename T, int DIM, int PDIM>
bool
BVHTree<T, DIM, PDIM>::default_discrete_test(const T* V, const int* p0, const int* p1) const
{
    if (prims_share_vertex<PDIM>(p0, p1)) {
        return false;
    }

    auto p = get_vertices<T, DIM, PDIM>(V, p0);
    auto q = get_vertices<T, DIM, PDIM>(V, p1);

    if constexpr (DIM == 3) {
        if constexpr (PDIM == 3) { // 3D triangles
            return NarrowPhase<T>::discrete_tri_tri(p[0], p[1], p[2], q[0], q[1], q[2]);
        }
        if constexpr (PDIM == 4) { // 3D tets
            // NOTE: Point-in-tet is not a full tet-tet intersection test.
            // But, point-in-tet is usually what I want if doing tet collisions :)
            for (int i = 0; i < 4; ++i) {
                if (NarrowPhase<T>::point_in_tet(p[i], q[0], q[1], q[2], q[3])) {
                    return true;
                }
                if (NarrowPhase<T>::point_in_tet(q[i], p[0], p[1], p[2], p[3])) {
                    return true;
                }
            }
        }
    } else if constexpr (DIM == 2) {
        if constexpr (PDIM == 2) { // 2D edges
            return NarrowPhase<T>::discrete_edge_edge(p[0], p[1], q[0], q[1]);
        } else if constexpr (PDIM == 3) { // 2D triangles
            Eigen::Vector2<T> v0[2];
            Eigen::Vector2<T> v1[2];
            for (int i = 0; i < 3; ++i) {
                get_vertex(p0[i], V, v0[0]);
                get_vertex(p0[(i + 1) % 3], V, v0[1]);
                for (int j = 0; j < 3; ++j) {
                    get_vertex(p1[j], V, v1[0]);
                    get_vertex(p1[(j + 1) % 3], V, v1[1]);
                    if (NarrowPhase<T>::discrete_edge_edge(v0[0], v0[1], v1[0], v1[1])) {
                        return true;
                    }
                }
            }
            return false;
        }
    }

    return false;
}

} // end namespace ccd
} // end namespace mcl

template class mcl::ccd::BVHTree<double, 3, 3>;
template class mcl::ccd::BVHTree<double, 2, 3>;
template class mcl::ccd::BVHTree<double, 3, 4>;
template class mcl::ccd::BVHTree<float, 3, 3>;
template class mcl::ccd::BVHTree<float, 2, 3>;
template class mcl::ccd::BVHTree<float, 3, 4>;
