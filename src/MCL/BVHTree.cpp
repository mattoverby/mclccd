// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "BVHTree.hpp"
#include "NarrowPhase.hpp"

// We use our own edited version of KdBVH
// so it must be included before Eigen/BVH.
#include "KdBVH.hpp"
#include <unsupported/Eigen/BVH>

#include <tbb/parallel_for.h>
#include <unordered_set>
#include <set>
#include <string>

namespace mcl
{
namespace ccd
{

template <typename T, int DIM, int PDIM>
BVHTree<T,DIM,PDIM>::BVHTree()
{
    tree = std::make_shared<Eigen::KdBVH<T,DIM,LeafType> >();
}

template <typename T, int DIM, int PDIM>
void BVHTree<T,DIM,PDIM>::update(const T* V0, const T* V1, const int *P, int np, const int *active)
{
    if (np == 0)
    {
        leaves.clear();
        tree = std::make_shared<Eigen::KdBVH<T,DIM,LeafType> >();
        return;
    }

    assert(options.vf_ccd_eta >= 0);
    assert(options.ee_ccd_eta >= 0);
    T box_eta = options.vf_ccd_eta;
    if (DIM == 3)
    {
        box_eta = std::max(box_eta, options.ee_ccd_eta);
    }

    bool update_reptri = false;

    if ((int)leaves.size() != np)
    {
        leaves.resize(np);
        update_reptri = PDIM == 2 || PDIM == 3;
    }

    // Update representative triangles
    if (update_reptri)
    {
        std::unordered_set<int> seen_verts;
        std::set<std::string> seen_edges;

        for(int i=0; i<np; ++i)
        {
            BVHLeaf<T,DIM> &leaf = leaves[i];
            leaf.v.setZero();
            leaf.e.setZero();
            for (int j=0; j<PDIM; ++j)
            {
                int vi = P[i*PDIM+j];
                bool v_not_seen = seen_verts.emplace(vi).second;
                if (v_not_seen)
                    leaf.v[j]=1;

                if (PDIM != 3)
                    continue;

                int e0 = vi;
                int e1 = P[i*PDIM+((j+1)%3)];
                if (e1 < e0) { std::swap(e0, e1); }
                std::string h = std::to_string(e0)+' '+std::to_string(e1);
                bool e_not_seen = seen_edges.emplace(h).second;
                if (e_not_seen)
                    leaf.e[j]=1;
            }
        }
    }

    // Update bounding boxes
    tbb::parallel_for(tbb::blocked_range<int>(0,np), [&](const tbb::blocked_range<int>& r)
    {
      for (int i=r.begin(); i<r.end(); ++i)
      {
        BVHLeaf<T,DIM> &leaf = leaves[i];

        // Is active if all active, or any vertex in p active
        leaf.box.active = active==nullptr;
        leaf.idx = i;
        leaf.box.setEmpty();
        leaf.box.t0.setEmpty();
        leaf.box.t1.setEmpty();
        for (int j=0; j<PDIM; ++j)
        {
            int vi = P[i*PDIM + j];
            VecType xi_t0 = VecType::Zero();
            VecType xi_t1 = VecType::Zero();
            for (int k=0; k<DIM; ++k)
            {
                xi_t0[k] = V0[vi*DIM+k];
                xi_t1[k] = V1[vi*DIM+k];
            }
            leaf.box.t0.extend(xi_t0);
            leaf.box.t1.extend(xi_t1);
            leaf.box.extend(xi_t0);
            leaf.box.extend(xi_t1);

            if (active != nullptr)
                if (active[vi])
                    leaf.box.active = true;
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

template <typename T, int DIM, int PDIM>
void BVHTree<T,DIM,PDIM>::traverse(const T* V0, const T* V1, const int *P) const
{
    if (leaves.size() <= 1)
        return;

    // Create a list of initial BVH intersection lists
    std::vector<std::pair<NodeIndex,NodeIndex> > queue;
    queue.reserve(leaves.size());
    make_frontlist(std::make_pair(tree->getRootIndex(), false), queue);
    assert((int)queue.size() > 0);
    std::atomic<int> stop(0);

    int nq = queue.size();
    if (options.parallel)
    {
        tbb::parallel_for(tbb::blocked_range<int>(0,nq), [&](const tbb::blocked_range<int>& r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                collide(V0, V1, P, queue[i].first, queue[i].second, stop);
            }
        });
    }
    else
    {
        for (int i=0; i<nq; ++i)
        {
            collide(V0, V1, P, queue[i].first, queue[i].second, stop);
        }
    }
}

template <typename T, int DIM, int PDIM>
void BVHTree<T,DIM,PDIM>::traverse(BVHTraverse<T,DIM> *traverser) const
{
    const Eigen::KdBVH<T,DIM,LeafType> &tree_ref = *tree.get();
    Eigen::BVIntersect(tree_ref, *traverser);
}

template <typename T, int DIM, int PDIM>
void BVHTree<T,DIM,PDIM>::make_frontlist(const NodeIndex &idx,
    std::vector<std::pair<NodeIndex,NodeIndex> > &queue) const
{
    if (idx.second) // is leaf node
        return;

    NodeIndex l, r;
    get_children(idx, l, r);
    assert(l.first >= 0 && r.first >= 0);
    make_frontlist(l, queue);
    make_frontlist(r, queue);

    // if both leaf, should not have same idx
    assert((l.second==r.second) ? (l.first!=r.first) : true);
    
    queue.emplace_back(l, r);
}

template <typename T, int DIM, int PDIM>
void BVHTree<T,DIM,PDIM>::get_children(const NodeIndex &idx, NodeIndex &l, NodeIndex &r) const
{
    typedef typename Eigen::KdBVH<T, DIM, BVHLeaf<T,DIM> >::VolumeIterator VolumeIter;
    typedef typename Eigen::KdBVH<T, DIM, BVHLeaf<T,DIM> >::ObjectIterator ObjectIter;
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
    int num_children = 0;
    bool is_l = true;
    for (; vBegin != vEnd; ++vBegin)
    {
        if (is_l) { l = NodeIndex(*vBegin, false); is_l = false; }
        else { r = NodeIndex(*vBegin, false); }
        num_children++;
    }
    for (; oBegin != oEnd; ++oBegin)
    {
        if (is_l) { l = NodeIndex(oBegin->idx, true); is_l = false; }
        else { r = NodeIndex(oBegin->idx, true); }
        num_children++;
    }
    assert(num_children == 2); // assumes binary tree
}

template <typename T, int DIM, int PDIM>
const typename BVHLeaf<T,DIM>::BoxType& BVHTree<T,DIM,PDIM>::get_box(const NodeIndex &idx) const
{
    if (idx.second) { return leaves[idx.first].box; }
    return tree->getVolume(idx.first);
}

template <typename T, int DIM, int PDIM>
void BVHTree<T,DIM,PDIM>::collide(
    const T* V0, const T* V1, const int *P,
    const NodeIndex &left, const NodeIndex &right,
    std::atomic<int> &stop) const
{
    // Check if we should stop traversing, only possible if discrete test
    if (options.discrete && append_discrete != nullptr)
    {
        if (stop.load() > 0) { return; }
    }

    // if both leaf, not eq
    assert(left.first >= 0 && right.first >= 0);
    assert((left.second == right.second) ? (left.first != right.first) : true);

    // Check if nodes intersect or inactive
    if (!boxes_intersect(left,right))
        return;

    // Both are leaf nodes
    if (left.second && right.second)
    {
        if (options.continuous)
        {
            std::array<PairType, NumCandidates> pairs;
            get_candidates(left.first, right.first, P, pairs);
            for (size_t i=0; i<NumCandidates; ++i)
            {
                if (pairs[i].second < 0) {
                    continue;
                }

                T toi = -1;
                if (narrow_phase != nullptr) { toi = narrow_phase(pairs[i].first, pairs[i].second); }
                else { toi = default_narrow_phase(V0, V1, pairs[i].first, pairs[i].second); }
                if (toi >= 0 && append_pair != nullptr)
                {
                    append_pair(pairs[i].first, pairs[i].second, toi);
                }
            }
        }

        if (options.discrete)
        {
            int p0[DIM], p1[DIM];
            for (int i=0; i<DIM; ++i)
            {
                p0[i] = P[left.first*DIM+i];
                p1[i] = P[right.first*DIM+i];
            }
            bool d_hit = default_discrete_test(V1, p0, p1);
            if (d_hit && append_discrete != nullptr)
            {
                bool stop_traverse = append_discrete(left.first, right.first);
                if (stop_traverse) { stop++; }
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
    }
    else
    {
        NodeIndex l, r;
        get_children(left, l, r);
        collide(V0, V1, P, l, right, stop);
        collide(V0, V1, P, r, right, stop);
    }
}

template <typename T, int DIM, int PDIM>
void BVHTree<T,DIM,PDIM>::get_candidates(int p0, int p1, const int *P,
    std::array<PairType, NumCandidates> &pairs) const
{
    const LeafType& l0 = leaves[p0];
    const LeafType& l1 = leaves[p1];

    for (size_t i=0; i<NumCandidates; ++i)
    {
        pairs[i].first.setZero();
        pairs[i].second = -1;
    }

    int f0[PDIM], f1[PDIM];
    for (int i=0; i<PDIM; ++i)
    {
        f0[i] = P[p0*PDIM+i];
        f1[i] = P[p1*PDIM+i];
    }

    auto vf_shared_vertex = [](int vi, const int *fi)->bool
    {
        for (int i=0; i<PDIM; ++i)
            if (vi == fi[i])
                return true;

        return false;
    };

    int pair_index = 0;

    // VF, f0 -> f1
    for (int i=0; i<PDIM; ++i)
    {
        if (!l0.v[i]) // v not represented
        {
            pairs[pair_index++].second = -1;
            continue;
        }

        if (vf_shared_vertex(f0[i], f1))
        {
            pairs[pair_index++].second = -1;
            continue;
        }
    
        PairType &pair = pairs[pair_index++];
        pair.second = 1; // is_vf

        pair.first = Eigen::Vector4i(-1,-1,-1,-1);
        pair.first[0] = f0[i];
        for (int j=0; j<PDIM; ++j)
            pair.first[j+1] = f1[j];

        if (filter_pair != nullptr)
            if (filter_pair(pair.first, pair.second))
                pair.second = -1;
    }

    // VF, f1 -> f0
    for (int i=0; i<PDIM; ++i)
    {
        if (!l1.v[i]) // v not represented
        {
            pairs[pair_index++].second = -1;
            continue;
        }

        if (vf_shared_vertex(f1[i], f0))
        {
            pairs[pair_index++].second = -1;
            continue;
        }

        PairType &pair = pairs[pair_index++];
        pair.second = 1; // is_vf

        pair.first = Eigen::Vector4i(-1,-1,-1,-1);
        pair.first[0] = f1[i];
        for (int j=0; j<PDIM; ++j)
            pair.first[j+1] = f0[j];

        if (filter_pair != nullptr)
            if (filter_pair(pair.first, pair.second))
                pair.second = -1;
    }

    if (PDIM != 3)
        return;

    // EE
    for (int i=0; i<3; ++i)
    {
        if (!l0.e[i]) { continue; }
        for (int j=0; j<3; ++j)
        {
            if (!l1.e[j])
            {
                pairs[pair_index++].second = -1;
                continue;
            }

            Eigen::Vector4i sten(i, (i+1)%3, j, (j+1)%3);

            // shares vertex?
            if (f0[sten[0]]==f1[sten[2]] || f0[sten[0]]==f1[sten[3]] ||
                f0[sten[1]]==f1[sten[2]] || f0[sten[1]]==f1[sten[3]])
            {
                pairs[pair_index++].second = -1;
                continue;
            }

            PairType &pair = pairs[pair_index++];
            pair.second = 0; // is not vf
            pair.first = Eigen::Vector4i(f0[sten[0]], f0[sten[1]], f1[sten[2]], f1[sten[3]]);

            if (filter_pair != nullptr)
                if (filter_pair(pair.first, pair.second))
                    pair.second = -1;
        }
    }
}

template <typename T, int DIM, int PDIM>
bool BVHTree<T,DIM,PDIM>::boxes_intersect(const NodeIndex &left, const NodeIndex &right) const
{
    const typename BVHLeaf<T,DIM>::BoxType &lbox = get_box(left);
    const typename BVHLeaf<T,DIM>::BoxType &rbox = get_box(right);

    // Check if both branches inactive
    if (!lbox.active && !rbox.active)
        return false;

    // Check nodes intersect
    if (!lbox.intersects(rbox))
        return false;

    // Check time varying boxes if not empty
    if (lbox.t0.isEmpty() || rbox.t0.isEmpty()) { return true; }
    if (lbox.t1.isEmpty() || rbox.t1.isEmpty()) { return true; }

    for (int i=0; i<DIM; ++i)
    {
        if((lbox.t0.min()[i] > rbox.t0.max()[i] && lbox.t1.min()[i] > rbox.t1.max()[i]) ||
        (lbox.t0.max()[i] < rbox.t0.min()[i] && lbox.t1.max()[i] < rbox.t1.min()[i]))
            return false;
    }

    return true;
}

template <typename T, int DIM, int PDIM>
T BVHTree<T,DIM,PDIM>::default_narrow_phase(const T* V0, const T* V1, const Eigen::Vector4i &sten, bool is_vf) const
{
    VecType verts0[4], verts1[4];
    for (int i=0; i<4; ++i)
    {
        assert(sten[i] >= 0);
        for (int j=0; j<DIM; ++j)
        {
            verts0[i][j] = V0[sten[i]*DIM+j];
            verts1[i][j] = V1[sten[i]*DIM+j];
        }
    }

    T toi = -1;
    int hit = 0;
    if (is_vf)
    {
        hit = NarrowPhaseACCD<T,DIM>::query_ccd_vf(verts0, verts1, options.vf_ccd_eta, toi);
    }
    else
    {
        hit = NarrowPhaseACCD<T,DIM>::query_ccd_ee(verts0, verts1, options.ee_ccd_eta, toi);
    }

    if (hit == 1)
        return toi;

    return -1;
}

template <typename T, int DIM, int PDIM>
bool BVHTree<T,DIM,PDIM>::default_discrete_test(const T* V, const int *p0, const int *p1) const
{
    // Shares vertex?
    for (int i=0; i<DIM; ++i)
        for (int j=0; j<DIM; ++j)
            if (p0[i] == p1[j])
                return false;

    if (DIM == 3)
    {
        Eigen::Matrix<T,3,1> v0[DIM], v1[DIM];
        for (int i=0; i<3; ++i)
        {
            for (int j=0; j<3; ++j)
            {
                v0[i][j] = V[p0[i]*DIM+j];
                v1[i][j] = V[p1[i]*DIM+j];
            }
        }
        return NarrowPhase<T,3>::discrete_tri_tri(v0, v1);
    }
    else if (DIM == 2)
    {
        Eigen::Matrix<T,2,1> v0[DIM], v1[DIM];
        for (int i=0; i<2; ++i)
        {
            for (int j=0; j<2; ++j)
            {
                v0[i][j] = V[p0[i]*DIM+j];
                v1[i][j] = V[p1[i]*DIM+j];
            }
        }
        return NarrowPhase<T,2>::discrete_edge_edge(v0, v1); 
    }

    return false;
}

} // end namespace ccd
} // end namespace mcl

template class mcl::ccd::BVHTree<double,3,3>;
template class mcl::ccd::BVHTree<double,2,3>;
template class mcl::ccd::BVHTree<float,3,3>;
template class mcl::ccd::BVHTree<float,2,3>;
