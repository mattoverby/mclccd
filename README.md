# mclccd

Continuous collision detection for 2D and 3D triangular meshes.

## Build

```sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

## Example

mcl::BVHTree is constructed given a list of vertices and triangles.  Collisions are gathered through
callbacks (append_pair and append_discrete), the former being CCD and latter being triangle-triangle (3D) or edge-edge (2D). These callbacks are called in parallel during the BVH traversal, and can be disabled in mcl::BVHTree::options.


```
Eigen::MatrixXd V0 = ... // vertices at t=0
Eigen::MatrixXd V1 = ... // vertices at t=1
Eigen::MatrixXi F = ... // nf x 3 triangles
mcl::BVHTree<double,3> tree;
tree.update(V0, V1, F); // updates BVH
tbb::concurrent_vector<MyContactType> contacts;
tbb::concurrent_unordered_set<std::pair<int,int>> discrete;
tree.append_pair = [&](const Eigen::Vector4i &sten, int type, const double &toi)->void
{
    // Gather CCD contacts through callback:
    contacts.emplace_back(sten, type, toi);
};
tree.append_discrete = [&](int p0, int p1)->bool
{
    // Gather intersected triangle-triangle at t=1
    // Return true to stop traversing BVH
    discrete.emplace(p0, p1);
    return false;
};
tree.traverse(V0, V1, F); // perform collision detection
```

Narrow phase kernels use [ACCD](https://doi.org/10.1145/3450626.3459767), a form of conservative advancement, by default. The CCD test can be pre-filtered and replaced using optional callbacks.

```
tree.filter_pair = [&](const Eigen::Vector4i &sten, bool is_vf)->bool
{
    // Return true to skip CCD test
    return my_skip_narrow_phase_test(...);
};
tree.narrow_phase = [&](const Eigen::Vector4i &sten, bool is_vf)->double
{
    // Use custom narrow phase kernels to compute time of impact (toi)
    // and return toi if intersected, otherwise return < 0.
    return my_narrow_phase_test(...);
};

```

Traversers can be used to walk down the BVH.

```
class PointInTriangle : public mcl::BVHTraverse<double,3>
{
public:
    using mcl::BVHTraverse<double,2>::VolumeType; // mcl::BoundingBox
    using mcl::BVHTraverse<double,2>::ObjectType; // mcl::BVHLeaf
    PointInTriangle(const Eigen::Vector2d &pt) : point(pt) {}
    bool intersectVolume(const VolumeType &v)
    {
        return v.intersects(point);
    }
    bool intersectObject(const ObjectType &o)
    {
        int tri_idx = o.idx;
        if (my_point_in_triangle_test(point, tri_idx))
        {
            isects.emplace_back(tri_idx);
        }
        return false; // Return true to stop traversing
    }
    Eigen::Vector2d point; // input point
    std::vector<int> isects; // output intersections
};
```

Which are invoked using the traverse function

```
Eigen::MatrixXd V = ... // nv x 2 vertices
Eigen::MatrixXi F = ... // nf x 3 triangles
mcl::BVHTree<double,2> tree;
tree.update(V, V, F);
PointInTriangle pt_in_tri(Eigen::Vector3d(0,0,0));
tree.traverse(&pt_in_tri);
for (size_t i=0; i<pt_in_tri.isects.size(); ++i)
{
    // process intersections...
}

```
