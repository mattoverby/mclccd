# mclccd

Continuous collision detection for 2D and 3D triangular meshes.

## Build

```sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

## Example Collision Detection

Collisions are gathered through callbacks (append_pair and append_discrete), the former being CCD and latter being triangle-triangle (3D) or edge-edge (2D). These callbacks are called in parallel during the BVH traversal.


```cpp
Eigen::MatrixXd V0 = ... // nv x 2 or 3 vertices at t=0
Eigen::MatrixXd V1 = ... // nv x 2 or 3 vertices at t=1
Eigen::MatrixXi F = ... // nf x 3 triangles
mcl::BVHTree<double,3> tree;
tree.update(V0, V1, F); // creates or updates BVH
tree.append_pair = [&](const Eigen::Vector4i &sten, int type, const double &toi)->void
{
    contacts.emplace_back(sten, type, toi);
};
tree.append_discrete = [&](int p0, int p1)->bool
{
    discrete.emplace(p0, p1); // tri-tri or edge-edge
    return false; // return true to stop traversing
};
tree.traverse(V0, V1, F); // perform collision detection
```

By default, CCD uses [ACCD](https://doi.org/10.1145/3450626.3459767) narrow phase kernels, a form of conservative advancement that robustly supports large thickness/gaps. The CCD test can be pre-filtered and replaced using (optional) custom functions via callbacks:

```cpp
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

### General BVH Traversal

Apart from continuous and discrete collision detections, traversers can be used to walk down the BVH. Two functions need to be implemented in the derived class to determine if the branch should be traversed (intersectVolume) and to process the leaf (intersectObject). See [examples/traverser.cpp](examples/traverser.cpp) for an example.