// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "MCL/BVHTree.hpp"
#include "MCL/NarrowPhase.hpp"
#include "MCL/ccd_internal/Timer.hpp"

#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_vector.h>
#include <atomic>
#include <iostream>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h>
#include <igl/readOBJ.h>

// Collects AABB nodes for rendering
class NodeCollector : public mcl::BVHTraverse<double,3>
{
public:
    using mcl::BVHTraverse<double,3>::VolumeType;
    using mcl::BVHTraverse<double,3>::ObjectType;
    std::vector<Eigen::Vector3d> edges0, edges1;
    bool intersectVolume(const VolumeType &volume);
    bool intersectObject(const ObjectType&) { return false; }
};

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> RowMatrixXd;
typedef Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> RowMatrixXi;

// Load and render a dillo
int main(int, char**)
{
    mcl::Timer timer;

    using namespace Eigen;
    std::cout << "Loading the mesh: " << std::flush;
    MatrixXd V_;
    MatrixXi F_;
    if (!igl::readPLY(MCLCCD_ROOT_DIR "/examples/data/armadillo.ply", V_, F_)) // no intersections but small elements
    //if (!igl::readOBJ(MCLCCD_ROOT_DIR "/examples/data/hand.obj", V_, F_)) // has discrete self-intersections
        return EXIT_FAILURE;

    std::cout << "Num faces: " << F_.rows() << ", num verts: " << V_.rows() << std::endl;
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_, F_);
    viewer.data().set_face_based(true);
    viewer.core().is_animating = true;

    RowMatrixXd V = V_;
    RowMatrixXi F = F_;

    // Initialize the tree
    std::cout << "Building the tree: " << std::flush; 
    timer.start();
    mcl::BVHTree<double,3> tree;
    tree.options.box_eta = std::numeric_limits<float>::epsilon();
    tree.update(V, V, F);
    std::cout << timer.get_ms() << " ms" << std::endl;

    // Traverse BVH for intersections
    std::cout << "Performing CCD: " << std::flush; 
    timer.start();
    tbb::concurrent_vector<std::pair<Eigen::Vector4i,int>> pairs;
    tbb::concurrent_unordered_set<int> discrete;
    tree.append_pair = [&](const Eigen::Vector4i &sten, int type, const double &toi)->void
    {
        (void)(toi);
        pairs.emplace_back(sten, type);
    };
    tree.append_discrete = [&](int p0, int p1)->bool
    {
        discrete.emplace(p0);
        discrete.emplace(p1);
        return false; // don't stop traversing
    };
    tree.traverse(V, V, F);
    std::cout << timer.get_ms() << " ms" << std::endl;
    std::cout << "Num discrete isect: " << discrete.size() << std::endl;
    std::cout << "Num ccd isect: " << pairs.size() << std::endl;

    // Render the tree
    std::cout << "Gathering boxes to render: " << std::flush;
    timer.start();
    NodeCollector collector;
    tree.traverse(&collector);
    std::cout << timer.get_ms() << " ms" << std::endl;
    MatrixXd e0(collector.edges0.size(), 3);
    MatrixXd e1(collector.edges1.size(), 3);
    int ne = e0.rows();
    for (int i=0; i<ne; ++i)
    {
        e0.row(i) = collector.edges0[i];
        e1.row(i) = collector.edges1[i];
    }

    // Per-face colors
    MatrixXd C = MatrixXd::Ones(F.rows(), 3) * 0.7;
    for (tbb::concurrent_unordered_set<int>::const_iterator it = discrete.begin();
        it != discrete.end(); ++it) { C.row(*it) = RowVector3d(1,0,0); }

    viewer.data().add_edges(e0, e1, RowVector3d(1,0,0));
    viewer.data().set_colors(C);

    // Launch viewer
    std::cout << "Launching Viewer" << std::endl;
    viewer.launch();
    return EXIT_SUCCESS;
}


bool NodeCollector::intersectVolume(const VolumeType &volume)
{
    using namespace Eigen;
    Vector3d min = volume.min();
    Vector3d max = volume.max();

    // Bottom quad
    Vector3d a = min;
    Vector3d b(max[0], min[1], min[2]);
    Vector3d c(max[0], min[1], max[2]);
    Vector3d d(min[0], min[1], max[2]);
    // Top quad
    Vector3d e(min[0], max[1], min[2]);
    Vector3d f(max[0], max[1], min[2]);
    Vector3d g = max;
    Vector3d h(min[0], max[1], max[2]);

    // make edges
    // bottom
    edges0.emplace_back(a); edges1.emplace_back(b);
    edges0.emplace_back(a); edges1.emplace_back(d);
    edges0.emplace_back(c); edges1.emplace_back(b);
    edges0.emplace_back(c); edges1.emplace_back(d);
    // top
    edges0.emplace_back(e); edges1.emplace_back(f);
    edges0.emplace_back(e); edges1.emplace_back(h);
    edges0.emplace_back(g); edges1.emplace_back(f);
    edges0.emplace_back(g); edges1.emplace_back(h);
    // columns
    edges0.emplace_back(d); edges1.emplace_back(h);
    edges0.emplace_back(min); edges1.emplace_back(e);
    edges0.emplace_back(b); edges1.emplace_back(f);
    edges0.emplace_back(c); edges1.emplace_back(max);
    return true;
}

