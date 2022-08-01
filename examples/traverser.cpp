// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "MCL/BVHTree.hpp"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>

Eigen::Vector3d point_triangle_barys(
    const Eigen::Vector2d &p,
    const Eigen::Vector2d &p0,
    const Eigen::Vector2d &p1,
    const Eigen::Vector2d &p2)
{
	Eigen::Vector2d v0 = p1 - p0, v1 = p2 - p0, v2 = p - p0;
	float d00 = v0.dot(v0);
	float d01 = v0.dot(v1);
	float d11 = v1.dot(v1);
	float d20 = v2.dot(v0);
	float d21 = v2.dot(v1);
	float invDenom = 1.0 / (d00 * d11 - d01 * d01);
	Eigen::Vector3d r;
	r[1] = (d11 * d20 - d01 * d21) * invDenom;
	r[2] = (d00 * d21 - d01 * d20) * invDenom;
	r[0] = 1.0 - r[1] - r[2];
	return r;
}

class PointInTriangle : public mcl::BVHTraverse<double,2>
{
public:
    using mcl::BVHTraverse<double,2>::VolumeType;
    using mcl::BVHTraverse<double,2>::ObjectType;
    Eigen::Vector2d pt;
    const Eigen::MatrixXd &V;
    const Eigen::MatrixXi &F;
    std::vector<int> in_tri;

    PointInTriangle(
        const Eigen::Vector2d &pt_,
        const Eigen::MatrixXd &V_,
        const Eigen::MatrixXi &F_) :
        pt(pt_), V(V_), F(F_) {}

    bool intersectVolume(const VolumeType &v)
    {
        return v.contains(pt);
    }

    bool intersectObject(const ObjectType &o)
    {
        Eigen::Vector3i f = F.row(o.idx);
        Eigen::Vector3d b = point_triangle_barys(
            pt, V.row(f[0]), V.row(f[1]), V.row(f[2]));
        if (b.minCoeff() >= 0) { in_tri.emplace_back(o.idx); }
        return false; // return true to stop traversing
    }
};

// Load and render a dillo
int main(int, char**)
{
    using namespace Eigen;
    MatrixXd V, TC, CN, FN;
    MatrixXi F, FTC;
    std::string obj = MCLCCD_ROOT_DIR "/examples/data/cathead.obj";
    igl::readOBJ(obj, V, TC, CN, F, FTC, FN);
    if (!TC.rows() || !FTC.rows()) {
        return EXIT_FAILURE;
    }

    AlignedBox<double,2> box;
    for (int i=0; i<(int)V.rows(); ++i) {;
        box.extend(Vector2d(TC(i,0), TC(i,1)));
    }

    // Random vertices inside the bounds
    // of the input mesh.
    MatrixXd pts = MatrixXd::Random(10,3);
    MatrixXd ptsC = MatrixXd::Zero(10,3);
    pts.array() += 1.0;
    pts.array() /= 2.0;
    pts.col(2).array() = 0;
    pts.col(0).array() *= box.sizes()[0];
    pts.col(1).array() *= box.sizes()[1];
    pts.col(0).array() += box.min()[0];
    pts.col(1).array() += box.min()[1];

    // Build the tree
    mcl::BVHTree<double,2> tree;
    tree.update(TC, TC, FTC);

    // Traverse and color intersected faces
    MatrixXd FC = MatrixXd::Zero(FTC.rows(), 3);
    FC.col(2).array() = 1;
    for (int i=0; i<(int)pts.rows(); ++i)
    {
        Vector2d pt = Vector2d(pts(i,0), pts(i,1));
        PointInTriangle t(pt, TC, FTC);
        tree.traverse(&t);
        for (int j=0; j<(int)t.in_tri.size(); ++j)
        {
            FC.row(t.in_tri[j]) = RowVector3d(1,0,0);
        }
    }

    // Launch viewer. Render faces that
    // contain a point in red.
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(TC, FTC);
    viewer.data().set_colors(FC);
    viewer.data().add_points(pts, ptsC);
    viewer.launch();
    return EXIT_SUCCESS;
}