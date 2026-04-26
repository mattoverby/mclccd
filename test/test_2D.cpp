// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include "MCL/BVHTree.hpp"
#include <stdexcept>

void
assert_true(bool condition, const std::string& msg)
{
    if (!condition) {
        throw std::runtime_error(msg);
    }
}

int
main(int, char**)
{
    // Three triangles, discrete
    {
        mcl::ccd::BVHTree<double, 2> tree;
        tree.options.parallel = false;
        Eigen::MatrixXd V(9, 2);
        V << 0, 0, 0, 1, 1, 0, 0.5, 0.5, 0.5, 1.5, 1.5, 0.5, 10, 10, 10, 11, 11, 10;
        Eigen::MatrixXi F(3, 3);
        F << 0, 1, 2, 3, 4, 5, 6, 7, 8;

        std::vector<Eigen::Vector2i> pairs;
        tree.update(V, V, F);
        tree.append_discrete = [&](int p0, int p1) -> bool {
            pairs.emplace_back((p0 > p1) ? Eigen::Vector2i(p0, p1) : Eigen::Vector2i(p1, p0));
            return false;
        };
        tree.traverse(V, V, F);

        assert_true(pairs.size() == 1, "didn't detect single intersection: " + std::to_string(pairs.size()));
        assert_true(pairs[0][0] == 1 && pairs[0][1] == 0, "bad intersection prim idx");
    }

    return EXIT_SUCCESS;
}