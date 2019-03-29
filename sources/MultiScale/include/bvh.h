#ifndef BVH_H
#define BVH_H

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct BVHnode
{
    Vector3f corner1;//Min of axes x, y and z
    Vector3f corner2;//Max of axes x, y and z
    int left;//if -1, no left son
    int right;//if -1, no right son
    vector<unsigned int> t;//triangles indices if it is a leaf
};

class BVH
{
public:
    BVH(vector<Vector3i>& triangles, vector<Vector3f>& points);

    unsigned int build_node(vector<BVHnode>& nodes, vector<Vector3i> &triangles, const unsigned int begin, const unsigned int end);

private:
    Vector3f meanTrianglePoint(int triangle);
    vector<Vector3i> m_triangles;
    vector<Vector3f> m_points;
    vector<BVHnode> m_nodes;
};

struct predicate
{
    unsigned int axe;
    float cut;
    vector<Vector3f> points;
    predicate(const unsigned int _axe, const float _cut, vector<Vector3f>& _points) : axe(_axe), cut(_cut), points(_points) {}
    bool operator() (const Vector3i t) const
    {
        int center = 0;
        for (unsigned int j=0; j<3; j++)
            center+=points[t(j)](axe);
        return (center/3>cut);
    }
};

#endif // BVH_H
