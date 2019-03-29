#include "../include/bvh.h"

BVH::BVH(vector<Vector3i> &triangles, vector<Vector3f> &points)
{
    cout << "BVH" << endl;
    m_triangles = triangles;
    m_points = points;

    cout << "Creation of the BVH..." << endl;
    build_node(m_nodes, m_triangles, 0, m_triangles.size()-1);
    cout << "      BVH built" << endl;
}

Vector3f BVH::meanTrianglePoint(int triangle)
{
    Vector3f center = Vector3f(0);
    for (unsigned int j=0; j<3; j++)
        center+=m_points[m_triangles[triangle](j)];
    return center/3.0;
}

unsigned int BVH::build_node(vector<BVHnode> &nodes, vector<Vector3i>& triangles, const unsigned int begin, const unsigned int end)
{
    BVHnode node;

    // construire la boite englobante des centres des triangles d'indices [begin .. end[

    //Contruction of the bounding box of the triangles from indices begin to end
    Vector3f mini = meanTrianglePoint(begin);
    Vector3f maxi = meanTrianglePoint(begin);
    for (unsigned int i=begin+1; i<end; i++)
    {
        Vector3f currentTriangleCenter = meanTrianglePoint(i);
        for (unsigned int j=0; j<3; j++)
        {
            mini(j) = (mini(j)<currentTriangleCenter(j)) ? mini(j) : currentTriangleCenter(j);
            maxi(j) = (maxi(j)>currentTriangleCenter(j)) ? maxi(j) : currentTriangleCenter(j);
        }
    }
    node.corner1=mini;
    node.corner2=maxi;

    //cout << "mini : " << mini << endl << "maxi : " << maxi << endl;
    //cout << "end - begin : " << end-begin << endl;
    if (end - begin <= 1)//One triangle left : leaf
    {
        node.left = -1;
        node.right = -1;
        for (unsigned int i=begin; i<end; i++)
            node.t.push_back(i);
        nodes.push_back(node);
        return nodes.size()-1;
    }

    Vector3f gap = maxi-mini;
    unsigned int axe = gap(1)>gap(0) ? (gap(1)>gap(2) ? 1 : 2) : (gap(0)>gap(2) ? 0 : 2);
    //cout << "axe : " << axe << endl;

    float cut = (maxi(axe) + mini(axe))/2.0;

    Vector3i * middle = partition(triangles.data() + begin, triangles.data() + end, predicate(axe, cut, m_points));
    unsigned int mid = distance(triangles.data(), middle);

    //cout << cut << " " << end-begin << " " << end << " " << begin << " " << mid << endl;
    if (mid==begin || mid==end)//Leaf
    {
        node.left = -1;
        node.right = -1;
        for (unsigned int i=begin; i<end; i++)
            node.t.push_back(i);
        nodes.push_back(node);
        return nodes.size()-1;
    }
    node.left = build_node(nodes, triangles, begin, mid);
    node.right = build_node(nodes, triangles, mid, end);

    nodes.push_back(node);
    return nodes.size()-1;
}
