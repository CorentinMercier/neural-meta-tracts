#ifndef KDTREE_H
#define KDTREE_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

//#include <Eigen/Core>
//#include <Eigen/Sparse>
//#include <Eigen/Geometry>
//#include <Eigen/Dense>

//using namespace Eigen;
using namespace std;

//Structure pour les tris
struct compareElements
{
    compareElements(int dimension)
    {
        this->m_dimension=dimension;
    }
    bool operator() (vector<float> i, vector<float> j)
    {
        return i[m_dimension]<j[m_dimension];
    }
    int m_dimension;

};

class kdtree
{
public:
    kdtree(vector<vector<float> > elements, unsigned int depth=0);
    void addElement(vector<float> &element);
    void getNeighbours(vector<float> &objectif, vector<unsigned int> &listNeighbours, unsigned int numberOfNeighbours = 5);
    void getNeighboursByDistance(vector<float> &objectif, vector<unsigned int> &listNeighbours, float distanceCible);
    unsigned int count();
    void getElementAtDepth(vector<int> & element, unsigned int depth);

private:
    void getNeighboursRecursive(vector<float> &objectif, vector<vector<float>> & neighbours);
    void getNeighboursByDistanceRecursive(vector<float> &objectif, vector<vector<float> > &neighbours, float distanceCible);
    float distance(vector<float> & objectif);
    float maxDistance3DPoints(vector<float> &objectif);
    void addElementRecursive(vector<float> &element, unsigned int depth);
    unsigned int countRec(unsigned int actualCount);
    void getElementAtDepthRecursive(vector<int> & element, unsigned int depth, unsigned int currentDepth);
    float getMinAxis() {return m_axisMin;}
    float getMaxAxis() {return m_axisMax;}

    vector<float> m_location;
    kdtree* m_leftChild;
    kdtree* m_rightChild;
    unsigned int m_depth;
    float m_axisMin;
    float m_axisMax;
};

#endif // KDTREE_H
