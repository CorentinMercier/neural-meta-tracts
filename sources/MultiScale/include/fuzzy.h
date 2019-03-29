#ifndef FUZZY_H
#define FUZZY_H

#include <iostream>
#include <fstream>
#include <limits>
#include "fiber.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/Dense>

#define WEIGHTEDCURRENTSSIMILARITY 0
#define QBMETRIC 1
#define WEIGHTEDCURRENTS 2
#define MINOFCLOSESTDISTANCES 3
#define VARIFOLD 4
#define WEIGHTEDCURRENTSFUZZY 5


using namespace std;
using namespace Eigen;

//static unsigned int metric;

class Fuzzy
{
public:
    Fuzzy(string data);
    ~Fuzzy();

    unsigned int getValueFromPoint(Vector3f p) const;
    float getValueFromSegment(Vector3f p1, Vector3f p2) const;
    float getValueFromVolume(Vector3f pointBefore, Transfo T1, Transfo T2, Vector3f pointAfter) const;
    float getDistanceToMask(Vector3f pt);
    float getDistanceToMaskFromSurface(Transfo T, Vector3f tangent);
    void computeDistanceField();

private:
    string m_filenameData;

    Matrix4f m_affineTransfo;
    Matrix4f m_affineTransfoInverse;
    unsigned int m_dimX;
    unsigned int m_dimY;
    unsigned int m_dimZ;
    Vector3f voxelDimensionsInCm;

    vector<int> m_data;
    vector<Vector3f> m_contourPts;
    vector<float> m_distance;
};

#endif // FUZZY_H
