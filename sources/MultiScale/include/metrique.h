#ifndef METRIQUE_H
#define METRIQUE_H

#include "iostream"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/Dense>

#include "fiber.h"
#include "fuzzy.h"

//static unsigned int metric;

using namespace Eigen;
using namespace std;

class Metrique
{
public:
    Metrique(unsigned int metric);
    ~Metrique();

    float MC(Fiber &fiber1, Fiber &fiber2);
    float SC(Fiber &fiber1, Fiber &fiber2);
    float LC(Fiber &fiber1, Fiber &fiber2);
    float MDF(Fiber &fiber1, Fiber &fiber2);
    float PDM(Fiber &fiber1, Fiber &fiber2, float sigma);
    float VarifoldDistance(Fiber &fiber1, Fiber &fiber2, float sigma);
    float CurrentDistance(Fiber &fiber1, Fiber &fiber2, float sigma);
    float WeightedCurrentsDistance(Fiber &fiber1, Fiber &fiber2, float lambdac, float lambdab, float lambdag);
    float scalarProductWeightedCurrentsParallel(Fiber &fiber1, Fiber &fiber2, float lambdac, float lambdab, float lambdag);
    float averageOfPointwiseEuclideanMetric(Fiber &fiber1, Fiber &fiber2); //Metric used by Quickbundles
    float metriqueTest(Fiber fiber1, Fiber fiber2);

    float endpointsDistance(Fiber &fiber1, Fiber &fiber2);
    Fiber mergeOfFibers(Fiber& fiber1, Fiber& fiber2);
    Fiber mergeOfFibersV2(Fiber& fiber1, Fiber& fiber2);

    void setMetric(unsigned int value){m_metric = value;}

private:
    unsigned int m_metric;
};




#endif // METRIQUE_H
