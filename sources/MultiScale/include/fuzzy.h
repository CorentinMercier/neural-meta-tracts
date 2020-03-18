#ifndef FUZZY_H
#define FUZZY_H

#include <iostream>
#include <fstream>
#include <limits>
#include "fiber.h"

#define WEIGHTEDCURRENTSSIMILARITY 0
#define QBMETRIC 1
#define WEIGHTEDCURRENTS 2
#define MINOFCLOSESTDISTANCES 3
#define VARIFOLD 4
#define WEIGHTEDCURRENTSFUZZY 5


using namespace std;

class Fuzzy
{
public:
	Fuzzy(string data);
	~Fuzzy();

	unsigned int getValueFromPoint(glm::vec3 p) const;
	float getValueFromSegment(glm::vec3 p1, glm::vec3 p2) const;
	float getValueFromVolume(glm::vec3 pointBefore, Transfo T1, Transfo T2, glm::vec3 pointAfter) const;
	float getDistanceToMask(glm::vec3 pt);
	float getDistanceToMaskFromSurface(Transfo T, glm::vec3 tangent);
	void computeDistanceField();

private:
	string m_filenameData;

	glm::mat4 m_affineTransfo;
	glm::mat4 m_affineTransfoInverse;
	unsigned int m_dimX;
	unsigned int m_dimY;
	unsigned int m_dimZ;
	glm::vec3 voxelDimensionsInCm;

	vector<int> m_data;
	vector<glm::vec3> m_contourPts;
	vector<float> m_distance;
};

#endif // FUZZY_H
