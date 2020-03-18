#ifndef FIBER_H
#define FIBER_H

#include "iostream"
#include "priorityqueue.h"
#include <set>
#include "surface.h"

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>

#define PROFILE_PER_POINT true
#define MIN_RADIUS 0.05

using namespace std;

struct Transfo
{
	glm::vec3 center;
	glm::vec3 rotation;
	glm::vec3 ellipseOrientation;
	float a;
	float b;
};

struct Profile
{
	Transfo ellipse;
};

//class Surface;

class Fiber
{
public:
	Fiber(std::vector<glm::vec3> points, unsigned int size);
	Fiber();
	//Fiber& operator=(Fiber f);
	~Fiber();

	unsigned int size() const {return m_data.size();}

	glm::vec3 operator[](unsigned int element) const;

	void setTangents(std::vector<glm::vec3> t);
	float length();
	float cylinderLength();
	void PutInSameDirection(const Fiber &fiber);

	float getWidth(unsigned int i);
	void setWidth(vector<float> & width);

	Transfo getProfileTransform(unsigned int i);
	void setProfileTransform(vector<Transfo> & profileTransformation);

	void setProfile(Profile profile){m_profile=profile;}
	Profile getProfile(){return m_profile;}
	void setTransform(int i, glm::mat3 matrice){m_transfoMatrix[i]=matrice;}
	glm::mat3 getTransform(int i){return m_transfoMatrix[i];}

	void setValue(float value){m_value=value;}
	void setCylinderValue(float value){m_Cvalue=value;}
	float getPoint(int i, int j);
	glm::vec3 getCenter(int i);
	glm::vec3 getTangent(int i);
	float getValue(){return m_value;}
	float getValueOfCylinder(){return m_Cvalue;}
	float getMaxDistEndPoints(Fiber &fiber);

	void resample(unsigned int nb_points);
	void resampleV2(unsigned int nb_points);
	bool is_flipped(){return m_flipped;}
	void flip(){m_flipped=!m_flipped;}

	void computeCentersAndTangents();

//    void addPaire(paire p) {m_paires.push(p);}
//    paire getPaire(){return m_paires.top();}
//    void removePaire(){m_paires.pop();}
//    bool thereIsAPaire(){return !(m_paires.empty());}

	bool is_outlier(){return m_isOutlier;}
	void setOutlier(bool outlier){m_isOutlier=outlier;}
	bool is_valid(){return (m_valid && ((m_fiberMode && m_value>=m_fuzzy_limit) || (!m_fiberMode && m_Cvalue>=m_fuzzy_limit)));}
	bool is_validCylinder(){return (m_valid && m_Cvalue>=m_fuzzy_limit);}
	void setValid(bool v){m_valid=v;}// if (!m_valid) clearQueue(m_paires);}
	void setToFiberMode(bool mode){m_fiberMode = mode;}

	void addNeighbours(vector<unsigned int> new_neighbours);
	set<unsigned int> m_neighbours;

	void addToHierarchy(unsigned int new_element) {m_hierarchy.insert(new_element);}
	void addToHierarchy(vector<unsigned int> new_elements);
	vector<unsigned int> getHierarchy();

	unsigned int getNbFibHierarchy() {return m_NbFibHierarchy;}
	void addFibHierarchy(unsigned int nbFibHierarchy) {m_NbFibHierarchy+=nbFibHierarchy;}

	float computeVolume();
	float getSelfScalarProduct(float lambdag);
	float getSelfVarifoldScalarProduct(float sigma);

	void setColor(unsigned char r, unsigned char g, unsigned char b){m_red=r; m_green=g; m_blue=b;}
	glm::vec3 getColor(){return glm::vec3((float)m_red/255, (float)m_green/255, (float)m_blue/255);}
	void setCylinderColor(unsigned char r, unsigned char g, unsigned char b){m_Cred=r; m_Cgreen=g; m_Cblue=b;}
	glm::vec3 getCylinderColor(){return glm::vec3((float)m_Cred/255, (float)m_Cgreen/255, (float)m_Cblue/255);}

	void geometricSmoothing();
	void pointsSmoothing();

	void moveExtremities(glm::vec3 deviation);

	/* cut the fiber according to the surface given as argument */
	Fiber cut(Surface *surface);

	void setFuzzyLimit(float limit){m_fuzzy_limit = limit;}

	void clearCentersAndTangents(){m_centers.clear(); m_tangents.clear();}

private:
	//Data
	float m_length = 0;
	float m_Clength = 0;
	float m_value = 0.0f;
	float m_Cvalue = 0.0f;
	bool m_fiberMode = true;
	/*Fiber color*/
	unsigned char m_red;
	unsigned char m_green;
	unsigned char m_blue;
	/*Cylinder color*/
	unsigned char m_Cred;
	unsigned char m_Cgreen;
	unsigned char m_Cblue;
	//unsigned int m_nb_points;
	std::vector<glm::vec3> m_data;
	std::vector<glm::vec3> m_centers;
	std::vector<glm::vec3> m_tangents;
	bool m_computedCenterAndTangents=false;
	vector<float> m_width;//Generalized cylinders with single radius (circles)
	//m_profileTransformation starts with sx, sy, sz (scaling parameter), then alpha, beta, gamma (angle of rotation on the three axis)
	vector<Transfo> m_profileTransformation;//Generalized cylinders built with a transformation matrix to change the circle profile

	//Profile along the fiber, defined by a profile and a vector of transformation matrices
	Profile m_profile;
	vector<glm::mat3> m_transfoMatrix;

	bool m_flipped=false;
	bool m_valid=true;
	bool m_isOutlier=false;
	//priority_queue<paire, vector<paire>, comparisonPaires> m_paires;
	set<unsigned int> m_hierarchy;
	unsigned int m_NbFibHierarchy=1;
	float m_volume=-1;
	float m_selfScalarProduct;
	bool m_selfScalarProductComputed=false;

	//Functions
	vector<float> cumulateLength();

	//Fuzzy limit
	float m_fuzzy_limit = 0.0f;
};

#endif // FIBER_H
