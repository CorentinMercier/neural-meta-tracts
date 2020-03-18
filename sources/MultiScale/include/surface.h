#ifndef SURFACE_H
#define SURFACE_H

#include <iostream>
#include "../../opengl/openglutils.h"
#include "../include/bvh.h"
#include "libraries.h"

#ifdef USE_VTK
#include <vtkGenericDataObjectReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkMapper.h>
#endif

#include <glm/glm.hpp>

//#include <embree3/rtcore.h>
//#include <embree3/rtcore_ray.h>

#ifdef USE_VTK
using namespace vtk;
#endif
using namespace std;

class Fiber;

class Surface
{
public:
	/* Load surface from file surface_file */
	Surface(string surface_file);
	/* Destructor of Surface objects */
	~Surface();
	/* Return the vaoSurface pointer */
	GLuint getVaoSurface() {return m_vaoSurface;}
	/* Return the eboSurface pointer */
	GLuint getEboSurface() {return m_eboSurface;}
	/* Return the vboSurface pointer */
	GLuint getVboSurface() {return m_vboSurface;}
	/* Return the number of points of the surface*/
	unsigned int getNbPoints() {return static_cast<unsigned int>(m_surfacePoints.size());}
	/* Return the number of triangles of the surface*/
	unsigned int getNbTriangles() {return static_cast<unsigned int>(m_surfaceTriangles.size());}
	/* Return true if the point is inside the surface */
	bool isInside(Vector3f point);
	/* Return true if there is an intersection with the surface, returning the intersection point */
	bool intersect(Vector3f origin, Vector3f direction, Vector3f &intersection);

private:
	/* Compute a BVH with the triangles to make intersection computation faster */
	void computeBVH();

	vector<Vector3f> m_surfacePoints;
	vector<Vector3i> m_surfaceTriangles;
	vector<Vector3f> m_surfaceNormals;

	GLuint m_vaoSurface;
	GLuint m_eboSurface;
	GLuint m_vboSurface;

	//For raytracing
//    bool launchRayInside(float dirX, float dirY, float dirZ, float tfar, Vector3f point);
//    bool launchRay(Vector3f direction, float tfar, Vector3f origin, Vector3f &intersection);
//    void rtcInitialize();
//    RTCDevice m_rtcdevice;
//    RTCScene m_rtcscene;
//    RTCGeometry m_rtcgeometry;
//    RTCIntersectContext m_rtccontext;// = RTC_INTERSECT_CONTEXT_FLAG_COHERENT;
//    BVH *bvh;
};

#endif // SURFACE_H
