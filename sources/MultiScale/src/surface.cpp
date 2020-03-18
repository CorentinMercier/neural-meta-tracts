#include "../include/surface.h"
#include "../include/fiber.h"

Surface::Surface(string surface_file)
{
#ifdef USE_VTK
	vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName(surface_file.c_str());
	reader->Update();
	vtkSmartPointer<vtkPolyData> mesh  = vtkSmartPointer<vtkPolyData>::New();
	mesh=reader->GetPolyDataOutput();
	unsigned int nbVertices = static_cast<unsigned int>(mesh->GetNumberOfPoints());
	unsigned int nbTriangles = static_cast<unsigned int>(mesh->GetNumberOfPolys());
	if (nbTriangles==0)
	{
		vtkSmartPointer<vtkTriangleFilter> filter  = vtkSmartPointer<vtkTriangleFilter>::New();
		filter->SetInputConnection(reader->GetOutputPort());
		filter->Update();
		mesh = filter->GetOutput();
	}
	nbTriangles = static_cast<unsigned int>(mesh->GetNumberOfPolys());

	cout << "This surface is composed of " << nbTriangles << " triangles and " << nbVertices << " vertices" << endl;

	//Acquiring the points of the surface
	m_surfacePoints.resize(nbVertices);
	m_surfaceNormals.resize(nbVertices);
	vector<int> nbNormals(nbVertices);
	for (unsigned int i=0; i<nbVertices; i++)
	{
		double p[3];
		mesh->GetPoint(i,p);
		m_surfacePoints[i] = Vector3f(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
		m_surfaceNormals[i] = Vector3f(0);
		nbNormals[i]=0;
	}
	vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	triangles = mesh->GetPolys();
	m_surfaceTriangles.resize(nbTriangles);
	for (unsigned int i=0; i<nbTriangles; i++)
	{
		triangles->GetNextCell(idList);
		m_surfaceTriangles[i] = Vector3i(static_cast<int>(idList->GetId(0)), static_cast<int>(idList->GetId(1)), static_cast<int>(idList->GetId(2)));
	}
	//Creation of the normal per vertex
	for (unsigned int i=0; i<nbTriangles; i++)
	{
		Vector3f currentNormal = (m_surfacePoints[static_cast<unsigned int>(m_surfaceTriangles[i](1))]-m_surfacePoints[static_cast<unsigned int>(m_surfaceTriangles[i](0))])
								 .cross(m_surfacePoints[static_cast<unsigned int>(m_surfaceTriangles[i](2))]-m_surfacePoints[static_cast<unsigned int>(m_surfaceTriangles[i](1))]);
		for (unsigned int j=0; j<3; j++)
		{
			nbNormals[static_cast<unsigned int>(m_surfaceTriangles[i](j))]++;
			m_surfaceNormals[static_cast<unsigned int>(m_surfaceTriangles[i](j))]+=currentNormal;
		}
	}
	for (unsigned int i=0; i<nbVertices; i++)
	{
		m_surfaceNormals[i]/=nbNormals[i];
		m_surfaceNormals[i].normalized();
	}
	glGenVertexArrays(1, &m_vaoSurface);
	glCreateBuffers(1, &m_vboSurface);
	glCreateBuffers(1, &m_eboSurface);
	glBindVertexArray(m_vaoSurface);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboSurface); PRINT_OPENGL_ERROR();
	glBufferData(GL_ARRAY_BUFFER, m_surfacePoints.size()*3*sizeof(float)*2, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, 0, m_surfacePoints.size()*3*sizeof(float), m_surfacePoints.data()); PRINT_OPENGL_ERROR();
	glBufferSubData(GL_ARRAY_BUFFER, m_surfacePoints.size()*3*sizeof(float), m_surfacePoints.size()*3*sizeof(float), m_surfaceNormals.data()); PRINT_OPENGL_ERROR();
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_eboSurface); PRINT_OPENGL_ERROR();
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_surfaceTriangles.size()*3*sizeof(int), m_surfaceTriangles.data(), GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

	cout << "Surface " << surface_file << " loaded" << endl;
	//bvh = new BVH(m_surfaceTriangles, m_surfacePoints);
	//rtcInitialize();
	cout << "Initialization complete" << endl;
#endif
}

Surface::~Surface()
{
//    rtcReleaseScene(m_rtcscene);
//    rtcReleaseDevice(m_rtcdevice);
}

void Surface::computeBVH()
{

}

//void Surface::rtcInitialize()
//{
//    m_rtcdevice = rtcNewDevice(nullptr);
//    m_rtcscene = rtcNewScene(m_rtcdevice);
//    m_rtcgeometry = rtcNewGeometry(m_rtcdevice, RTC_GEOMETRY_TYPE_TRIANGLE);
//    rtcSetSharedGeometryBuffer(m_rtcgeometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, m_surfacePoints.data(), 0, sizeof(Vector3f), m_surfacePoints.size());
//    rtcSetSharedGeometryBuffer(m_rtcgeometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, m_surfaceTriangles.data(),0, sizeof(Vector3i), m_surfaceTriangles.size());
//    rtcCommitGeometry(m_rtcgeometry);
//    rtcAttachGeometry(m_rtcscene, m_rtcgeometry);
//    //rtcReleaseGeometry(m_rtcgeometry);
//    rtcCommitScene(m_rtcscene);

//    rtcInitIntersectContext(&m_rtccontext);
//    //Vector3f* vertices = (Vector3f*) rtcSetNewGeometryBuffer(rtcgeometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(Vector3f), m_surfacePoints.size());
//    //Vector3i* triangles = (Vector3i*) rtcSetNewGeometryBuffer(rtcgeometry,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT3,sizeof(Vector3i),m_surfaceTriangles.size());

//}

//bool Surface::launchRayInside(float dirX, float dirY, float dirZ, float tfar, Vector3f point)
//{
//    RTCRayHit ray;
//    ray.ray.org_x = point(0);
//    ray.ray.org_y = point(1);
//    ray.ray.org_z = point(2);
//    ray.ray.dir_x = dirX;
//    ray.ray.dir_y = dirY;
//    ray.ray.dir_z = dirZ;
//    ray.ray.tnear = 0.0f;
//    ray.ray.tfar = tfar;
//    ray.ray.flags = 0;
//    ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
//    ray.hit.geomID = RTC_INVALID_GEOMETRY_ID;
//    rtcInitIntersectContext(&m_rtccontext);
//    rtcIntersect1(m_rtcscene, &m_rtccontext, &ray);
//    bool intersection = (ray.ray.tfar < tfar);
//    if (intersection)
//    {
//        unsigned int object = ray.hit.primID;
//        Vector3f orientation = (m_surfacePoints[m_surfaceTriangles[object](1)]-point).normalized();
//        Vector3f normal = (m_surfaceNormals[m_surfaceTriangles[object](0)]+m_surfaceNormals[m_surfaceTriangles[object](1)]+m_surfaceNormals[m_surfaceTriangles[object](2)])/3.f;
//        //Vector3f normal = Vector3f(ray.hit.Ng_x, ray.hit.Ng_y, ray.hit.Ng_z);
//        normal.normalize();
//        intersection = (orientation.dot(normal)>0);
//    }
//    return intersection;
//}

//bool Surface::isInside(Vector3f point)
//{
//    float tfar = 500;
//    return (launchRayInside(1,0,0,tfar, point) && launchRayInside(0,1,0,tfar, point));
//}

//bool Surface::launchRay(Vector3f direction, float tfar, Vector3f origin, Vector3f &intersection)
//{
//    RTCRayHit ray;
//    ray.ray.org_x = origin(0);
//    ray.ray.org_y = origin(1);
//    ray.ray.org_z = origin(2);
//    ray.ray.dir_x = direction(0);
//    ray.ray.dir_y = direction(1);
//    ray.ray.dir_z = direction(2);
//    ray.ray.tnear = 0.0f;
//    ray.ray.tfar = tfar;
//    ray.ray.flags = 0;
//    ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
//    ray.hit.geomID = RTC_INVALID_GEOMETRY_ID;
//    rtcInitIntersectContext(&m_rtccontext);
//    rtcIntersect1(m_rtcscene, &m_rtccontext, &ray);
//    bool intersectionFound = (ray.ray.tfar < tfar);
//    if (intersectionFound)
//    {
//        unsigned int object = ray.hit.primID;
//        //intersection = m_surfacePoints[m_surfaceTriangles[object](0)]*(1-ray.hit.u-ray.hit.v)+m_surfacePoints[m_surfaceTriangles[object](1)]*ray.hit.u+m_surfacePoints[m_surfaceTriangles[object](2)]*ray.hit.v;
//        intersection = origin + direction*ray.ray.tfar;
//        //intersection = origin;
//        //tfar = ray.ray.tfar;
//    }
//    return intersectionFound;
//}

//bool Surface::intersect(Vector3f origin, Vector3f direction/*its length determines the ending of intersection detection*/, Vector3f &intersection)
//{
//    float tfar = direction.norm();
//    direction.normalize();
//    //cout << "tfar before intersection : " << tfar << endl;
//    return launchRay(direction, tfar, origin, intersection);;
//}
