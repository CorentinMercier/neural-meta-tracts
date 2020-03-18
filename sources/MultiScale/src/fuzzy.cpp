#include "../include/fuzzy.h"

Fuzzy::Fuzzy(string data)
{
	m_filenameData = data;
	fstream dataFile;
	dataFile.open(m_filenameData, ios_base::in);
	dataFile >> m_dimX;
	dataFile >> m_dimY;
	dataFile >> m_dimZ;
	for (unsigned int i=0; i<4; i++)
		for (unsigned int j=0; j<4; j++)
			dataFile >> m_affineTransfo[j][i];
	voxelDimensionsInCm = glm::vec3(m_affineTransfo[0][0], m_affineTransfo[1][1], m_affineTransfo[2][2]);
	m_affineTransfoInverse = glm::inverse(m_affineTransfo);
	cout << "Dimensions of the grid: " << m_dimX << "x" << m_dimY << "x" << m_dimZ << endl;
	cout << "Affine transformation: " << endl << glm::to_string(m_affineTransfo) << endl << "Voxels dimensions in cm: " << endl << glm::to_string(voxelDimensionsInCm) << endl;

	m_data.reserve(m_dimX*m_dimY*m_dimZ);
	m_distance.reserve(m_dimX*m_dimY*m_dimZ);
	for (unsigned int i=0; i<m_dimX; i++)
		for (unsigned int j=0; j<m_dimY; j++)
			for (unsigned int k=0; k<m_dimZ; k++)
			{
				//float a;
				int a;
				dataFile >> a;
				m_data.push_back((int)a);
				m_distance.push_back(-1);
			}
}

Fuzzy::~Fuzzy()
{
	m_data.clear();
	m_contourPts.clear();
	m_distance.clear();
}

unsigned int Fuzzy::getValueFromPoint(glm::vec3 p) const
{
	glm::vec3 gridPosition = p * glm::mat3(m_affineTransfoInverse) + glm::vec3(m_affineTransfoInverse[3][0], m_affineTransfoInverse[3][1] , m_affineTransfoInverse[3][2]);
	return m_data[(int)floor(gridPosition[2]) + m_dimZ * (int)floor(gridPosition[1]) + m_dimY * m_dimZ * (int)floor(gridPosition[0])];
}

unsigned int getSign(float value)
{
	if (value>0)
		return 1;
	else
		return 0;
}

float k(float x, float dirX, float dimV)
{
	return (getSign(dirX) + floor(x/dimV)) * dimV;
}

float Fuzzy::getValueFromSegment(glm::vec3 p1, glm::vec3 p2)const
{
	glm::vec3 direction = p2-p1;
	glm::vec3 directionN = glm::normalize(direction);
	float totalDist = glm::length(direction);
	//Number of voxels changes along x axis
	unsigned int voxelsX = static_cast<unsigned>(fabs(floor(p2.x / voxelDimensionsInCm.x) - floor(p1.x / voxelDimensionsInCm.x)));
	unsigned int voxelsY = static_cast<unsigned>(fabs(floor(p2.y / voxelDimensionsInCm.y) - floor(p1.y / voxelDimensionsInCm.y)));
	unsigned int voxelsZ = static_cast<unsigned>(fabs(floor(p2.z / voxelDimensionsInCm.z) - floor(p1.z / voxelDimensionsInCm.z)));
	unsigned int totalVoxelsCrossed = 1 + voxelsX + voxelsY + voxelsZ;
	glm::vec3 pTemp = p1;
	float value = 0;
	float epsilon = 0.0001f;
	float gammaX, gammaY, gammaZ, gamma;
	float sommeGamma = 0;
	for (unsigned int i=0; i<totalVoxelsCrossed - 1; i++)
	{
		if (fabs(directionN.x) < epsilon)
			gammaX =  numeric_limits<float>::max();
		else
			gammaX = (k(pTemp.x, direction.x, voxelDimensionsInCm.x) - pTemp.x) / directionN.x;
		if (fabs(directionN.y) < epsilon)
			gammaY =  numeric_limits<float>::max();
		else
			gammaY = (k(pTemp.y, direction.y, voxelDimensionsInCm.y) - pTemp.y) / directionN.y;
		if (fabs(directionN.z) < epsilon)
			gammaZ =  numeric_limits<float>::max();
		else
			gammaZ = (k(pTemp.z, direction.z, voxelDimensionsInCm.z) - pTemp.z) / directionN.z;
		if (gammaX == 0) gammaX =  numeric_limits<float>::max();
		if (gammaY == 0) gammaY =  numeric_limits<float>::max();
		if (gammaZ == 0) gammaZ =  numeric_limits<float>::max();

		gamma = min(min(gammaX, gammaY), gammaZ);
		sommeGamma += gamma;
		pTemp += (gamma+epsilon) * directionN;
		value += gamma/totalDist * getValueFromPoint(pTemp);
	}
	value += glm::length(p2 - pTemp) / totalDist * getValueFromPoint(p2);
	sommeGamma += glm::length(p2-pTemp);

	return value;
}

float Fuzzy::getValueFromVolume(glm::vec3 pointBefore, Transfo T1, Transfo T2, glm::vec3 pointAfter) const
{
	float value = 0;
	//Get center fiber value
	value += getValueFromSegment(T1.center, T2.center);
	unsigned int nberPerEllipse = 6;
	unsigned int nberOfRadii = 2;
	for (unsigned int j=1; j<nberOfRadii+1; j++)
	{
		float r=(float)j/nberOfRadii;
		for (unsigned int i=0; i<nberPerEllipse; i++)
		{
			float p = (float)i/nberPerEllipse;
			glm::vec3 ellipse1 = glm::vec3(r*T1.a * cos(p*2*M_PI), r*T1.b * sin(p*2*M_PI), 0.0f);
			glm::vec3 ellipse2 = glm::vec3(r*T2.a * cos(p*2*M_PI), r*T2.b * sin(p*2*M_PI), 0.0f);
			glm::vec3 tangentCurve1 = glm::normalize(T2.center - pointBefore);
			glm::vec3 tangentCurve2 = glm::normalize(pointAfter - T1.center);
			glm::vec3 axeY1 = glm::cross(tangentCurve1, T1.ellipseOrientation);
			glm::vec3 axeY2 = glm::cross(tangentCurve2, T2.ellipseOrientation);
			glm::mat3 rotation1 = glm::mat3(T1.ellipseOrientation, axeY1, tangentCurve1);
			glm::mat3 rotation2 = glm::mat3(T2.ellipseOrientation, axeY2, tangentCurve2);
			glm::vec3 p1 = rotation1 * ellipse1 + T1.center;
			glm::vec3 p2 = rotation2 * ellipse2 + T2.center;
			value += getValueFromSegment(p1, p2);
		}
	}
	value /= (nberPerEllipse * nberOfRadii + 1);
	return value;
}

float Fuzzy::getDistanceToMaskFromSurface(Transfo T, glm::vec3 tangent)
{
	float value=0;
	//Get center point value
	value+=getDistanceToMask(T.center);
	unsigned int nberPerEllipse = 6;
	unsigned int nberOfRadii = 2;
	for (unsigned int j=1; j<nberOfRadii+1; j++)
	{
		float r=(float)j/nberOfRadii;
		for (unsigned int i=0; i<nberPerEllipse; i++)
		{
			float p = (float)i/nberPerEllipse;
			glm::vec3 ellipse = glm::vec3(r*T.a * cos(p*2*M_PI), r*T.b * sin(p*2*M_PI), 0.0f);
			glm::vec3 tangentCurve = glm::normalize(tangent);
			glm::vec3 axeY = glm::cross(tangentCurve, T.ellipseOrientation);
			glm::mat3 rotation = glm::mat3(T.ellipseOrientation, axeY, tangentCurve);
			glm::vec3 pt = rotation * ellipse + T.center;
			value += getDistanceToMask(pt);
		}
	}
	value /= (nberPerEllipse * nberOfRadii + 1);
	return value;
}

float Fuzzy::getDistanceToMask(glm::vec3 pt)
{
	if (getValueFromPoint(pt)==1)
		return 0;
	glm::vec3 gridPosition = pt * glm::mat3(m_affineTransfoInverse) + glm::vec3(m_affineTransfoInverse[3][0], m_affineTransfoInverse[3][1] , m_affineTransfoInverse[3][2]);
	unsigned int position = (int)floor(gridPosition[2]) + m_dimZ * (int)floor(gridPosition[1]) + m_dimY * m_dimZ * (int)floor(gridPosition[0]);
	if (m_distance[position]!=-1)
		return m_distance[position];
	float distance=numeric_limits<float>::max();
	float temp;
	for (unsigned int i=0; i<m_contourPts.size(); i++)
	{
		temp = glm::length(m_contourPts[i] - pt);
		if (temp < distance)
			distance = temp;
	}
	m_distance[position] = distance;
	return distance;
}

void Fuzzy::computeDistanceField()
{
	//Contour
	vector<vector<glm::vec3>> contours(omp_get_max_threads());
	m_contourPts.clear();
#pragma omp parallel for
	for (unsigned int i=0; i<m_dimX; i++)
	{
		for (unsigned int j=0; j<m_dimY; j++)
		{
			for (unsigned int k=0; k<m_dimZ; k++)
			{
				unsigned int value = 0;
				if (m_data[k+j*m_dimZ+i*m_dimZ*m_dimY]==1)
				{
					if (i==0 || j==0 || k==0 || i==m_dimX-1 || j==m_dimY-1 || k==m_dimZ-1)
						contours[omp_get_thread_num()].push_back(glm::vec3(i,j,k) * glm::mat3(m_affineTransfo) + glm::vec3(m_affineTransfo[3][0], m_affineTransfo[3][1] , m_affineTransfo[3][2]));
					else
					{
						for (unsigned int l=i-1; l<i+2; l++)
							for (unsigned int m=j-1;m<j+2;m++)
								for (unsigned int n=k-1;n<k+2;n++)
									value+=m_data[n+m*m_dimZ+l*m_dimZ*m_dimY];
						if (value!=27)
							contours[omp_get_thread_num()].push_back(glm::vec3(i,j,k) * glm::mat3(m_affineTransfo) + glm::vec3(m_affineTransfo[3][0], m_affineTransfo[3][1] , m_affineTransfo[3][2]));
					}
				}
			}
		}
	}
	unsigned int sizeOfContour=0;
	for (int i=0; i<omp_get_max_threads(); i++)
		sizeOfContour+=contours[i].size();
	m_contourPts.reserve(sizeOfContour);
	for (int i=0; i<omp_get_max_threads(); i++)
	{
		for (unsigned int j=0; j<contours[i].size(); j++)
			m_contourPts.push_back(contours[i][j]);
	}
	cout << "Nber of contours: " << sizeOfContour << endl;
}
