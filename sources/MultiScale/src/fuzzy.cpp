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
            dataFile >> m_affineTransfo(i,j);
    voxelDimensionsInCm = /*Vector3f::Ones()**/m_affineTransfo.block(0,0,3,3).diagonal();
    m_affineTransfoInverse=m_affineTransfo.inverse();
    cout << "Dimensions of the grid: " << m_dimX << "x" << m_dimY << "x" << m_dimZ << endl;
    cout << "Affine transformation: " << endl << m_affineTransfo << endl << "Voxels dimensions in cm: " << endl << voxelDimensionsInCm << endl;

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
    m_affineTransfo.resize(0,0);
    m_affineTransfoInverse.resize(0,0);
    m_data.clear();
    m_contourPts.clear();
    m_distance.clear();
}

unsigned int Fuzzy::getValueFromPoint(Vector3f p) const
{
    Vector3f gridPosition = p*m_affineTransfoInverse.block(0,0,3,3)+m_affineTransfoInverse.block(0,3,3,1);
    return m_data[(int)floor(gridPosition(2))+m_dimZ*(int)floor(gridPosition(1))+m_dimY*m_dimZ*(int)floor(gridPosition(0))];
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

float Fuzzy::getValueFromSegment(Vector3f p1, Vector3f p2)const
{
    Vector3f direction = p2-p1;
    Vector3f directionN = direction.normalized();
    float totalDist = direction.norm();
    //Number of voxels changes along x axis
    unsigned int voxelsX = static_cast<unsigned>(fabs(floor(p2.x()/voxelDimensionsInCm.x())-floor(p1.x()/voxelDimensionsInCm.x())));
    unsigned int voxelsY = static_cast<unsigned>(fabs(floor(p2.y()/voxelDimensionsInCm.y())-floor(p1.y()/voxelDimensionsInCm.y())));
    unsigned int voxelsZ = static_cast<unsigned>(fabs(floor(p2.z()/voxelDimensionsInCm.z())-floor(p1.z()/voxelDimensionsInCm.z())));
    unsigned int totalVoxelsCrossed = 1 + voxelsX + voxelsY + voxelsZ;
    Vector3f pTemp = p1;
    float value = 0;
    float epsilon = 0.0001f;
    float gammaX, gammaY, gammaZ, gamma;
    float sommeGamma = 0;
    for (unsigned int i=0; i<totalVoxelsCrossed - 1; i++)
    {
        if (fabs(directionN.x())<epsilon)
            gammaX =  numeric_limits<float>::max();
        else
            gammaX = (k(pTemp.x(), direction.x(), voxelDimensionsInCm.x()) - pTemp.x()) / directionN.x();
        if (fabs(directionN.y())<epsilon)
            gammaY =  numeric_limits<float>::max();
        else
            gammaY = (k(pTemp.y(), direction.y(), voxelDimensionsInCm.y()) - pTemp.y()) / directionN.y();
        if (fabs(directionN.z())<epsilon)
            gammaZ =  numeric_limits<float>::max();
        else
            gammaZ = (k(pTemp.z(), direction.z(), voxelDimensionsInCm.z()) - pTemp.z()) / directionN.z();
        if (gammaX == 0) gammaX =  numeric_limits<float>::max();
        if (gammaY == 0) gammaY =  numeric_limits<float>::max();
        if (gammaZ == 0) gammaZ =  numeric_limits<float>::max();

        gamma = min(min(gammaX, gammaY), gammaZ);
        sommeGamma += gamma;
        pTemp += (gamma+epsilon) * directionN;
        value += gamma/totalDist * getValueFromPoint(pTemp);
    }
    value += (p2-pTemp).norm()/totalDist * getValueFromPoint(p2);
    sommeGamma += (p2-pTemp).norm();

    return value;
}

float Fuzzy::getValueFromVolume(Vector3f pointBefore, Transfo T1, Transfo T2, Vector3f pointAfter) const
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
            Vector3f ellipse1 = Vector3f(r*T1.a * cos(p*2*M_PI), r*T1.b * sin(p*2*M_PI), 0.0f);
            Vector3f ellipse2 = Vector3f(r*T2.a * cos(p*2*M_PI), r*T2.b * sin(p*2*M_PI), 0.0f);
            Vector3f tangentCurve1 = (T2.center - pointBefore).normalized();
            Vector3f tangentCurve2 = (pointAfter - T1.center).normalized();
            Vector3f axeY1 = tangentCurve1.cross(T1.ellipseOrientation);
            Vector3f axeY2 = tangentCurve2.cross(T2.ellipseOrientation);
            Matrix3f rotation1;
            rotation1 << T1.ellipseOrientation.transpose(), axeY1.transpose(), tangentCurve1.transpose();
            Matrix3f rotation2;
            rotation2 << T2.ellipseOrientation.transpose(), axeY2.transpose(), tangentCurve2.transpose();
            Vector3f p1 = rotation1 * ellipse1 + T1.center;
            Vector3f p2 = rotation2 * ellipse2 + T2.center;
            value += getValueFromSegment(p1, p2);
        }
    }
    value /= (nberPerEllipse * nberOfRadii + 1);
    return value;
}

float Fuzzy::getDistanceToMaskFromSurface(Transfo T, Vector3f tangent)
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
            Vector3f ellipse = Vector3f(r*T.a * cos(p*2*M_PI), r*T.b * sin(p*2*M_PI), 0.0f);
            Vector3f tangentCurve = tangent.normalized();
            Vector3f axeY = tangentCurve.cross(T.ellipseOrientation);
            Matrix3f rotation;
            rotation << T.ellipseOrientation.transpose(), axeY.transpose(), tangentCurve.transpose();
            Vector3f pt = rotation * ellipse + T.center;
            value+=getDistanceToMask(pt);
        }
    }
    value /= (nberPerEllipse * nberOfRadii + 1);
    return value;
}

float Fuzzy::getDistanceToMask(Vector3f pt)
{
    if (getValueFromPoint(pt)==1)
        return 0;
    Vector3f gridPosition = pt*m_affineTransfoInverse.block(0,0,3,3)+m_affineTransfoInverse.block(0,3,3,1);
    unsigned int position=(int)floor(gridPosition(2))+m_dimZ*(int)floor(gridPosition(1))+m_dimY*m_dimZ*(int)floor(gridPosition(0));
    if (m_distance[position]!=-1)
        return m_distance[position];
    float distance=numeric_limits<float>::max();
    float temp;
    for (unsigned int i=0; i<m_contourPts.size(); i++)
    {
        temp=(m_contourPts[i]-pt).norm();
        if (temp<distance)
            distance=temp;
    }
    m_distance[position]=distance;
    return distance;
}

void Fuzzy::computeDistanceField()
{
    //Contour
    vector<vector<Vector3f>> contours(omp_get_max_threads());
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
                        contours[omp_get_thread_num()].push_back(Vector3f(i,j,k)*m_affineTransfo.block(0,0,3,3)+m_affineTransfo.block(0,3,3,1));
                    else
                    {
                        for (unsigned int l=i-1; l<i+2; l++)
                            for (unsigned int m=j-1;m<j+2;m++)
                                for (unsigned int n=k-1;n<k+2;n++)
                                    value+=m_data[n+m*m_dimZ+l*m_dimZ*m_dimY];
                        if (value!=27)
                            contours[omp_get_thread_num()].push_back(Vector3f(i,j,k)*m_affineTransfo.block(0,0,3,3)+m_affineTransfo.block(0,3,3,1));
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
