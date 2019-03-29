#include "camera.h"

camera::camera()
{
//    Quaternion<float> q1=AngleAxisf(-M_PI_2, Vector3f::UnitX());
//    Quaternion<float> q2=AngleAxisf(M_PI_2, Vector3f::UnitY());
//    m_quat=q2*q1;
    m_quat.setIdentity();
}

void camera::setupCamera()
{
    computeProjection();
    computeModelview();
}

void camera::computeProjection()
{
    float ratio=static_cast<float>((float)m_screenSizeX/m_screenSizeY);
    const float f=1.0f/tan(m_fov/2.0f);
    const float fp=f/ratio;

    const float L=m_znear-m_zfar;
    //const float epsilon=1e-6f;

    //ASSERT(abs(L)>epsilon,"z_far-z_near too small");

    const float C=(m_zfar+m_znear)/L;
    const float D=(2.0f*m_zfar*m_znear)/L;

    m_projection << fp, 0,  0, 0,
                     0, f,  0, 0,
                     0, 0,  C, D,
                     0, 0, -1, 0;
}

void camera::computeModelview()
{
    Matrix4f worldMatrixZoom;
    worldMatrixZoom << 1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, m_dist,
                       0, 0, 0, 1;
    Matrix4f worldMatrixTranslation;
    worldMatrixTranslation << 1, 0, 0, m_translation[0],
                              0, 1, 0, m_translation[1],
                              0, 0, 1, m_translation[2],
                              0, 0, 0, 1;
    Matrix4f worldMatrixRotation=Matrix4f::Identity();
    Matrix3f m;
    m=m_quat.toRotationMatrix();//AngleAxisf(-M_PI_2, Vector3f::UnitX())*AngleAxisf(M_PI_2, Vector3f::UnitY());
    worldMatrixRotation.block<3,3>(0,0)=m;
    m_modelview=worldMatrixZoom*worldMatrixRotation*worldMatrixTranslation;
}


void camera::goUp(float dL)
{
    Vector3f const y(0,-1,0);
    m_translation += dL*(m_quat.conjugate()*y);
}

void camera::goRight(float dL)
{
    Vector3f const x(-1,0,0);
    m_translation += dL*(m_quat.conjugate()*x);
}

void camera::goForward(float dL)
{
    Vector3f const z(0,0,1);
    m_translation += dL*(m_quat.conjugate()*z);
}

float camera::projectToDisc(float const x, float const y)
{
    float const n=sqrt(x*x+y*y);
    if(n<m_discRadius*0.707107f)
        return sqrt(pow(m_discRadius,2)-n*n);
    else
        return pow(m_discRadius,2)/(2.0f*n);
}

void camera::rotation(int const x, int const y)
{
    float const xOld=m_xPrevious;
    float const yOld=m_yPrevious;

    float const x0=(2.0f*xOld-m_screenSizeX)/m_screenSizeX;
    float const y0=(m_screenSizeY-2.0f*yOld)/m_screenSizeY;
    float const x1=(2.0f*x-m_screenSizeX)/m_screenSizeX;
    float const y1=(m_screenSizeY-2.0f*y)/m_screenSizeY;
    if(sqrt(pow(x0-x1,2)+pow(y0-y1,2)>1e-6))
    {
        Vector3f const p1=Vector3f(x0, y0, projectToDisc(x0, y0));
        Vector3f const p2=Vector3f(x1, y1, projectToDisc(x1, y1));
        Vector3f const axis=(p1.cross(p2)).normalized();
        Vector3f const u=p1-p2;
        float t=u.norm()/(2.0f*m_discRadius);
        t=min(max(t,-1.0f),1.0f); //clamp
        float const phi = 2.0f*asin(t);
        //compute quaternion to apply
        m_quat=AngleAxisf(phi, axis)*m_quat;
    }
    m_xPrevious=x;
    m_yPrevious=y;
}

void camera::rotateAlongY(float x)
{
    float const phi = 2.0f*asin(x);
    //compute quaternion to apply
    m_quat=AngleAxisf(phi, Vector3f::UnitY())*m_quat;
}

void camera::zoom(int const y)
{
    m_dist += (fabs(m_dist)+1.0f)*(y-m_yPrevious)/500.0f;
    m_dist = min(m_dist,0.0f);
    m_yPrevious=y;
}
