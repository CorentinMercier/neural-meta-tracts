#ifndef CAMERA_H
#define CAMERA_H

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class camera
{
public:
    camera();

    const Matrix4f getModelview() {return m_modelview;}
    const Matrix4f getProjection() {return m_projection;}
    const Matrix4f getNormal() {return m_normal;}

    void setModelview(Matrix4f modelview) {m_modelview=modelview;}
    void setProjection(Matrix4f projection) {m_projection=projection;}
    void getNormal(Matrix4f normal) {m_normal=normal;}

    void setupCamera();
    void computeProjection();
    void computeModelview();

    void setScreenSize(unsigned int width, unsigned int height){m_screenSizeX=width; m_screenSizeY=height;}
    void setDist(float dist){m_dist=dist;}
    void setQuaternion(Quaternion<float> q){m_quat=q;}
    void setTranslation(Vector3f t){m_translation=t;}

    unsigned int getScreenSizeX(){return m_screenSizeX;}
    unsigned int getScreenSizeY(){return m_screenSizeY;}
    float getDist(){return m_dist;}
    Quaternion<float> getQuat(){return m_quat;}
    Vector3f getTranslation(){return m_translation;}

    int &xPrevious(){return m_xPrevious;}
    int &yPrevious(){return m_yPrevious;}

    //Movements of camera
    void goUp(float dL);
    void goRight(float dL);
    void goForward(float dL);

    void rotation(const int x, const int y);
    void rotateAlongY(float x);
    void zoom(int const y);

private:

    float projectToDisc(float const x, float const y);

    //Matrices for 3D transformations
    Matrix4f m_modelview;
    Matrix4f m_projection;
    Matrix4f m_normal;

    //Camera's parameters
    float m_fov=55.0f*M_PI/180.0f;
    float m_znear=1e-1f;
    float m_zfar=500.0f;
    float m_dist=-10.0f;
    unsigned int m_screenSizeX=800;
    unsigned int m_screenSizeY=800;
    Vector3f m_translation=Vector3f(0.0f, 0.0f, -1.0f);

    //Rotations
    Quaternion<float> m_quat;
    float m_discRadius=0.8;

    //Mouse positions
    int m_xPrevious;
    int m_yPrevious;
};

#endif // CAMERA_H
