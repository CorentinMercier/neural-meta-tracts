#include "fiber_helper.h"

vector<Vector3f> build_circular_profile(float radius,int N_sample)
{
    vector<Vector3f> profile(N_sample);
    for(int k=0; k<N_sample; ++k)
    {
        const float u = k/static_cast<float>(N_sample);

        const float x = radius*cos(u*2*M_PI);
        const float y = radius*sin(u*2*M_PI);
        const float z = 0.0f;

        profile[k]=(Vector3f(x,y,z));
    }
    return profile;
}

vector<Vector3f> getOrientation(Fiber &fib)
{
    vector<Vector3f> orientation(fib.size());
    for (unsigned int i=0; i<fib.size(); i++)
        orientation[i] = fib.getProfileTransform(i).ellipseOrientation;
    return orientation;
}

vector<vector<Vector3f> > build_elliptical_profile(Fiber &fib,int N_sample)
{
    vector<vector<Vector3f> > profile(fib.size());
    for (unsigned int i=0; i<fib.size(); i++)
    {
        float a=fib.getProfileTransform(i).a;
        float b=fib.getProfileTransform(i).b;
        profile[i].resize(N_sample);
        for(int k=0; k<N_sample; ++k)
        {
            const float u = k/static_cast<float>(N_sample);
            const float x = a*cos(u*2*M_PI);
            const float y = b*sin(u*2*M_PI);
            const float z = 0.0f;
            profile[i][k]=Vector3f(x,y,z);
        }
    }
    return profile;
}

vector<Matrix3f> getRotations(Fiber &fib)
{
    vector<Matrix3f> rotations(fib.size());
    for (unsigned int i=0; i<fib.size(); i++)
        rotations[i] = fib.getTransform(i);
    return rotations;
}

vector<Vector3f> build_profile(Fiber &fib, int N_sample)
{
    vector<Vector3f> profile(N_sample);
    float a=fib.getProfile().ellipse.a;
    float b=fib.getProfile().ellipse.b;
    for(int k=0; k<N_sample; ++k)
    {
        const float u = k/static_cast<float>(N_sample);
        const float x = a*cos(u*2*M_PI);
        const float y = b*sin(u*2*M_PI);
        const float z = 0.0f;
        profile[k]=Vector3f(x,y,z);
    }
    return profile;
}

vector<Vector3f> merge_fiber_into_contiguous_vertices_lines(const vector<vector<Vector3f> > &fibers)
{
    vector<Vector3f> lines;
    for( size_t k_fiber=0; k_fiber<fibers.size(); ++k_fiber )
        for( size_t k_sample=0; k_sample<fibers[k_fiber].size()-1; ++k_sample )
        {
            lines.push_back(fibers[k_fiber][k_sample]);
            lines.push_back(fibers[k_fiber][k_sample+1]);
        }
    return lines;
}

