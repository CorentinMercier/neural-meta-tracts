#include "fiber_helper.h"

vector<glm::vec3> build_circular_profile(float radius,int N_sample)
{
	vector<glm::vec3> profile(N_sample);
	for(int k=0; k<N_sample; ++k)
	{
		const float u = k/static_cast<float>(N_sample);

		const float x = radius*cos(u*2*M_PI);
		const float y = radius*sin(u*2*M_PI);
		const float z = 0.0f;

		profile[k]=(glm::vec3(x,y,z));
	}
	return profile;
}

vector<glm::vec3> getOrientation(Fiber &fib)
{
	vector<glm::vec3> orientation(fib.size());
	for (unsigned int i=0; i<fib.size(); i++)
		orientation[i] = fib.getProfileTransform(i).ellipseOrientation;
	return orientation;
}

vector<vector<glm::vec3> > build_elliptical_profile(Fiber &fib,int N_sample)
{
	vector<vector<glm::vec3> > profile(fib.size());
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
			profile[i][k]=glm::vec3(x,y,z);
		}
	}
	return profile;
}

vector<glm::mat3> getRotations(Fiber &fib)
{
	vector<glm::mat3> rotations(fib.size());
	for (unsigned int i=0; i<fib.size(); i++)
		rotations[i] = fib.getTransform(i);
	return rotations;
}

vector<glm::vec3> build_profile(Fiber &fib, int N_sample)
{
	vector<glm::vec3> profile(N_sample);
	float a=fib.getProfile().ellipse.a;
	float b=fib.getProfile().ellipse.b;
	for(int k=0; k<N_sample; ++k)
	{
		const float u = k/static_cast<float>(N_sample);
		const float x = a*cos(u*2*M_PI);
		const float y = b*sin(u*2*M_PI);
		const float z = 0.0f;
		profile[k]=glm::vec3(x,y,z);
	}
	return profile;
}

vector<glm::vec3> merge_fiber_into_contiguous_vertices_lines(const vector<vector<glm::vec3> > &fibers)
{
	vector<glm::vec3> lines;
	for( size_t k_fiber=0; k_fiber<fibers.size(); ++k_fiber )
		for( size_t k_sample=0; k_sample<fibers[k_fiber].size()-1; ++k_sample )
		{
			lines.push_back(fibers[k_fiber][k_sample]);
			lines.push_back(fibers[k_fiber][k_sample+1]);
		}
	return lines;
}

