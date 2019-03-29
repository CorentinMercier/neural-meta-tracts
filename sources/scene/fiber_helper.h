#pragma once

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include "../MultiScale/include/fiber.h"



using namespace std;

vector<Vector3f> build_circular_profile(float radius,int N_sample);
vector<vector<Vector3f> > build_elliptical_profile(Fiber &fib,int N_sample);
vector<Vector3f> build_profile(Fiber &fib, int N_sample);
vector<Vector3f> getOrientation(Fiber &fib);
vector<Matrix3f> getRotations(Fiber &fib);

vector<Vector3f> merge_fiber_into_contiguous_vertices_lines(vector<vector<Vector3f> > const& fibers);
