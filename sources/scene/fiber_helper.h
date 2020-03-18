#pragma once

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include "../MultiScale/include/fiber.h"



using namespace std;

vector<glm::vec3> build_circular_profile(float radius,int N_sample);
vector<vector<glm::vec3> > build_elliptical_profile(Fiber &fib,int N_sample);
vector<glm::vec3> build_profile(Fiber &fib, int N_sample);
vector<glm::vec3> getOrientation(Fiber &fib);
vector<glm::mat3> getRotations(Fiber &fib);

vector<glm::vec3> merge_fiber_into_contiguous_vertices_lines(vector<vector<glm::vec3> > const& fibers);
