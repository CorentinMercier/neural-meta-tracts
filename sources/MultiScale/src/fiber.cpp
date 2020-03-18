#include "../include/fiber.h"
//#include "../include/surface.h"

Fiber::Fiber(std::vector<glm::vec3> points, unsigned int size)
{
	m_data = std::move(points);
	m_width.resize(m_data.size());
	m_profileTransformation.resize(m_data.size());
	m_transfoMatrix.resize(size);
	//    Transfo test;
	//    test.a=MIN_RADIUS;
	//    test.b=MIN_RADIUS;
	m_width[0]=0.f;
	glm::vec3 orient = glm::normalize(m_data[1] - m_data[0]);
	//    {
	//        glm::vec3 orientation=(points[1)-points[0)).normalized();
	//        test.rotation=orientation;
	//        test.ellipseOrientation=glm::vec3(-orientation[1], orientation[0], 0).normalized();
	//        test.center=m_data[0);
	//    }
	Transfo test = {.center=m_data[0], .rotation=orient, .ellipseOrientation = glm::normalize(glm::vec3(-orient[1], orient[0], 0)), .a=MIN_RADIUS, .b=MIN_RADIUS};
	m_profile.ellipse = test;
	m_profileTransformation[0] = test;
	for (unsigned int i=1; i<size-1; i++)
	{
		m_width[i] = 0.f;
		glm::vec3 orientation = glm::normalize(m_data[i+1] - m_data[i-1]);
		test.rotation = orientation;
		test.ellipseOrientation = glm::normalize(glm::vec3(-orientation[1], orientation[0], 0));
		test.center = m_data[i];
		m_profileTransformation[i] = test;
	}
	m_width[size-1] = 0.f;
	{
		glm::vec3 orientation = glm::normalize(m_data[size - 1] - m_data[size - 2]);
		test.rotation = orientation;
		test.ellipseOrientation = glm::normalize(glm::vec3(-orientation[1], orientation[0], 0));
		test.center = m_data[size - 1];
	}
	m_profileTransformation[size - 1] = test;

	glm::vec3 lastOrientation = m_profileTransformation[0].ellipseOrientation;
	for (unsigned int i=0; i<m_profileTransformation.size(); i++)
	{
		glm::mat3 matrice;
		if (glm::dot(m_profileTransformation[i].ellipseOrientation, lastOrientation) < 0)
			lastOrientation = -m_profileTransformation[i].ellipseOrientation;
		else
			lastOrientation = m_profileTransformation[i].ellipseOrientation;
		glm::vec3 axeY = glm::cross(m_profileTransformation[i].rotation, lastOrientation);
		matrice = glm::transpose(glm::mat3(lastOrientation[0], axeY[0], m_profileTransformation[i].rotation[0],
				lastOrientation[1], axeY[1], m_profileTransformation[i].rotation[1],
				lastOrientation[2], axeY[2], m_profileTransformation[i].rotation[2]));
		m_transfoMatrix[i] = matrice;
	}
	//Compute a random color for the fiber
	m_red=(float)rand()/RAND_MAX*255;
	m_green=(float)rand()/RAND_MAX*255;
	m_blue=(float)rand()/RAND_MAX*255;
}

Fiber::Fiber()
{
	m_valid = false;
}

Fiber::~Fiber()
{
	m_data.clear();
	m_centers.clear();
	m_hierarchy.clear();
	m_neighbours.clear();
	m_profileTransformation.clear();
	m_tangents.clear();
	m_transfoMatrix.clear();
	m_width.clear();
}

glm::vec3 Fiber::operator[](unsigned int element) const
{
	if (element<m_data.size())
		if (m_flipped)//If fiber is flipped, points are run through the opposite way
			return m_data[m_data.size() - element - 1];
		else
			return m_data[element];
	else
	{
		cerr << "Erreur, accès à un élément de la fibre non existant : " << element << "/" << m_data.size() << endl;
		return (glm::vec3)(0);
	}
}

float Fiber::getPoint(int i, int j)
{
	if (m_flipped)
		return m_data[m_data.size() - 1 - i][j];
	else
		return m_data[i][j];
}

glm::vec3 Fiber::getCenter(int i)
{
	if (m_flipped)
		return m_centers[m_centers.size() - 1 - i];
	else
		return m_centers[i];
}

glm::vec3 Fiber::getTangent(int i)
{
	if (m_flipped)
		return m_tangents[m_tangents.size() - 1 - i];
	else
		return m_tangents[i];
}

void Fiber::setTangents(std::vector<glm::vec3> t)
{
	m_tangents = std::move(t);
}

float Fiber::length()
{
	//No computing if already done
	if (m_length != 0) return m_length;
	m_length = 0;
	for (unsigned int i=0; i<m_data.size()-1; i++)
		m_length += glm::length(m_data[i + 1] - m_data[i]);
	return m_length;
}

float Fiber::cylinderLength()
{
	//No computing if already done
	if (m_Clength != 0) return m_Clength;
	m_Clength = 0;
	for (unsigned int i=0; i<m_data.size()-1; i++)
		m_Clength += glm::length(m_profileTransformation[i + 1].center - m_profileTransformation[i].center);
	return m_Clength;
}

vector<float> Fiber::cumulateLength()
{
	vector<float> cumulatedLength(m_data.size());
	cumulatedLength[0] = 0;
	for (unsigned int i=1; i<m_data.size(); i++)
		cumulatedLength[i] = glm::length(m_data[i] - m_data[i - 1]) + cumulatedLength[i - 1];
	m_length = cumulatedLength[m_data.size() - 1];
	return cumulatedLength;
}

void Fiber::resample(unsigned int nb_points)
{
	if (nb_points<2)
	{
		cerr << "Error, impossible to resample with less than 2 points" << endl;
		return;
	}
	//No resample if fiber already has the number of points asked to prevent computational approximations leading to errors
	if (nb_points == this->size())
		return;
	vector<float> cumulatedLength = this->cumulateLength();


	std::vector<glm::vec3> new_points;
	float next_point = 0;
	float ratio;
	glm::vec3 delta;
	unsigned int i = 0, j = 0, k = 0;
	float step = cumulatedLength[m_data.size() - 1] / (nb_points - 1);

	new_points.resize(nb_points);

	while (next_point < m_length)
	{
		if (next_point == cumulatedLength[k])
		{
			new_points[i] = m_data[j];
			next_point += step;
			i++; j++; k++;
		}
		else if (next_point < cumulatedLength[k])
		{
			ratio = 1 - (cumulatedLength[k] - next_point) / (cumulatedLength[k] - cumulatedLength[k - 1]);
			delta = (m_data[j] - m_data[j - 1]);
			new_points[i] = m_data[j] + (glm::vec3)(ratio * delta);
			next_point += step;
			i++;
		}
		else
		{
			j++; k++;
		}
	}
	//Adding of the last point
	new_points[nb_points - 1] = m_data[m_data.size() - 1];
	m_data.resize(nb_points);
	m_width.resize(m_data.size());
	m_profileTransformation.resize(m_data.size());
	m_data[0] = new_points[0];
	glm::vec3 orient = glm::normalize(new_points[1] - new_points[0]);
	m_width[0] = 0.f;
	Transfo test = {.center=m_data[0], .rotation=orient, .ellipseOrientation=glm::normalize(glm::vec3(-orient[1], orient[0], 0)), .a=MIN_RADIUS, .b=MIN_RADIUS};
	m_profile.ellipse = test;
	m_profileTransformation[0] = test;
	for (unsigned int i=1; i<m_data.size()-1; i++)
	{
		m_data[i] = new_points[i];
		m_width[i] = 0.f;
		glm::vec3 orientation = glm::normalize(new_points[i + 1] - new_points[i - 1]);
		test.rotation = orientation;
		test.ellipseOrientation = glm::normalize(glm::vec3(-orientation[1], orientation[0], 0));
		test.center = m_data[i];
		m_profileTransformation[i] = test;
	}
	m_data[m_data.size() - 1] = new_points[m_data.size() - 1];
	m_width[m_data.size() - 1] = 0.f;
	{
		glm::vec3 orientation = glm::normalize(new_points[m_data.size() - 1] - new_points[m_data.size() - 2]);
		test.rotation = orientation;
		test.ellipseOrientation = glm::normalize(glm::vec3(-orientation[1], orientation[0], 0));
		test.center = m_data[m_data.size() - 1];
	}
	m_profileTransformation[m_data.size() - 1] = test;
	m_length = 0;

	m_transfoMatrix.resize(m_data.size());

	glm::vec3 lastOrientation = m_profileTransformation[0].ellipseOrientation;
	for (unsigned int i=0; i<m_profileTransformation.size(); i++)
	{
		glm::mat3 matrice;
		if (glm::dot(m_profileTransformation[i].ellipseOrientation, lastOrientation) < 0)
			lastOrientation = -m_profileTransformation[i].ellipseOrientation;
		else
			lastOrientation = m_profileTransformation[i].ellipseOrientation;
		glm::vec3 axeY = glm::cross(m_profileTransformation[i].rotation, lastOrientation);
		matrice = glm::transpose(glm::mat3(lastOrientation[0], axeY[0], m_profileTransformation[i].rotation[0],
				lastOrientation[1], axeY[1], m_profileTransformation[i].rotation[1],
				lastOrientation[2], axeY[2], m_profileTransformation[i].rotation[2]));
		m_transfoMatrix[i] = matrice;
	}
}

void Fiber::resampleV2(unsigned int nb_points)
{
	if (nb_points<2)
	{
		cerr << "Error, impossible to resample with less than 2 points" << endl;
		return;
	}
	//No resample if fiber already has the number of points asked to prevent computational approximations leading to errors
	if (nb_points == this->size())
		return;
	vector<float> cumulatedLength = this->cumulateLength();

	std::vector<glm::vec3> new_points;
	float next_point = 0;
	float ratio;
	glm::vec3 delta;
	unsigned int i = 0, j = 0, k = 0;
	float step = cumulatedLength[m_data.size() - 1] / (nb_points - 1);

	new_points.resize(nb_points);
	while (next_point<m_length)
	{
		if (next_point == cumulatedLength[j])
		{
			new_points[i] = m_data[j];
			next_point += step;
			i++; j++;
		}
		else if (next_point < cumulatedLength[j])
		{
			ratio = 1 - (cumulatedLength[j] - next_point) / (cumulatedLength[j] - cumulatedLength[k]);
			delta = (m_data[j] - m_data[k]);
			new_points[i] = m_data[j] + (glm::vec3)(ratio * delta);
			next_point += step;
			i++;
			k = j - 1;
		}
		else
			j++;
	}
	//Adding of the last point
	new_points[nb_points - 1] = m_data[m_data.size() - 1];

	m_data.resize(nb_points);
	m_width.resize(m_data.size());
	m_profileTransformation.resize(m_data.size());
	m_data[0] = new_points[0];
	glm::vec3 orient = glm::normalize(new_points[1] - new_points[0]);
	m_width[0] = 0.f;
	Transfo test = {.center=m_data[0], .rotation=orient, .ellipseOrientation=glm::normalize(glm::vec3(-orient[1], orient[0], 0)), .a=MIN_RADIUS, .b=MIN_RADIUS};
	m_profile.ellipse = test;
	m_profileTransformation[0] = test;
	for (unsigned int i=1; i<m_data.size() - 1; i++)
	{
		m_data[i] = new_points[i];
		m_width[i] = 0.f;
		glm::vec3 orientation = glm::normalize(new_points[i + 1] - new_points[i - 1]);
		test.rotation = orientation;
		test.ellipseOrientation = glm::normalize(glm::vec3(-orientation[1], orientation[0], 0));
		test.center = m_data[i];
		m_profileTransformation[i] = test;
	}
	m_data[m_data.size() - 1] = new_points[m_data.size() - 1];
	m_width[m_data.size() - 1] = 0.f;
	{
		glm::vec3 orientation = glm::normalize(new_points[m_data.size() - 1] - new_points[m_data.size() - 2]);
		test.rotation = orientation;
		test.ellipseOrientation = glm::normalize(glm::vec3(-orientation[1], orientation[0], 0));
		test.center = m_data[m_data.size() - 1];
	}
	m_profileTransformation[m_data.size() - 1] = test;
	m_length = 0;

	m_transfoMatrix.resize(m_data.size());

	glm::vec3 lastOrientation = m_profileTransformation[0].ellipseOrientation;
	for (unsigned int i=0; i<m_profileTransformation.size(); i++)
	{
		glm::mat3 matrice;
		if (glm::dot(m_profileTransformation[i].ellipseOrientation, lastOrientation) < 0)
			lastOrientation = -m_profileTransformation[i].ellipseOrientation;
		else
			lastOrientation = m_profileTransformation[i].ellipseOrientation;
		glm::vec3 axeY = glm::cross(m_profileTransformation[i].rotation, lastOrientation);
		matrice = glm::transpose(glm::mat3(lastOrientation[0], axeY[0], m_profileTransformation[i].rotation[0],
				lastOrientation[1], axeY[1], m_profileTransformation[i].rotation[1],
				lastOrientation[2], axeY[2], m_profileTransformation[i].rotation[2]));
		m_transfoMatrix[i] = matrice;
	}
}

float Fiber::getWidth(unsigned int i)
{
	if (m_flipped)
		return m_width[m_data.size() - 1 - i];
	else
		return m_width[i];
}


void Fiber::setWidth(vector<float> & width)
{
	m_width.resize(width.size());
#pragma omp parallel for
	for (unsigned int i=0; i< width.size(); i++)
		m_width[i] = width[i];
}


Transfo Fiber::getProfileTransform(unsigned int i)
{
	if (m_flipped)
		return m_profileTransformation[m_data.size() - 1 - i];
	else
		return m_profileTransformation[i];
}

void Fiber::setProfileTransform(vector<Transfo> &profileTransformation)
{
	m_profileTransformation.resize(profileTransformation.size());
#pragma omp parallel for
	for (unsigned int i=0; i<profileTransformation.size(); i++)
	{
		m_profileTransformation[i].rotation = profileTransformation[i].rotation;
		m_profileTransformation[i].ellipseOrientation = profileTransformation[i].ellipseOrientation;
		m_profileTransformation[i].a = profileTransformation[i].a;
		m_profileTransformation[i].b = profileTransformation[i].b;
		m_profileTransformation[i].center = profileTransformation[i].center;
	}
}


void Fiber::PutInSameDirection(Fiber const &fiber)
{
	float distanceEndpoints1 = 0;
	float distanceEndpoints2 = 0;
	//    distanceEndpoints1=((*this)[0]-fiber[0]).squaredNorm()+((*this)[m_data.size()-1]-fiber[fiber.size()-1]).squaredNorm();
	//    distanceEndpoints2=((*this)[0]-fiber[fiber.size()-1]).squaredNorm()+((*this)[m_data.size()-1]-fiber[0]).squaredNorm();
	distanceEndpoints1 = min(glm::length2((*this)[0] - fiber[0]), glm::length2((*this)[m_data.size() - 1] - fiber[fiber.size() - 1]));
	distanceEndpoints2 = min(glm::length2((*this)[0] - fiber[fiber.size() - 1]), glm::length2((*this)[m_data.size()-1]-fiber[0]));
	//    cout << "Points 0 : " << (*this)[0] << " " << fiber[0] << endl;
	//    cout << "Points de fin : " << (*this)[m_data.size()-1] << " " << fiber[fiber.size()-1] << endl;
	//    cout << "Distances : " << distanceEndpoints1 << " " << distanceEndpoints2 << endl;
	//The shortest distance determines if the fibers are oriented the same way or not
	if (distanceEndpoints2 < distanceEndpoints1)
		this->flip();
}

void Fiber::computeCentersAndTangents()
{
	if (m_computedCenterAndTangents)
		return;
	m_centers.resize(m_data.size() - 1);
	m_tangents.resize(m_data.size() - 1);
#pragma omp parallel for
	for (unsigned int i=0; i<m_data.size() - 1; i++)
	{
		m_centers[i] = (m_data[i] + m_data[i + 1]) / 2.f;
		m_tangents[i] = m_data[i + 1] - m_data[i];
	}
	m_computedCenterAndTangents = true;
}

//PutInSameDirection(fiber) should be called before
//Compute the max euclidian distance between the endpoints of the two fibers
float Fiber::getMaxDistEndPoints(Fiber &fiber)
{
	float distanceEndpoints1 = glm::length2((*this)[0]-fiber[0]);
	float distanceEndpoints2 = glm::length2((*this)[m_data.size()-1]-fiber[fiber.size()-1]);
	return (max(distanceEndpoints1, distanceEndpoints2));
}

void Fiber::addNeighbours(vector<unsigned int> new_neighbours)
{
	for (unsigned int i=0; i<new_neighbours.size(); i++)
		m_neighbours.insert(new_neighbours[i]);
}

void Fiber::addToHierarchy(vector<unsigned int> new_elements)
{
	for (unsigned int i=0; i<new_elements.size(); i++)
		m_hierarchy.insert(new_elements[i]);
}

vector<unsigned int> Fiber::getHierarchy()
{
	vector<unsigned int> hierarchy;
	for (set<unsigned int>::iterator ite=m_hierarchy.begin(); ite!=m_hierarchy.end(); ++ite)
		hierarchy.push_back((unsigned int)*ite);
	return hierarchy;
}

float Fiber::computeVolume()
{
	if (m_volume != -1)
		return m_volume;
	vector<float>theta(m_data.size());
	theta[0] = M_PI_2;
	theta[m_data.size() - 1] = M_PI_2;
	for (unsigned int i=1; i<m_data.size() - 1; i++)
	{
		float alpha = acos(min(max(glm::dot(glm::normalize((*this)[i] - (*this)[i - 1]), glm::normalize((*this)[i + 1] - (*this)[i])), -1.f), 1.f));
		theta[i] = (M_PI - alpha) / 2.f;
		if (isnan(theta[i]))
			cout << "Acos nan : " << glm::dot(glm::normalize((*this)[i] - (*this)[i - 1]), glm::normalize((*this)[i + 1] - (*this)[i])) << endl;
	}
	float cumulVolume = 0;
	for (unsigned int i=0; i<m_data.size() - 1; i++)
	{
		float h1 = sin(theta[i]) * m_width[i];
		float h2 = sin(theta[i + 1] * m_width[i + 1]);
		float L = cos(theta[i]) * m_width[i] + cos(theta[i + 1]) * m_width[i + 1] + glm::length((*this)[i + 1] - (*this)[i]);
		cumulVolume += pow(L, 3) / 3 * pow(h2 - h1, 2) + pow(L, 2) * h1 * (h2 - h1) + L * pow(h1, 2);
	}
	cumulVolume *= M_PI;
	m_volume = cumulVolume;
	return m_volume;
}

float Fiber::getSelfScalarProduct(float lambdag)
{
	if (m_selfScalarProductComputed)
		return m_selfScalarProduct;
	computeCentersAndTangents();
	float scalarProduct = 0;
	float kg;
#pragma omp parallel for reduction (+:scalarProduct)
	for (unsigned int i=0; i<m_data.size()-1; i++)
		for (unsigned int j=0; j<m_data.size()-1; j++)
		{
			kg = exp(-glm::length2(this->getCenter(i) - this->getCenter(j)) / (lambdag * lambdag));
			scalarProduct += glm::dot(this->getTangent(i), (kg * this->getTangent(j)));
		}
	m_selfScalarProduct = scalarProduct;
	m_selfScalarProductComputed = true;
	return m_selfScalarProduct;
}

float Fiber::getSelfVarifoldScalarProduct(float sigma)
{
	if (m_selfScalarProductComputed)
		return m_selfScalarProduct;
	this->computeCentersAndTangents();
	float scalarProduct = 0;
	float ksigma = 0, kn = 0;
#pragma omp parallel for reduction (+:scalarProduct)
	for (unsigned int i=0; i<m_data.size() - 1; i++)
		for (unsigned int j=0; j<m_data.size() - 1; j++)
		{
			ksigma = exp(-glm::length2(this->getCenter(i) - this->getCenter(j)) / (sigma * sigma));
			kn = glm::dot(this->getTangent(i), this->getTangent(j));
			scalarProduct += ksigma * kn * kn / (glm::length(this->getTangent(i)) * glm::length(this->getTangent(j)));
		}
	m_selfScalarProduct = scalarProduct;
	m_selfScalarProductComputed = true;
	return m_selfScalarProduct;
}

float smooth(float center, float left, float right)
{
	return (center + (left + right) / 2.f) / 2.f;
}

Transfo smooth(Transfo t1, Transfo t2, Transfo t3)
{
	t1.a = smooth(t1.a, t2.a, t3.a);
	t1.b = smooth(t1.b, t2.b, t3.b);
	t1.ellipseOrientation = (t1.ellipseOrientation+(t2.ellipseOrientation+t3.ellipseOrientation) / 2.f) / 2.f;
	t1.rotation = (t1.rotation + (t2.rotation + t3.rotation) / 2.f) / 2.f;
	return t1;
}

void Fiber::geometricSmoothing()
{
	for (unsigned int i=0; i<size() - 3; i++)
	{
		m_profileTransformation[i + 1] = smooth(m_profileTransformation[i + 1], m_profileTransformation[i], m_profileTransformation[i + 2]);
	}
}

void Fiber::pointsSmoothing()
{
	for (unsigned int i=0; i<size()-3; i++)
	{
		m_data[i + 1] = (m_data[i+1] + (m_data[i] + m_data[i + 2]) / 2.f) / 2.f;
	}
}

void Fiber::moveExtremities(glm::vec3 deviation)
{
	m_data[0] += deviation;
	m_data[m_data.size() - 1] += deviation;
}

//Fiber Fiber::cut(Surface *surface)
//{
//    float lengthOfLastRay = 10;
//    //Start at the middle of the fiber
//    unsigned int startingPoint = size()/2;
//    if (!surface->isInside((*this)[startingPoint]))
//    {
//        Fiber noFib;
//        return noFib;
//    }
//    unsigned int beginning = startingPoint;
//    unsigned int ending = startingPoint;
//    //Version starting at the middle point and launching a ray between two consecutive points, checking for the absence/presence of the surface between the points
//    glm::vec3 endPoint;
//    while (ending<size()-1 && !surface->intersect((*this)[ending], (*this)[ending+1]-(*this)[ending], endPoint))
//        ending++;
//    //if (ending == size()-1 && !surface->intersect((*this)[ending], (*this)[ending+1]-(*this)[ending], endPoint))
//    if (ending == size()-1)// && !surface->intersect((*this)[ending-1], (*this)[ending]-(*this)[ending-1], endPoint))
//    {
//        if (!surface->intersect((*this)[ending], ((*this)[ending]-(*this)[ending-1])*lengthOfLastRay, endPoint))
//        {
//            Fiber noFib;
//            return noFib;
//        }
//        //        else
//        //        {
//        //            //surface->intersect((*this)[ending], ((*this)[ending]-(*this)[ending-1])*100, endPoint);
//        //            endPoint = m_data[ending);
//        //        }
//    }
//    //    endPoint = m_data[ending-2);
//    //    ending-=3;
//    glm::vec3 beginPoint;
//    while (beginning!=0 && !surface->intersect((*this)[beginning], (*this)[beginning-1]-(*this)[beginning], beginPoint))
//        beginning--;
//    //if (beginning == 1 && !surface->intersect((*this)[beginning], (*this)[beginning-1]-(*this)[beginning], beginPoint))
//    if (beginning == 0)// && !surface->intersect((*this)[beginning+1], (*this)[beginning]-(*this)[beginning+1], beginPoint))
//    {
//        if (!surface->intersect((*this)[beginning], ((*this)[beginning]-(*this)[beginning+1])*lengthOfLastRay, beginPoint))
//        {
//            Fiber noFib;
//            return noFib;
//        }
//        //        else
//        //        {
//        //            beginPoint = m_data[0);
//        //        }
//    }
//    //    cout << beginning << " " << ending << endl;
//    //    beginPoint = m_data[beginning);
//    //    beginning++;

//    //Version checking that every consecutive point is inside, starting from the middle one
//    //    while(ending<size() && surface->isInside((*this)[ending]))
//    //        ending++;
//    //    while(beginning!=0 && surface->isInside((*this)[beginning]))
//    //        beginning--;


//    if (static_cast<int>(ending)-static_cast<int>(beginning)<2)
//    {
//        Fiber noFib;
//        return noFib;
//    }

//    MatrixXf data;
//    data.resize(ending-beginning+3, 3);
//    data.block(1, 0, (ending-beginning+1), 3) = m_data.block(beginning, 0, (ending-beginning+1), 3);
//    data[0) = beginPoint;
//    data[ending-beginning+2) = endPoint;
//    //    cout << "Nouvelle fibre" << endl;
//    unsigned int nbInside = 0;
//    for (unsigned int i=0; i<data.size(); i++)
//        if (surface->isInside(data[i)))
//            nbInside++;
//    if ((float)nbInside/data.size()<0.5)
//    {
//        Fiber noFib;
//        return noFib;
//    }
//    //        cout << i << endl << data[i) << endl << endl;


//    //    MatrixXf data = m_data.block(beginning, 0, (ending-beginning), 3);
//    Fiber newFib(data, data.size());
//    return newFib;
//}
