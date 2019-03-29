#include "../include/fiber.h"
//#include "../include/surface.h"

Fiber::Fiber(MatrixXf points, unsigned int size)
{
    //m_data.rows()=size;
    m_data.resizeLike(points);
    m_width.resize(m_data.rows());
    m_profileTransformation.resize(m_data.rows());
    m_transfoMatrix.resize(size);
    //    Transfo test;
    //    test.a=MIN_RADIUS;
    //    test.b=MIN_RADIUS;
    m_data.row(0)=points.row(0);
    m_width[0]=0.f;
    Vector3f orient=(points.row(1)-points.row(0)).normalized();
    //    {
    //        Vector3f orientation=(points.row(1)-points.row(0)).normalized();
    //        test.rotation=orientation;
    //        test.ellipseOrientation=Vector3f(-orientation[1], orientation[0], 0).normalized();
    //        test.center=m_data.row(0);
    //    }
    Transfo test = {.center=m_data.row(0), .rotation=orient, .ellipseOrientation=Vector3f(-orient[1], orient[0], 0).normalized(), .a=MIN_RADIUS, .b=MIN_RADIUS};
    m_profile.ellipse=test;
    m_profileTransformation[0]=test;
    for (unsigned int i=1; i<size-1; i++)
    {
        m_data.row(i)=points.row(i);
        m_width[i]=0.f;
        Vector3f orientation=(points.row(i+1)-points.row(i-1)).normalized();
        test.rotation=orientation;
        test.ellipseOrientation=Vector3f(-orientation[1], orientation[0], 0).normalized();
        test.center=m_data.row(i);
        m_profileTransformation[i]=test;
    }
    m_data.row(size-1)=points.row(size-1);
    m_width[size-1]=0.f;
    {
        Vector3f orientation=(points.row(size-1)-points.row(size-2)).normalized();
        test.rotation=orientation;
        test.ellipseOrientation=Vector3f(-orientation[1], orientation[0], 0).normalized();
        test.center=m_data.row(size-1);
    }
    m_profileTransformation[size-1]=test;

    Vector3f lastOrientation=m_profileTransformation[0].ellipseOrientation;
    for (unsigned int i=0; i<m_profileTransformation.size(); i++)
    {
        Matrix3f matrice;
        if (m_profileTransformation[i].ellipseOrientation.dot(lastOrientation)<0)
            lastOrientation=-m_profileTransformation[i].ellipseOrientation;
        else
            lastOrientation=m_profileTransformation[i].ellipseOrientation;
        Vector3f axeY=m_profileTransformation[i].rotation.cross(lastOrientation);
        matrice << lastOrientation(0), axeY(0), m_profileTransformation[i].rotation(0),
                lastOrientation(1), axeY(1), m_profileTransformation[i].rotation(1),
                lastOrientation(2), axeY(2), m_profileTransformation[i].rotation(2);
        m_transfoMatrix[i]=matrice;
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
    m_data.resize(0,0);
    m_centers.resize(0,0);
    m_hierarchy.clear();
    m_neighbours.clear();
    m_profileTransformation.clear();
    m_tangents.resize(0,0);
    m_transfoMatrix.clear();
    m_width.clear();
}

//Fiber& Fiber::operator=(Fiber f)
//{

//}

Vector3f Fiber::operator[](unsigned int element) const
{
    if (element<m_data.rows())
        if (m_flipped)//Si la fibre est retournée, on la parcourt dans l'autre sens
            return (Vector3f)m_data.row(m_data.rows()-element-1);
        else
            return (Vector3f)m_data.row(element);
    else
    {
        cerr << "Erreur, accès à un élément de la fibre non existant : " << element << "/" << m_data.rows() << endl;
        return (Vector3f)(0);
    }
}

float Fiber::getPoint(int i, int j)
{
    if (m_flipped)
        return m_data(m_data.rows()-1-i, j);
    else
        return m_data(i,j);
}

Vector3f Fiber::getCenter(int i)
{
    if (m_flipped)
        return m_centers.row(m_centers.rows()-1-i);
    else
        return m_centers.row(i);
}

Vector3f Fiber::getTangent(int i)
{
    if (m_flipped)
        return m_tangents.row(m_tangents.rows()-1-i);
    else
        return m_tangents.row(i);
}

void Fiber::setTangents(MatrixXf t)
{
    m_tangents.resizeLike(t);
    for (int i=0; i<t.rows(); i++)
        m_tangents.row(i)=t.row(i);
}

float Fiber::length()
{
    //No computing if already done
    if (m_length!=0) return m_length;
    m_length=0;
    for (unsigned int i=0; i<m_data.rows()-1; i++)
        m_length+=(m_data.row(i+1)-m_data.row(i)).norm();
    return m_length;
}

float Fiber::cylinderLength()
{
    //No computing if already done
    if (m_Clength!=0) return m_Clength;
    m_Clength=0;
    for (unsigned int i=0; i<m_data.rows()-1; i++)
        m_Clength+=(m_profileTransformation[i+1].center-m_profileTransformation[i].center).norm();
    return m_Clength;
}

vector<float> Fiber::cumulateLength()
{
    vector<float> cumulatedLength(m_data.rows());
    cumulatedLength[0]=0;
    for (unsigned int i=1; i<m_data.rows(); i++)
        cumulatedLength[i]=(m_data.row(i)-m_data.row(i-1)).norm()+cumulatedLength[i-1];
    m_length=cumulatedLength[m_data.rows()-1];
    return cumulatedLength;
}

void Fiber::resample(unsigned int nb_points)
{
    //    cout << "Affichage matrice" << endl;
    //    for (unsigned int i=0; i< m_data.rows(); i++)
    //    {
    //        cout << m_data.row(i) << endl;
    //    }
    //    cout << "Fin affichage matrice" << endl;
    //    cout << "resample" << endl;
    if (nb_points<2)
    {
        cerr << "Error, impossible to resample with less than 2 points" << endl;
        return;
    }
    //No resample if fiber already has the number of points asked to prevent computational approximations leading to errors
    if (nb_points==this->size())
        return;
    vector<float> cumulatedLength = this->cumulateLength();

    //cout << "resample2" << endl;

    MatrixXf new_points;
    //vector<Vector3f> newpoint;
    //Vector3f temp;
    float next_point=0;
    float ratio;
    Vector3f delta;
    unsigned int i=0, j=0, k=0;
    float step=cumulatedLength[m_data.rows()-1]/(nb_points-1);//(m_data.rows()-1);

    new_points.resize(nb_points,3);
    //newpoint.resize(nb_points);

    //cout << "resample3" << endl;
    while (next_point<m_length)
    {
        if (next_point==cumulatedLength[k])
        {
            new_points.row(i)=m_data.row(j);
            next_point+=step;
            i++; j++; k++;
        }
        else if (next_point<cumulatedLength[k])
        {
            ratio=1-(cumulatedLength[k]-next_point)/(cumulatedLength[k]-cumulatedLength[k-1]);
            delta=(m_data.row(j)-m_data.row(j-1));
            new_points.row(i)=(Vector3f)m_data.row(j)+(Vector3f)(ratio*delta);
            next_point+=step;
            i++;
        }
        else
        {
            j++; k++;
        }
    }
    //cout << "resample4" << endl;
    //Adding of the last point
    new_points.row(nb_points-1)=m_data.row(m_data.rows()-1);
    //temp = m_data.row(m_data.rows()-1);
    //newpoint[nb_points-1] = temp;
    //cout << "resample4" << endl;

    //Current fibers gets the data
    //m_data.rows()=nb_points;
    //delete(m_data);
    //m_data.resize(nb_points, 3);
    //    cout << new_points.size() << endl;
    //    cout << "Affichage matrice" << endl;
    //    for (unsigned int i=0; i< m_data.rows(); i++)
    //    {
    //        cout << m_data.row(i) << endl;
    //    }
    //    cout << "Fin affichage matrice" << endl;
    //m_data.resizeLike(new_points);
    //m_data = new_points;
    m_data.conservativeResize(nb_points, 3);
    //    cout << m_data.rows() << endl;
    //    cout << "resample4" << endl;
    m_width.resize(m_data.rows());
    //    cout << "resample4" << endl;
    m_profileTransformation.resize(m_data.rows());
    //Transfo test;
    //test.a=MIN_RADIUS;
    //test.b=MIN_RADIUS;
    //    cout << "resample5" << endl;
    m_data.row(0)=new_points.row(0);
    Vector3f orient=(new_points.row(1)-new_points.row(0)).normalized();
    m_width[0]=0.f;
    //    {
    //        Vector3f orientation=(new_points.row(1)-new_points.row(0)).normalized();
    //        test.rotation=orientation;
    //        test.ellipseOrientation=Vector3f(-orientation[1], orientation[0], 0).normalized();
    //        test.center=m_data.row(0);
    //    }
    Transfo test = {.center=m_data.row(0), .rotation=orient, .ellipseOrientation=Vector3f(-orient[1], orient[0], 0).normalized(), .a=MIN_RADIUS, .b=MIN_RADIUS};
    m_profile.ellipse=test;
    m_profileTransformation[0]=test;
    for (unsigned int i=1; i<m_data.rows()-1; i++)
    {
        m_data.row(i)=new_points.row(i);
        m_width[i]=0.f;
        Vector3f orientation=(new_points.row(i+1)-new_points.row(i-1)).normalized();
        test.rotation=orientation;
        test.ellipseOrientation=Vector3f(-orientation[1], orientation[0], 0).normalized();
        test.center=m_data.row(i);
        m_profileTransformation[i]=test;
    }
    //    cout << "resample6" << endl;
    m_data.row(m_data.rows()-1)=new_points.row(m_data.rows()-1);
    m_width[m_data.rows()-1]=0.f;
    {
        Vector3f orientation=(new_points.row(m_data.rows()-1)-new_points.row(m_data.rows()-2)).normalized();
        test.rotation=orientation;
        test.ellipseOrientation=Vector3f(-orientation[1], orientation[0], 0).normalized();
        test.center=m_data.row(m_data.rows()-1);
    }
    m_profileTransformation[m_data.rows()-1]=test;
    m_length=0;

    //    cout << "resample7" << endl;
    m_transfoMatrix.resize(m_data.rows());

    Vector3f lastOrientation=m_profileTransformation[0].ellipseOrientation;
    for (unsigned int i=0; i<m_profileTransformation.size(); i++)
    {
        Matrix3f matrice;
        if (m_profileTransformation[i].ellipseOrientation.dot(lastOrientation)<0)
            lastOrientation=-m_profileTransformation[i].ellipseOrientation;
        else
            lastOrientation=m_profileTransformation[i].ellipseOrientation;
        Vector3f axeY=m_profileTransformation[i].rotation.cross(lastOrientation);
        matrice << lastOrientation(0), axeY(0), m_profileTransformation[i].rotation(0),
                lastOrientation(1), axeY(1), m_profileTransformation[i].rotation(1),
                lastOrientation(2), axeY(2), m_profileTransformation[i].rotation(2);
        m_transfoMatrix[i]=matrice;
    }
    //    cout << "resample8" << endl;

    //    m_data.resizeLike(new_points);
    //    m_width.resize(nb_points);
    //    for (int i=0; i<nb_points; i++)
    //    {
    //        m_data.row(i)=new_points.row(i);
    //        m_width[i]=0.f;
    //    }
    //    m_data.rows()=nb_points;
    //    m_length=0;
}

void Fiber::resampleV2(unsigned int nb_points)
{
    if (nb_points<2)
    {
        cerr << "Error, impossible to resample with less than 2 points" << endl;
        return;
    }
    //No resample if fiber already has the number of points asked to prevent computational approximations leading to errors
    if (nb_points==this->size())
        return;
    vector<float> cumulatedLength = this->cumulateLength();

    MatrixXf new_points;
    float next_point=0;
    float ratio;
    Vector3f delta;
    unsigned int i=0, j=0, k=0;
    float step=cumulatedLength[m_data.rows()-1]/(nb_points-1);

    new_points.resize(nb_points,3);
    while (next_point<m_length)
    {
        //cout << i << endl;
        //cout << m_data.row(0) << endl << endl;
        if (next_point==cumulatedLength[j])
        {
            new_points.row(i)=m_data.row(j);
            //temp = m_data.row(j);
            //newpoint[i] = temp;//m_data.row(j);
            next_point+=step;
            i++; j++; //k++;
        }
        else if (next_point<cumulatedLength[j])
        {
            ratio=1-(cumulatedLength[j]-next_point)/(cumulatedLength[j]-cumulatedLength[k]);
            delta=(m_data.row(j)-m_data.row(k));
            new_points.row(i)=(Vector3f)m_data.row(j)+(Vector3f)(ratio*delta);
            //temp = (Vector3f)m_data.row(j)+(Vector3f)(ratio*delta);
            //newpoint[i] = temp;
            next_point+=step;
            i++;
            k=j-1;
        }
        else
        {
            //j++; k++;
            j++;
        }
    }
    //Adding of the last point
    new_points.row(nb_points-1)=m_data.row(m_data.rows()-1);

    m_data.conservativeResize(nb_points, 3);
    m_width.resize(m_data.rows());
    m_profileTransformation.resize(m_data.rows());
    m_data.row(0)=new_points.row(0);
    Vector3f orient=(new_points.row(1)-new_points.row(0)).normalized();
    m_width[0]=0.f;
    Transfo test = {.center=m_data.row(0), .rotation=orient, .ellipseOrientation=Vector3f(-orient[1], orient[0], 0).normalized(), .a=MIN_RADIUS, .b=MIN_RADIUS};
    m_profile.ellipse=test;
    m_profileTransformation[0]=test;
    for (unsigned int i=1; i<m_data.rows()-1; i++)
    {
        m_data.row(i)=new_points.row(i);
        m_width[i]=0.f;
        Vector3f orientation=(new_points.row(i+1)-new_points.row(i-1)).normalized();
        test.rotation=orientation;
        test.ellipseOrientation=Vector3f(-orientation[1], orientation[0], 0).normalized();
        test.center=m_data.row(i);
        m_profileTransformation[i]=test;
    }
    m_data.row(m_data.rows()-1)=new_points.row(m_data.rows()-1);
    m_width[m_data.rows()-1]=0.f;
    {
        Vector3f orientation=(new_points.row(m_data.rows()-1)-new_points.row(m_data.rows()-2)).normalized();
        test.rotation=orientation;
        test.ellipseOrientation=Vector3f(-orientation[1], orientation[0], 0).normalized();
        test.center=m_data.row(m_data.rows()-1);
    }
    m_profileTransformation[m_data.rows()-1]=test;
    m_length=0;

    m_transfoMatrix.resize(m_data.rows());

    Vector3f lastOrientation=m_profileTransformation[0].ellipseOrientation;
    for (unsigned int i=0; i<m_profileTransformation.size(); i++)
    {
        Matrix3f matrice;
        if (m_profileTransformation[i].ellipseOrientation.dot(lastOrientation)<0)
            lastOrientation=-m_profileTransformation[i].ellipseOrientation;
        else
            lastOrientation=m_profileTransformation[i].ellipseOrientation;
        Vector3f axeY=m_profileTransformation[i].rotation.cross(lastOrientation);
        matrice << lastOrientation(0), axeY(0), m_profileTransformation[i].rotation(0),
                lastOrientation(1), axeY(1), m_profileTransformation[i].rotation(1),
                lastOrientation(2), axeY(2), m_profileTransformation[i].rotation(2);
        m_transfoMatrix[i]=matrice;
    }
}

float Fiber::getWidth(unsigned int i)
{
    if (m_flipped)
        return m_width[m_data.rows()-1-i];
    else
        return m_width[i];
}


void Fiber::setWidth(vector<float> & width)
{
    m_width.resize(width.size());
    for (unsigned int i=0; i< width.size(); i++)
        m_width[i]=width[i];
}


Transfo Fiber::getProfileTransform(unsigned int i)
{
    if (m_flipped)
        return m_profileTransformation[m_data.rows()-1-i];
    else
        return m_profileTransformation[i];
}

void Fiber::setProfileTransform(vector<Transfo> &profileTransformation)
{
    m_profileTransformation.resize(profileTransformation.size());
    for (unsigned int i=0; i<profileTransformation.size(); i++)
    {
        m_profileTransformation[i].rotation=profileTransformation[i].rotation;
        m_profileTransformation[i].ellipseOrientation=profileTransformation[i].ellipseOrientation;
        m_profileTransformation[i].a=profileTransformation[i].a;
        m_profileTransformation[i].b=profileTransformation[i].b;
        m_profileTransformation[i].center=profileTransformation[i].center;
    }
}


void Fiber::PutInSameDirection(Fiber const &fiber)
{
    float distanceEndpoints1=0;
    float distanceEndpoints2=0;
    //    distanceEndpoints1=((*this)[0]-fiber[0]).squaredNorm()+((*this)[m_data.rows()-1]-fiber[fiber.size()-1]).squaredNorm();
    //    distanceEndpoints2=((*this)[0]-fiber[fiber.size()-1]).squaredNorm()+((*this)[m_data.rows()-1]-fiber[0]).squaredNorm();
    distanceEndpoints1=min(((*this)[0]-fiber[0]).squaredNorm(),((*this)[m_data.rows()-1]-fiber[fiber.size()-1]).squaredNorm());
    distanceEndpoints2=min(((*this)[0]-fiber[fiber.size()-1]).squaredNorm(),((*this)[m_data.rows()-1]-fiber[0]).squaredNorm());
    //    cout << "Points 0 : " << (*this)[0] << " " << fiber[0] << endl;
    //    cout << "Points de fin : " << (*this)[m_data.rows()-1] << " " << fiber[fiber.size()-1] << endl;
    //    cout << "Distances : " << distanceEndpoints1 << " " << distanceEndpoints2 << endl;
    //The shortest distance determines if the fibers are oriented the same way or not
    if (distanceEndpoints2<distanceEndpoints1)
        this->flip();
}

void Fiber::computeCentersAndTangents()
{
    if (m_computedCenterAndTangents)
        return;
    m_centers.resize(m_data.rows()-1, 3);
    m_tangents.resize(m_data.rows()-1, 3);
#pragma omp parallel for
    for (unsigned int i=0; i<m_data.rows()-1; i++)
    {
        m_centers.row(i)=(m_data.row(i)+m_data.row(i+1))/2.0;
        m_tangents.row(i)=m_data.row(i+1)-m_data.row(i);
    }
    m_computedCenterAndTangents=true;
}

//PutInSameDirection(fiber) should be called before
//Compute the max euclidian distance between the endpoints of the two fibers
float Fiber::getMaxDistEndPoints(Fiber &fiber)
{
    float distanceEndpoints1=((*this)[0]-fiber[0]).squaredNorm();
    float distanceEndpoints2=((*this)[m_data.rows()-1]-fiber[fiber.size()-1]).squaredNorm();
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
    if (m_volume!=-1)
        return m_volume;
    vector<float>theta(m_data.rows());
    theta[0]=M_PI_2;
    theta[m_data.rows()-1]=M_PI_2;
    for (unsigned int i=1; i<m_data.rows()-1; i++)
    {
        float alpha=acos(min(max(Vector3f((*this)[i]-(*this)[i-1]).normalized().dot(Vector3f((*this)[i+1]-(*this)[i]).normalized()), -1.f), 1.f));
        //        if (isnan(alpha))
        //            alpha=0;
        theta[i]=(M_PI-alpha)/2;
        if (isnan(theta[i]))
            cout << "Acos nan : " << Vector3f((*this)[i]-(*this)[i-1]).normalized().dot(Vector3f((*this)[i+1]-(*this)[i]).normalized()) << endl;
    }
    float cumulVolume=0;
    for (unsigned int i=0; i<m_data.rows()-1; i++)
    {
        float h1=sin(theta[i])*m_width[i];
        float h2=sin(theta[i+1]*m_width[i+1]);
        float L=cos(theta[i])*m_width[i]+cos(theta[i+1])*m_width[i+1]+Vector3f((*this)[i+1]-(*this)[i]).norm();
        cumulVolume+=pow(L,3)/3*pow(h2-h1,2) + pow(L,2)*h1*(h2-h1) + L*pow(h1,2);
    }
    cumulVolume*=M_PI;
    m_volume=cumulVolume;
    return m_volume;
}

float Fiber::getSelfScalarProduct(float lambdag)
{
    if (m_selfScalarProductComputed)
        return m_selfScalarProduct;
    computeCentersAndTangents();
    float scalarProduct=0;
    float kg;
#pragma omp parallel for reduction (+:scalarProduct)
    for (unsigned int i=0; i<m_data.rows()-1; i++)
        for (unsigned int j=0; j<m_data.rows()-1; j++)
        {
            kg=exp(-(this->getCenter(i)-this->getCenter(j)).squaredNorm()/(lambdag*lambdag));
            scalarProduct+=this->getTangent(i).transpose()*kg*this->getTangent(j);
        }
    m_selfScalarProduct=scalarProduct;
    m_selfScalarProductComputed=true;
    return m_selfScalarProduct;
}

float Fiber::getSelfVarifoldScalarProduct(float sigma)
{
    if (m_selfScalarProductComputed)
        return m_selfScalarProduct;
    this->computeCentersAndTangents();
    float scalarProduct=0;
    float ksigma=0, kn=0;
#pragma omp parallel for reduction (+:scalarProduct)
    for (unsigned int i=0; i<m_data.rows()-1; i++)
        for (unsigned int j=0; j<m_data.rows()-1; j++)
        {
            ksigma=exp(-(this->getCenter(i)-this->getCenter(j)).squaredNorm()/(sigma*sigma));
            kn=(this->getTangent(i).transpose()*this->getTangent(j));
            scalarProduct+=ksigma*kn*kn/(this->getTangent(i).norm()*this->getTangent(j).norm());
        }
    m_selfScalarProduct=scalarProduct;
    m_selfScalarProductComputed=true;
    return m_selfScalarProduct;
}

float smooth(float center, float left, float right)
{
    return (center+(left+right)/2.f)/2.f;
}

Transfo smooth(Transfo t1, Transfo t2, Transfo t3)
{
    t1.a=smooth(t1.a, t2.a, t3.a);
    t1.b=smooth(t1.b, t2.b, t3.b);
    t1.ellipseOrientation=(t1.ellipseOrientation+(t2.ellipseOrientation+t3.ellipseOrientation)/2)/2;
    t1.rotation=(t1.rotation+(t2.rotation+t3.rotation)/2)/2;
    return t1;
}

void Fiber::geometricSmoothing()
{
    for (unsigned int i=0; i<size()-3; i++)
    {
        m_profileTransformation[i+1]=smooth(m_profileTransformation[i+1], m_profileTransformation[i], m_profileTransformation[i+2]);
    }
}

void Fiber::pointsSmoothing()
{
    for (unsigned int i=0; i<size()-3; i++)
    {
        m_data.row(i+1)=(m_data.row(i+1)+(m_data.row(i)+m_data.row(i+2))/2)/2;
    }
}

void Fiber::moveExtremities(Vector3f deviation)
{
    m_data.row(0)+=deviation;
    m_data.row(m_data.rows()-1)+=deviation;
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
//    Vector3f endPoint;
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
//        //            endPoint = m_data.row(ending);
//        //        }
//    }
//    //    endPoint = m_data.row(ending-2);
//    //    ending-=3;
//    Vector3f beginPoint;
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
//        //            beginPoint = m_data.row(0);
//        //        }
//    }
//    //    cout << beginning << " " << ending << endl;
//    //    beginPoint = m_data.row(beginning);
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
//    data.row(0) = beginPoint;
//    data.row(ending-beginning+2) = endPoint;
//    //    cout << "Nouvelle fibre" << endl;
//    unsigned int nbInside = 0;
//    for (unsigned int i=0; i<data.rows(); i++)
//        if (surface->isInside(data.row(i)))
//            nbInside++;
//    if ((float)nbInside/data.rows()<0.5)
//    {
//        Fiber noFib;
//        return noFib;
//    }
//    //        cout << i << endl << data.row(i) << endl << endl;


//    //    MatrixXf data = m_data.block(beginning, 0, (ending-beginning), 3);
//    Fiber newFib(data, data.rows());
//    return newFib;
//}
