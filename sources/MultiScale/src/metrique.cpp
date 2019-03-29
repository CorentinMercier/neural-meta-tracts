#include "../include/metrique.h"
#include "../include/fiber.h"

Metrique::Metrique(unsigned int metric):m_metric(metric)
{
}

Metrique::~Metrique()
{}

float fiberLength(const MatrixXf &fiber, int size)
{
    float length=0;
    for (int i=0; i<size-1; i++)
        length+=(fiber.row(i+1)-fiber.row(i)).squaredNorm();
    return length;
}

float minimum(const VectorXf &tableau, int size)
{
    float mini=tableau(0);
    for (int i=1; i<size; i++)
        if (mini>tableau(i)) mini=tableau(i);
    if (mini<0) cerr << "Négatif : " << mini << endl;
    return mini;
}

float Metrique::MC(Fiber &fiber1, Fiber &fiber2)
{
    MatrixXf distances;
    float dm12=0;
    float dm21=0;
    float length1;
    float length2;

    length1=fiber1.length();
    length2=fiber2.length();

    distances.resize(fiber1.size(), fiber2.size());
    for (unsigned int i=0; i<fiber1.size(); i++)
        for (unsigned int j=0; j<fiber2.size(); j++)
        {
            distances(i,j)= (fiber1[i]-fiber2[j]).squaredNorm();
            if (distances(i,j)<0) cerr << "Négatif : " << distances(i,j) << endl;

        }

    for (unsigned int i=0; i<fiber1.size(); i++)
        dm12+=minimum(distances.row(i), fiber2.size())/length1;
    for (unsigned int i=0; i<fiber2.size(); i++)
        dm21+=minimum(distances.col(i), fiber1.size())/length2;

    distances.resize(0,0);

    return (dm12+dm21)/2;
}

float Metrique::averageOfPointwiseEuclideanMetric(Fiber &fiber1, Fiber &fiber2)
{
    if (fiber1.size()!=fiber2.size())
    {
        if (fiber1.size()>fiber2.size())
            fiber2.resample(fiber1.size());
        else
            fiber1.resample(fiber2.size());
    }
    float sum=0;
    for (unsigned int i=0; i<fiber1.size(); i++)
        sum+=(fiber1[i]-fiber2[i]).squaredNorm();
    sum/=fiber1.size();
    return sum;
}

//Shorter mean of closest distances
float Metrique::SC(Fiber &fiber1, Fiber &fiber2)
{
    MatrixXf distances;
    float dm12=0;
    float dm21=0;
    float length1;
    float length2;

    length1=fiber1.length();
    length2=fiber2.length();

    distances.resize(fiber1.size(), fiber2.size());
    for (unsigned int i=0; i<fiber1.size(); i++)
        for (unsigned int j=0; j<fiber2.size(); j++)
        {
            distances(i,j)= (fiber1[i]-fiber2[j]).squaredNorm();
            if (distances(i,j)<0) cerr << "Négatif : " << distances(i,j) << endl;
        }

    for (unsigned int i=0; i<fiber1.size(); i++)
        dm12+=minimum(distances.row(i), fiber2.size())/length1;
    for (unsigned int i=0; i<fiber2.size(); i++)
        dm21+=minimum(distances.col(i), fiber1.size())/length2;

    distances.resize(0,0);

    return min(dm12,dm21);
}

//Longer mean of closest distances
float Metrique::LC(Fiber &fiber1, Fiber &fiber2)
{
    MatrixXf distances;
    float dm12=0;
    float dm21=0;
    float length1;
    float length2;

    length1=fiber1.length();
    length2=fiber2.length();

    distances.resize(fiber1.size(), fiber2.size());
    for (unsigned int i=0; i<fiber1.size(); i++)
        for (unsigned int j=0; j<fiber2.size(); j++)
        {
            distances(i,j)= (fiber1[i]-fiber2[j]).squaredNorm();
            if (distances(i,j)<0) cerr << "Négatif : " << distances(i,j) << endl;
        }

    for (unsigned int i=0; i<fiber1.size(); i++)
        dm12+=minimum(distances.row(i), fiber2.size())/length1;
    for (unsigned int i=0; i<fiber2.size(); i++)
        dm21+=minimum(distances.col(i), fiber1.size())/length2;

    distances.resize(0,0);

    return max(dm12,dm21);
}

//MDF distance
float Metrique::MDF(Fiber &fiber1, Fiber &fiber2)
{
    if (fiber1.size()!=fiber2.size())
    {
        cerr << "Error : fibers must be the same size, use the fiber::resample(unsigned int nb_points) function" << endl;
        return -1;
    }
    float ddirect=0;
    float dflipped=0;
    unsigned int size = fiber1.size();

    for (unsigned int i=0; i<size; i++)
    {
        ddirect+=(fiber1[i]-fiber2[i]).squaredNorm();
        dflipped+=(fiber1[i]-fiber2[size-i]).squaredNorm();
    }
    ddirect/=(float)size;
    dflipped/=(float)size;

    return min(ddirect, dflipped);
}

float scalarProductPDM(Fiber &fiber1, Fiber &fiber2, float sigma)
{
    float scalarProduct=0;
    float denominator=fiber1.size()*fiber2.size();
    for (unsigned int i=0; i<fiber1.size(); i++)
        for (unsigned int j=0; j<fiber2.size(); j++)
            scalarProduct+=exp(-(fiber1[i]-fiber2[j]).squaredNorm()/(sigma*sigma))/denominator;
    return scalarProduct;
}

//Point Density Model
float Metrique::PDM(Fiber &fiber1, Fiber &fiber2, float sigma)
{
    float scalarProduct11=scalarProductPDM(fiber1, fiber1, sigma);
    float scalarProduct22=scalarProductPDM(fiber2, fiber2, sigma);
    float scalarProduct12=scalarProductPDM(fiber1, fiber2, sigma);
    return sqrt(scalarProduct11+scalarProduct22-2*scalarProduct12);
}

float scalarProductVD(Fiber &fiber1, Fiber &fiber2, float sigma)
{
    float scalarProduct=0;
    float ksigma=0, kn=0;
    for (unsigned int i=0; i<fiber1.size()-1; i++)
    {
        for (unsigned int j=0; j<fiber2.size()-1; j++)
        {
            ksigma=exp(-(fiber1.getCenter(i)-fiber2.getCenter(j)).squaredNorm()/(sigma*sigma));
            kn=(fiber1.getTangent(i).transpose()*fiber2.getTangent(j));
            scalarProduct+=ksigma*kn*kn/(fiber1.getTangent(i).norm()*fiber2.getTangent(j).norm());
//            if (isnan(scalarProduct) && fiber1.getTangent(i).norm()*fiber2.getTangent(j).norm()==0)
//            {
//#pragma omp critical
//                {
//                cout << "Scalar product nan : " << fiber1.getTangent(i) << " " <<fiber2.getTangent(j) << endl;
//                cout << i << " " << j << " " << fiber1.size() << " " << fiber2.size() << endl;
//                }
//            }
        }
    }
    return scalarProduct;
}

//Varifold distance
float Metrique::VarifoldDistance(Fiber &fiber1, Fiber &fiber2, float sigma)
{
    fiber1.computeCentersAndTangents();
    fiber2.computeCentersAndTangents();
    float scalarProduct11=fiber1.getSelfVarifoldScalarProduct(sigma);//scalarProductVD(fiber1, fiber1, sigma);//
    float scalarProduct22=fiber2.getSelfVarifoldScalarProduct(sigma);//scalarProductVD(fiber2, fiber2, sigma);//
    float scalarProduct12=scalarProductVD(fiber1, fiber2, sigma);
    return sqrt(scalarProduct11+scalarProduct22-2*scalarProduct12);
}

float scalarProductCurrent(Fiber &fiber1, Fiber &fiber2, float sigma)
{
    float scalarProduct=0;
    float ksigma=0, kn=0;
    for (unsigned int i=0; i<fiber1.size()-1; i++)
    {
        for (unsigned int j=0; j<fiber2.size()-1; j++)
        {
            ksigma=exp(-(fiber1.getCenter(i)-fiber2.getCenter(j)).squaredNorm()/(sigma*sigma));
            kn=fiber1.getTangent(i).dot(fiber2.getTangent(j));
            scalarProduct+=ksigma*kn;
        }
    }
    return scalarProduct;
}

//Current distance
float Metrique::CurrentDistance(Fiber &fiber1, Fiber &fiber2, float sigma)
{
    fiber1.computeCentersAndTangents();
    fiber2.computeCentersAndTangents();
    float scalarProduct11=scalarProductCurrent(fiber1, fiber1, sigma);
    float scalarProduct22=scalarProductCurrent(fiber2, fiber2, sigma);
    float scalarProduct12=scalarProductCurrent(fiber1, fiber2, sigma);
    return sqrt(scalarProduct11+scalarProduct22-2*scalarProduct12);
}

float scalarProductWeightedCurrents(Fiber &fiber1, Fiber &fiber2, float lambdac, float lambdab, float lambdag)
{
    float scalarProduct=0;
    float kc, kb, kg;
    for (unsigned int i=0; i<fiber1.size()-1; i++)
    {
        for (unsigned int j=0; j<fiber2.size()-1; j++)
        {
            kg=exp(-(fiber1.getCenter(i)-fiber2.getCenter(j)).squaredNorm()/(lambdag*lambdag));
            scalarProduct+=fiber1.getTangent(i).transpose()*kg*fiber2.getTangent(j);
        }
    }
    kc=exp(-(fiber1[0]-fiber2[0]).squaredNorm()/(lambdac*lambdac));
    kb=exp(-(fiber1[fiber1.size()-1]-fiber2[fiber2.size()-1]).squaredNorm()/(lambdab*lambdab));
    scalarProduct*=kc*kb;
    return scalarProduct;
}

float Metrique::scalarProductWeightedCurrentsParallel(Fiber &fiber1, Fiber &fiber2, float lambdac, float lambdab, float lambdag)
{
    float scalarProduct=0;
    float kc, kb, kg;
#pragma omp parallel for reduction (+:scalarProduct)
    for (unsigned int i=0; i<fiber1.size()-1; i++)
    {
        for (unsigned int j=0; j<fiber2.size()-1; j++)
        {
            kg=exp(-(fiber1.getCenter(i)-fiber2.getCenter(j)).squaredNorm()/(lambdag*lambdag));
            scalarProduct+=fiber1.getTangent(i).transpose()*kg*fiber2.getTangent(j);
        }
    }
    kc=exp(-(fiber1[0]-fiber2[0]).squaredNorm()/(lambdac*lambdac));
    kb=exp(-(fiber1[fiber1.size()-1]-fiber2[fiber2.size()-1]).squaredNorm()/(lambdab*lambdab));
    scalarProduct*=kc*kb;
    return scalarProduct;
}

//Weighted currents distance
float Metrique::WeightedCurrentsDistance(Fiber &fiber1, Fiber &fiber2, float lambdac, float lambdab, float lambdag)
{        
    fiber1.computeCentersAndTangents();
    fiber2.computeCentersAndTangents();
    float scalarProduct11=fiber1.getSelfScalarProduct(lambdag);//scalarProductWeightedCurrents(fiber1, fiber1, lambdac, lambdab, lambdag);
    float scalarProduct22=fiber2.getSelfScalarProduct(lambdag);//scalarProductWeightedCurrents(fiber2, fiber2, lambdac, lambdab, lambdag);
    float scalarProduct12=scalarProductWeightedCurrents(fiber1, fiber2, lambdac, lambdab, lambdag);
    return sqrt(scalarProduct11+scalarProduct22-2*scalarProduct12);
}

float Metrique::metriqueTest(Fiber fiber1, Fiber fiber2)
{
    fiber1.PutInSameDirection(fiber2);

    if (m_metric==0)
    {
        fiber1.computeCentersAndTangents();
        fiber2.computeCentersAndTangents();
        float scalarProduct = fabs(scalarProductWeightedCurrents(fiber1, fiber2, 6,6,6));
        return (-scalarProduct);
    }
    else if (m_metric==1)
    {

        //float fuzzyValue = 1-fabs(fiber1.getValueOfCylinder()-fiber2.getValueOfCylinder());
        return averageOfPointwiseEuclideanMetric(fiber1, fiber2);///fuzzyValue;
    }
    else if (m_metric==2)
    {
        return WeightedCurrentsDistance(fiber1, fiber2, 6, 6, 8);
    }
    else if (m_metric==3)
    {
        return MC(fiber1, fiber2);
    }
    else if (m_metric==4)
    {
        return VarifoldDistance(fiber1, fiber2, 15);
    }
    else if (m_metric==5)//WC + Fuzzy
    {
        fiber1.computeCentersAndTangents();
        fiber2.computeCentersAndTangents();
        float scalarProduct = fabs(scalarProductWeightedCurrents(fiber1, fiber2, 7,7,7));
//        float limitValue = 0.8;
//        float fuzzyValue;
//        if (fiber1.getValueOfCylinder()>limitValue || fiber2.getValueOfCylinder()>limitValue)
//            fuzzyValue=min(fiber1.getValueOfCylinder(), fiber2.getValueOfCylinder());
//        else
//            fuzzyValue=limitValue;

//        if (fiber1.getValueOfCylinder()>1 || fiber2.getValueOfCylinder()>1)
//            cerr << "Values not normalized" << endl;
        float fuzzyValue = 1-fabs(fiber1.getValueOfCylinder()-fiber2.getValueOfCylinder());

        //float fuzzyValue = (1-pow(fiber1.getValueOfCylinder()-fiber2.getValueOfCylinder(), 2))/(1.f-limitValue)+limitValue;
        return -scalarProduct*fuzzyValue;
    }
    else
    {
        Fiber fiberTemp=mergeOfFibers(fiber1, fiber2);
        return fiberTemp.computeVolume();
    }
}

//Distance between the endpoints of fibers oriented in the same direction
float Metrique::endpointsDistance(Fiber &fiber1, Fiber &fiber2)
{
    return (fiber1[0]-fiber2[0]).squaredNorm()+(fiber1[fiber1.size()-1]-fiber2[fiber2.size()-1]).squaredNorm();
}

Fiber Metrique::mergeOfFibers(Fiber& fiber1, Fiber& fiber2)
{
    fiber1.PutInSameDirection(fiber2);
    unsigned int nbFib1H = fiber1.getNbFibHierarchy();
    unsigned int nbFib2H = fiber2.getNbFibHierarchy();
    float factor1=(float)nbFib1H/(nbFib1H+nbFib2H);
    float factor2=(float)nbFib2H/(nbFib1H+nbFib2H);

    //Calcul de la noucelle fibre
    //Recherche des points les plus proches au niveau des extrémités
    unsigned int extremite1sens1=0;
    unsigned int extremite2sens1=0;
    unsigned int extremite1sens2=fiber1.size()-1;
    unsigned int extremite2sens2=fiber2.size()-1;
    float dist1sens1=(fiber1[0]-fiber2[0]).norm();
    float dist2sens1=dist1sens1;
    float dist1sens2=(fiber1[extremite1sens2]-fiber2[extremite2sens2]).norm();
    float dist2sens2=dist1sens2;
    for (unsigned int i=1; i<min(fiber1.size(), fiber2.size()); i++)
    {
        float temp=(fiber1[i]-fiber2[0]).norm();
        if (temp<dist1sens1)
        {
            dist1sens1=temp;
            extremite1sens1=i;
        }
        temp=(fiber1[0]-fiber2[i]).norm();
        if (temp<dist2sens1)
        {
            dist2sens1=temp;
            extremite2sens1=i;
        }
        temp=(fiber1[fiber1.size()-1-i]-fiber2[fiber2.size()-1]).norm();
        if (temp<dist1sens2)
        {
            dist1sens2=temp;
            extremite1sens2=fiber1.size()-1-i;
        }
        temp=(fiber1[fiber1.size()-1]-fiber2[fiber2.size()-1-i]).norm();
        if (temp<dist2sens2)
        {
            dist2sens2=temp;
            extremite2sens2=fiber2.size()-1-i;
        }
    }
//        unsigned int extremitesens1, extremitesens2;
    if (dist1sens1<dist2sens1)
        extremite2sens1=0;
    else
        extremite1sens1=0;
    if (dist1sens2<dist2sens2)
        extremite2sens2=fiber2.size()-1;
    else
        extremite1sens2=fiber1.size()-1;

    vector<Vector3f> nvPoints;
    vector<float> width;
    if (extremite1sens1!=0 && extremite2sens1!=0)
        cout << extremite1sens1 << " " << extremite2sens1 << endl;
    nvPoints.push_back((fiber1[0]*factor1+fiber2[0]*factor2));
    width.push_back(1);
    //Extrémité 1
    for (unsigned int i=1; i<extremite1sens1; i++)
    {
        nvPoints.push_back((fiber1[i]*factor1+fiber2[0]*factor2));
        float scalarProduct1, scalarProduct2, scalarProduct3;
        int pt1=1, pt2=1;
        if (extremite1sens1==extremite1sens2)
            pt1=-1;
        scalarProduct1=1-fabs((fiber1[extremite1sens1+pt1]-fiber1[extremite1sens1]).normalized().dot((fiber2[extremite2sens1]-fiber1[extremite1sens1]).normalized()));
//        scalarProduct1=1-fabs((fiber1[extremite1sens1+pt1]-fiber1[extremite1sens1]).normalized().dot((nvPoints[nvPoints.size()-1]-fiber1[extremite1sens1]).normalized()));
        if (extremite2sens1==extremite2sens2)
            pt2=-1;
        scalarProduct2=1-fabs((fiber2[extremite2sens1+pt2]-fiber2[extremite2sens1]).normalized().dot((fiber1[extremite1sens1]-fiber2[extremite2sens1]).normalized()));
//        scalarProduct2=1-fabs((fiber2[extremite2sens1+pt2]-fiber2[extremite2sens1]).normalized().dot((nvPoints[nvPoints.size()-1]-fiber2[extremite2sens1]).normalized()));
        //scalarProduct3=1;//fabs((fiber1.getTangent(extremite1sens1)*factor1+fiber2.getTangent(extremite2sens1)*factor2).normalized().dot((fiber1[extremite1sens1]-fiber2[extremite2sens1]).normalized()));
        scalarProduct3=1;//fabs((fiber1.getTangent(extremite1sens1)-fiber2.getTangent(extremite2sens1)).normalized().dot((nvPoints[nvPoints.size()-1]-(nvPoints[nvPoints.size()-2])).normalized()));
        width.push_back(max((fiber1.getWidth(extremite1sens1)*factor1*scalarProduct1+fiber2.getWidth(extremite2sens1)*factor2*scalarProduct2+(fiber1[extremite1sens1]-fiber2[extremite2sens1]).norm()*scalarProduct3)/2, max(fiber1.getWidth(extremite1sens1)*scalarProduct1, fiber2.getWidth(extremite2sens1)*scalarProduct2)));
        //width.push_back((fiber2.getWidth(extremite2sens1)+Vector3f(fiber2[extremite2sens1]-nvPoints[nvPoints.size()-1]).norm())*scalarProduct2);

        //width.push_back(fiber1.getWidth(i));
    }
    for (unsigned int i=1; i<extremite2sens1; i++)
    {
        nvPoints.push_back(fiber2[i]*factor2+fiber1[0]*factor1);
        float scalarProduct1, scalarProduct2, scalarProduct3;
        int pt1=1, pt2=1;
        if (extremite1sens1==extremite1sens2)
            pt1=-1;
        scalarProduct1=1-fabs((fiber1[extremite1sens1+pt1]-fiber1[extremite1sens1]).normalized().dot((fiber2[extremite2sens1]-fiber1[extremite1sens1]).normalized()));
//        scalarProduct1=1-fabs((fiber1[extremite1sens1+pt1]-fiber1[extremite1sens1]).normalized().dot((nvPoints[nvPoints.size()-1]-fiber1[extremite1sens1]).normalized()));
        if (extremite2sens1==extremite2sens2)
            pt2=-1;
        scalarProduct2=1-fabs((fiber2[extremite2sens1+pt2]-fiber2[extremite2sens1]).normalized().dot((fiber1[extremite1sens1]-fiber2[extremite2sens1]).normalized()));
        //scalarProduct2=1-fabs((fiber2[extremite2sens1+pt2]-fiber2[extremite2sens1]).normalized().dot((nvPoints[nvPoints.size()-1]-fiber2[extremite2sens1]).normalized()));
        //scalarProduct3=1;//fabs((fiber1.getTangent(extremite1sens1)*factor1+fiber2.getTangent(extremite2sens1)*factor2).normalized().dot((fiber1[extremite1sens1]-fiber2[extremite2sens1]).normalized()));
        scalarProduct3=1;//fabs((fiber1.getTangent(extremite1sens1)-fiber2.getTangent(extremite2sens1)).normalized().dot((nvPoints[nvPoints.size()-1]-(nvPoints[nvPoints.size()-2])).normalized()));
        width.push_back(max((fiber1.getWidth(extremite1sens1)*factor1*scalarProduct1+fiber2.getWidth(extremite2sens1)*factor2*scalarProduct2+(fiber1[extremite1sens1]-fiber2[extremite2sens1]).norm()*scalarProduct3)/2, max(fiber1.getWidth(extremite1sens1)*scalarProduct1, fiber2.getWidth(extremite2sens1)*scalarProduct2)));
        //width.push_back((fiber1.getWidth(extremite1sens1)+Vector3f(fiber1[extremite1sens1]-nvPoints[nvPoints.size()-1]).norm())*scalarProduct1);

        //width.push_back(fiber2.getWidth(i));
    }

    float pas=((float)max(extremite1sens2-extremite1sens1, extremite2sens2-extremite2sens1)/min(extremite1sens2-extremite1sens1, extremite2sens2-extremite2sens1));
    bool fibre1Unitaire=(min(extremite1sens2-extremite1sens1, extremite2sens2-extremite2sens1)==extremite1sens2-extremite1sens1);
    float cumul1=extremite1sens1;float cumul2=extremite2sens1;
    if (extremite1sens1==0 && extremite2sens1==0)
    {
        if (fibre1Unitaire)
        {
            cumul1++;
            cumul2+=pas;
        }
        else
        {
            cumul2++;
            cumul1+=pas;
        }
        extremite1sens1=cumul1;
        extremite2sens1=cumul2;
    }
    while(extremite1sens1<(extremite1sens2+1) && extremite2sens1<(extremite2sens2+1))
    {
        nvPoints.push_back(fiber1[extremite1sens1]*factor1+fiber2[extremite2sens1]*factor2);
        float scalarProduct1, scalarProduct2, scalarProduct3;
        int pt1=1, pt2=1;
        if (extremite1sens1==extremite1sens2)
            pt1=-1;
        scalarProduct1=1-fabs((fiber1[extremite1sens1+pt1]-fiber1[extremite1sens1]).normalized().dot((fiber2[extremite2sens1]-fiber1[extremite1sens1]).normalized()));
//        scalarProduct1=fabs((fiber1[extremite1sens1+pt1]-fiber1[extremite1sens1]).normalized().dot((nvPoints[nvPoints.size()-1]-(nvPoints[nvPoints.size()-2])).normalized()));
        if (extremite2sens1==extremite2sens2)
            pt2=-1;
        scalarProduct2=1-fabs((fiber2[extremite2sens1+pt2]-fiber2[extremite2sens1]).normalized().dot((fiber1[extremite1sens1]-fiber2[extremite2sens1]).normalized()));
        //scalarProduct2=fabs((fiber2[extremite2sens1+pt2]-fiber2[extremite2sens1]).normalized().dot((nvPoints[nvPoints.size()-1]-(nvPoints[nvPoints.size()-2])).normalized()));
        //scalarProduct3=1;//-fabs((fiber1.getTangent(extremite1sens1)*factor1+fiber2.getTangent(extremite2sens1)*factor2).normalized().dot((fiber1[extremite1sens1]-fiber2[extremite2sens1]).normalized()));
        scalarProduct3=1;//fabs((fiber1.getTangent(extremite1sens1)-fiber2.getTangent(extremite2sens1)).normalized().dot((nvPoints[nvPoints.size()-1]-(nvPoints[nvPoints.size()-2])).normalized()));
        width.push_back(max((fiber1.getWidth(extremite1sens1)*factor1*scalarProduct1+fiber2.getWidth(extremite2sens1)*factor2*scalarProduct2+(fiber1[extremite1sens1]-fiber2[extremite2sens1]).norm()*scalarProduct3)/2, max(fiber1.getWidth(extremite1sens1)*scalarProduct1, fiber2.getWidth(extremite2sens1)*scalarProduct2)));
        //width.push_back(max((fiber1.getWidth(extremite1sens1)+(fiber1[extremite1sens1]-nvPoints[nvPoints.size()-1]).norm())*scalarProduct1, (fiber2.getWidth(extremite2sens1)+(fiber2[extremite2sens1]-nvPoints[nvPoints.size()-1]).norm())*scalarProduct2));
        if (fibre1Unitaire)
        {
            cumul1++;
            cumul2+=pas;
        }
        else
        {
            cumul2++;
            cumul1+=pas;
        }
        extremite1sens1=cumul1;
        extremite2sens1=cumul2;
    }

    //Extrémité 2
    for (unsigned int i=extremite1sens2+1; i<fiber1.size(); i++)
    {
        nvPoints.push_back(fiber1[i]*factor1+fiber2[fiber2.size()-1]*factor2);
        width.push_back(width[width.size()-1]);// fiber1.getWidth(i));
    }
    for (unsigned int i=extremite2sens2+1; i<fiber2.size(); i++)
    {
        nvPoints.push_back(fiber2[i]*factor2+fiber1[fiber1.size()-1]*factor1);
        width.push_back(width[width.size()-1]);//(fiber2.getWidth(i));
    }

    width[0]=width[1];

    MatrixXf points;
    points.resize(nvPoints.size(), 3);
    for (unsigned int i=0; i<nvPoints.size(); i++)
        points.row(i)=nvPoints[i];

    Fiber newFiber(points, nvPoints.size());
    newFiber.setWidth(width);
    newFiber.addFibHierarchy(nbFib1H+nbFib2H);
    return newFiber;
}

Transfo computeNewTransfo(Fiber& fiber1, unsigned int positionFib1, Fiber& fiber2, unsigned int positionFib2)
{
    //a : grand axe de l'ellipse
    //b : petit axe de l'ellipse
    Transfo newTransfo;
    Vector3f direction=fiber1[positionFib1]-fiber2[positionFib2];
    float distBetweenEllipses=direction.norm();
    direction.normalize();
    newTransfo.center=(fiber1[positionFib1]+fiber2[positionFib2])/2.;
    newTransfo.rotation=((fiber1.getProfileTransform(positionFib1).rotation+fiber2.getProfileTransform(positionFib2).rotation)/2.).normalized();
    /// Non englobage de toutes les fibres
//    float angle1=acos(fiber1.getProfileTransform(positionFib1).ellipseOrientation.dot(direction.normalized()));
//    float value1=fiber1.getProfileTransform(positionFib1).b/sqrt(1-(1-pow(fiber1.getProfileTransform(positionFib1).b/fiber1.getProfileTransform(positionFib1).a,2))*pow(cos(angle1),2));
//    float angle2=acos(fiber2.getProfileTransform(positionFib2).ellipseOrientation.dot(direction.normalized()));
//    float value2=fiber2.getProfileTransform(positionFib2).b/sqrt(1-(1-pow(fiber2.getProfileTransform(positionFib2).b/fiber2.getProfileTransform(positionFib2).a,2))*pow(cos(angle2),2));

//    float value1=fiber1.getProfileTransform(positionFib1).a*fabs(fiber1.getProfileTransform(positionFib1).ellipseOrientation.dot(direction.normalized()));
//    float value2=fiber2.getProfileTransform(positionFib2).a*fabs(fiber2.getProfileTransform(positionFib2).ellipseOrientation.dot(direction.normalized()));

    //Repère de l'ellipse
    float p2=pow(tan(acos(min(1.f,fabs((fiber1.getProfileTransform(positionFib1).ellipseOrientation.cross(fiber1.getProfileTransform(positionFib1).rotation).normalized()).dot(direction))))),2);
    //cout << p2 << " " << sqrt(p2) << " ";
    float value1=2*sqrt(pow(fiber1.getProfileTransform(positionFib1).a,2)*p2+pow(fiber1.getProfileTransform(positionFib1).b,2))/sqrt(p2+1);
    p2=pow(tan(acos(min(1.f,fabs((fiber2.getProfileTransform(positionFib2).ellipseOrientation.cross(fiber2.getProfileTransform(positionFib2).rotation).normalized()).dot(direction))))),2);
    //cout << p2 << " " << sqrt(p2) << endl;
    float value2=2*sqrt(pow(fiber2.getProfileTransform(positionFib2).a,2)*p2+pow(fiber2.getProfileTransform(positionFib2).b,2))/sqrt(p2+1);

    newTransfo.center+=direction*(value1-value2)/4.;
    //cout << "value1: " << value1 << " value2: " << value2 << " ";
    float axe1, axe2;
    //if (value1<direction.norm() && value2<direction.norm())
    axe1=(distBetweenEllipses+value1/2.+value2/2.)/2.;
    //else if (value1>value2)
    //    axe1=value1;
    //else
    //    axe1=value2;
    //    value1=fiber1.getProfileTransform(positionFib1).b/sqrt(1-(1-pow(fiber1.getProfileTransform(positionFib1).b/fiber1.getProfileTransform(positionFib1).a,2))*pow(cos(angle1+M_PI_2),2));
    //    value2=fiber2.getProfileTransform(positionFib2).b/sqrt(1-(1-pow(fiber2.getProfileTransform(positionFib2).b/fiber2.getProfileTransform(positionFib2).a,2))*pow(cos(angle2+M_PI_2),2));
    //    value1=max((float)sqrt(max(-pow(value1,2)+pow(fiber1.getProfileTransform(positionFib1).a,2),0.)),
    //               fiber1.getProfileTransform(positionFib1).b*fabs(fiber1.getProfileTransform(positionFib1).ellipseOrientation.dot(direction.normalized().cross(newTransfo.rotation))));
    //    value2=max((float)sqrt(max(-pow(value2,2)+pow(fiber2.getProfileTransform(positionFib2).a,2),0.)),
    //               fiber2.getProfileTransform(positionFib2).b*fabs(fiber2.getProfileTransform(positionFib2).ellipseOrientation.dot(direction.normalized().cross(newTransfo.rotation))));

    p2=pow(tan(acos(min(1.f,fabs(fiber1.getProfileTransform(positionFib1).ellipseOrientation.dot(direction)))))/**fiber1.getProfileTransform(positionFib1).a*/,2);
    value1=2*sqrt(pow(fiber1.getProfileTransform(positionFib1).a,2)*p2+pow(fiber1.getProfileTransform(positionFib1).b,2))/sqrt(p2+1);
    p2=pow(tan(acos(min(1.f,fabs(fiber2.getProfileTransform(positionFib2).ellipseOrientation.dot(direction)))))/**fiber2.getProfileTransform(positionFib2).a*/,2);
    value2=2*sqrt(pow(fiber2.getProfileTransform(positionFib2).a,2)*p2+pow(fiber2.getProfileTransform(positionFib2).b,2))/sqrt(p2+1);
    axe2=max(value1, value2)/2.;


    newTransfo.a=max(axe1,axe2);
    newTransfo.b=min(axe1,axe2);
    //cout << newTransfo.a << " " << axe1 << " " << newTransfo.b << " " << axe2 << endl;
    //if (isnan(axe2))
    //cout << fiber1.getProfileTransform(positionFib1).a*fabs(fiber1.getProfileTransform(positionFib1).ellipseOrientation.dot(direction.normalized())) << " " << fiber1.getProfileTransform(positionFib1).a << endl;
    //cout << "value1 :" << value1 << " value2: " << value2 << " a: " << newTransfo.a << " b: " << newTransfo.b << endl;
    if (newTransfo.a==axe2)
        newTransfo.ellipseOrientation=(direction.normalized().cross(newTransfo.rotation)).normalized();
    else
        newTransfo.ellipseOrientation=direction.normalized();
    return newTransfo;
}

Fiber Metrique::mergeOfFibersV2(Fiber& fiber1, Fiber& fiber2)
{
    //cout << endl;
    fiber1.PutInSameDirection(fiber2);
    unsigned int nbFib1H = fiber1.getNbFibHierarchy();
    unsigned int nbFib2H = fiber2.getNbFibHierarchy();
    float factor1=(float)nbFib1H/(nbFib1H+nbFib2H);
    float factor2=(float)nbFib2H/(nbFib1H+nbFib2H);

    //Calcul de la noucelle fibre
    //Recherche des points les plus proches au niveau des extrémités
    unsigned int extremite1sens1=0;
    unsigned int extremite2sens1=0;
    unsigned int extremite1sens2=fiber1.size()-1;
    unsigned int extremite2sens2=fiber2.size()-1;
    float dist1sens1=(fiber1[0]-fiber2[0]).norm();
    float dist2sens1=dist1sens1;
    float dist1sens2=(fiber1[extremite1sens2]-fiber2[extremite2sens2]).norm();
    float dist2sens2=dist1sens2;
    for (unsigned int i=1; i<min(fiber1.size(), fiber2.size()); i++)
    {
        float temp=(fiber1[i]-fiber2[0]).norm();
        if (temp<dist1sens1)
        {
            dist1sens1=temp;
            extremite1sens1=i;
        }
        temp=(fiber1[0]-fiber2[i]).norm();
        if (temp<dist2sens1)
        {
            dist2sens1=temp;
            extremite2sens1=i;
        }
        temp=(fiber1[fiber1.size()-1-i]-fiber2[fiber2.size()-1]).norm();
        if (temp<dist1sens2)
        {
            dist1sens2=temp;
            extremite1sens2=fiber1.size()-1-i;
        }
        temp=(fiber1[fiber1.size()-1]-fiber2[fiber2.size()-1-i]).norm();
        if (temp<dist2sens2)
        {
            dist2sens2=temp;
            extremite2sens2=fiber2.size()-1-i;
        }
    }
    if (dist1sens1<dist2sens1)
        extremite2sens1=0;
    else
        extremite1sens1=0;
    if (dist1sens2<dist2sens2)
        extremite2sens2=fiber2.size()-1;
    else
        extremite1sens2=fiber1.size()-1;

    vector<Vector3f> nvPoints;
    vector<float> width;
    vector<Transfo> profileTransform;
    if (extremite1sens1!=0 && extremite2sens1!=0)
        cout << extremite1sens1 << " " << extremite2sens1 << endl;
    //nvPoints.push_back((fiber1[0]*factor1+fiber2[0]*factor2));
    //width.push_back(1);
    //profileTransform.push_back(computeNewTransfo(fiber1, 0, fiber2, 0));
    //Extrémité 1
    for (unsigned int i=0; i<extremite1sens1; i++)
    {
        nvPoints.push_back((fiber1[i]*factor1+fiber2[0]*factor2));
        profileTransform.push_back(fiber1.getProfileTransform(i));
        //profileTransform.push_back(computeNewTransfo(fiber1, i, fiber2, 0));
    }
    for (unsigned int i=0; i<extremite2sens1; i++)
    {
        nvPoints.push_back(fiber2[i]*factor2+fiber1[0]*factor1);
        profileTransform.push_back(fiber2.getProfileTransform(i));
        //profileTransform.push_back(computeNewTransfo(fiber1, 0, fiber2, i));
    }

    float pas=((float)max(extremite1sens2-extremite1sens1, extremite2sens2-extremite2sens1)/min(extremite1sens2-extremite1sens1, extremite2sens2-extremite2sens1));
    bool fibre1Unitaire=(min(extremite1sens2-extremite1sens1, extremite2sens2-extremite2sens1)==extremite1sens2-extremite1sens1);
    float cumul1=extremite1sens1;float cumul2=extremite2sens1;
//    if (extremite1sens1==0 && extremite2sens1==0)
//    {
//        if (fibre1Unitaire)
//        {
//            cumul1++;
//            cumul2+=pas;
//        }
//        else
//        {
//            cumul2++;
//            cumul1+=pas;
//        }
//        extremite1sens1=cumul1;
//        extremite2sens1=cumul2;
//    }
    //cout << "pas : " << pas << endl;
    while(extremite1sens1<(extremite1sens2) && extremite2sens1<(extremite2sens2))
    {
        nvPoints.push_back(fiber1[extremite1sens1]*factor1+fiber2[extremite2sens1]*factor2);
        profileTransform.push_back(computeNewTransfo(fiber1, extremite1sens1, fiber2, extremite2sens1));
        if (fibre1Unitaire)
        {
            cumul1++;
            extremite1sens1=cumul1;
            cumul2+=pas;
            float dist=(fiber1[extremite1sens1]-fiber2[cumul2]).norm();
            extremite2sens1=cumul2;
            for (int j=-20; j<30; j++)
            {
                if (cumul2+j>=0 && cumul2+j<fiber2.size())
                {
                    float distTemp=(fiber1[extremite1sens1]-fiber2[cumul2+j]).norm();
                    if (distTemp<dist)
                    {
                        extremite2sens1=cumul2+j;
                        dist=distTemp;
                    }
                }
            }
        }
        else
        {
            cumul2++;
            extremite2sens1=cumul2;
            cumul1+=pas;
            float dist=(fiber1[cumul1]-fiber2[extremite2sens1]).norm();
            extremite1sens1=cumul1;
            for (int j=-20; j<30; j++)
            {
                if (cumul1+j>=0 && cumul1+j<fiber1.size())
                {
                    float distTemp=(fiber1[cumul1+j]-fiber2[extremite2sens1]).norm();
                    if (distTemp<dist)
                    {
                        extremite1sens1=cumul1+j;
                        dist=distTemp;
                    }
                }
            }
        }
//        extremite1sens1=cumul1;
//        extremite2sens1=cumul2;
//        float theoreticalDistance=(fiber1[extremite1sens1]-fiber2[extremite2sens1]).norm();
//        if (fibre1Unitaire)
//        {
//            if (theoreticalDistance>((fiber1[extremite1sens1]-fiber2[extremite2sens1-1]).norm()))
//            {
//                extremite2sens1-=1;
//                cumul2-=1;
//            }
//        }
//        else
//        {
//            if (theoreticalDistance>((fiber1[extremite1sens1-1]-fiber2[extremite2sens1]).norm()))
//            {
//                extremite1sens1-=1;
//                cumul1-=1;
//            }
//        }
    }

    //Extrémité 2
    for (unsigned int i=extremite1sens2+1; i<fiber1.size(); i++)
    {
        nvPoints.push_back(fiber1[i]*factor1+fiber2[fiber2.size()-1]*factor2);
        profileTransform.push_back(fiber1.getProfileTransform(i));
        //profileTransform.push_back(computeNewTransfo(fiber1, i, fiber2, fiber2.size()-1));
    }
    for (unsigned int i=extremite2sens2+1; i<fiber2.size(); i++)
    {
        nvPoints.push_back(fiber2[i]*factor2+fiber1[fiber1.size()-1]*factor1);
        profileTransform.push_back(fiber2.getProfileTransform(i));
        //profileTransform.push_back(computeNewTransfo(fiber1, fiber1.size()-1, fiber2, i));
    }


    MatrixXf points;
    points.resize(nvPoints.size(), 3);
    for (unsigned int i=0; i<nvPoints.size(); i++)
        points.row(i)=nvPoints[i];

    Fiber newFiber(points, nvPoints.size());
    //Mode ellipse for every point
#if (PROFILE_PER_POINT)
    newFiber.setProfileTransform(profileTransform);

    //Mode ellipse with transformations along the fiber
#else
    Profile profile;
    profile.ellipse.a=0;
    profile.ellipse.b=0;
    for (unsigned int i=0; i<profileTransform.size(); i++)
    {
        profile.ellipse.a+=profileTransform[i].a;
        profile.ellipse.b+=profileTransform[i].b;
    }
    profile.ellipse.a/=(float)profileTransform.size();
    profile.ellipse.b/=(float)profileTransform.size();
    newFiber.setProfile(profile);
    Vector3f lastOrientation=profileTransform[0].ellipseOrientation;
    for (unsigned int i=0; i<profileTransform.size(); i++)
    {
        Matrix3f matrice;
        if (profileTransform[i].ellipseOrientation.dot(lastOrientation)<0)
            lastOrientation=-profileTransform[i].ellipseOrientation;
        else
            lastOrientation=profileTransform[i].ellipseOrientation;
        Vector3f axeY=profileTransform[i].rotation.cross(lastOrientation);
        matrice << lastOrientation(0), axeY(0), profileTransform[i].rotation(0),
                   lastOrientation(1), axeY(1), profileTransform[i].rotation(1),
                   lastOrientation(2), axeY(2), profileTransform[i].rotation(2);
        newFiber.setTransform(i, matrice);
    }
#endif

    newFiber.addFibHierarchy(nbFib1H+nbFib2H);
    return newFiber;
}
