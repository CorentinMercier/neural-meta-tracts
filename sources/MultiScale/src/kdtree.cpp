#include "../include/kdtree.h"

/// Par convention, les éléments sont sur les colonnes, la première ligne de chaque élément est un entier identifiant l'élément
/// Les lignes suivantes sont les différentes variables qui vont servir à la construction du kdtree
kdtree::kdtree(vector<vector<float> > elements, unsigned int depth)
{
    unsigned int dimension=elements[0].size();
    int numberOfElements=elements.size();
    unsigned int actualDimension=depth%(dimension-1);
    sort(elements.begin(), elements.end(), compareElements(actualDimension+1));
    m_location=elements[numberOfElements/2];
    m_axisMin=elements[0][actualDimension+1];
    m_axisMax=elements[elements.size()-1][actualDimension+1];
    m_depth=depth;
    if (numberOfElements>2)
    {
        vector<vector<float> > subvector (elements.begin()+numberOfElements/2+1, elements.end());
        m_rightChild=new kdtree(subvector, depth+1);
    }
    else
        m_rightChild=NULL;
    if(numberOfElements/2>0)
    {
        vector<vector<float> > subvector (elements.begin(), elements.begin()+numberOfElements/2);
        m_leftChild=new kdtree(subvector, depth+1);
    }
    else
        m_leftChild=NULL;

}

void kdtree::getNeighbours(vector<float> & objectif, vector<unsigned int> &listNeighbours, unsigned int numberOfNeighbours)
{
    listNeighbours.resize(numberOfNeighbours);
    vector<vector<float> > neighbours(numberOfNeighbours);
    //Aucun voisin n'est encore trouvé, donc la distance à l'objectif est mawimale
    for (unsigned int i=0; i<numberOfNeighbours; i++)
    {
        neighbours[i]=m_location;
        neighbours[i].push_back(std::numeric_limits<float>::max());
    }
    this->getNeighboursRecursive(objectif, neighbours);
    //Enregistrement des numéros des fibres ayant été déterminées comme voisines
    for (unsigned int i=0; i<numberOfNeighbours; i++)
        listNeighbours[i]=neighbours[i][0];
}

void kdtree::getNeighboursRecursive(vector<float> & objectif, vector<vector<float>> & neighbours)
{
    //Détermintaion du voisin le plus loin de l'objectif
    std::vector<vector<float> >::iterator result;
    result = std::max_element(neighbours.begin(), neighbours.end(), compareElements(neighbours[0].size()-1));
    float currentDistance=distance(objectif);
    float currentMax=neighbours[std::distance(neighbours.begin(), result)][neighbours[0].size()-1];
    //Si la distance du noeud actuel à l'objectif est plus petite que la distance maximale, l'ancien noeud est remplacé par le nouveau
    if (currentDistance<currentMax)//Le noeud est plus près de l'objectif qu'au moins un des éléments déjà parcourus
    {
        if (currentDistance!=0)//Le noeud n'est pas voisin de lui-même
        {
            neighbours[std::distance(neighbours.begin(), result)]=m_location;
            neighbours[std::distance(neighbours.begin(), result)].push_back(currentDistance);
            result = std::max_element(neighbours.begin(), neighbours.end(), compareElements(neighbours[0].size()-1));
            currentMax=neighbours[std::distance(neighbours.begin(), result)][neighbours[0].size()-1];;
        }
    }
    //Les sous-arbres ne sont parcourus que lorsque la distance de l'objectif est plus grande que la distance à la boite englobante du sous-arbre
    unsigned int nextDimension=((m_depth+1)%(objectif.size()-1))+1;
    if (m_leftChild && min(pow(objectif[nextDimension]-m_leftChild->getMinAxis(),2),pow(objectif[nextDimension]-m_leftChild->getMaxAxis(),2))<=currentMax)
        m_leftChild->getNeighboursRecursive(objectif, neighbours);
    if (m_rightChild && min(pow(objectif[nextDimension]-m_rightChild->getMinAxis(),2),pow(objectif[nextDimension]-m_rightChild->getMaxAxis(),2))<=currentMax)
        m_rightChild->getNeighboursRecursive(objectif, neighbours);
    //Si le noeud est plus loin que tous ceux déjà parcourus, on arrête le parcours
}

///
/// \brief kdtree::getNeighboursByDistance
/// Fonction renvoyant les éléments dont la distance selon les points 3D à l'objectif est inférieure à la cible
/// \param objectif
/// \param listNeighbours
/// \param distanceCible
///

void kdtree::getNeighboursByDistance(vector<float> &objectif, vector<unsigned int> &listNeighbours, float distanceCible)
{
    vector<vector<float> > neighbours;
    float currentDistance=maxDistance3DPoints(objectif);
    if (currentDistance<distanceCible)
        neighbours.push_back(m_location);
    unsigned int nextDimension=((m_depth+1)%(objectif.size()-1))+1;
    if (m_leftChild && min(fabs(objectif[nextDimension]-m_leftChild->getMinAxis()),fabs(objectif[nextDimension]-m_leftChild->getMaxAxis()))<=distanceCible)
        m_leftChild->getNeighboursByDistanceRecursive(objectif, neighbours, distanceCible);
    if (m_rightChild && min(fabs(objectif[nextDimension]-m_rightChild->getMinAxis()),fabs(objectif[nextDimension]-m_rightChild->getMaxAxis()))<=distanceCible)
        m_rightChild->getNeighboursByDistanceRecursive(objectif, neighbours, distanceCible);
    for (unsigned int i=0; i<neighbours.size(); i++)
        listNeighbours.push_back(neighbours[i][0]);
}

void kdtree::getNeighboursByDistanceRecursive(vector<float> &objectif, vector<vector<float> > &neighbours, float distanceCible)
{
    float currentDistance=maxDistance3DPoints(objectif);
    if (currentDistance<distanceCible)
        neighbours.push_back(m_location);
    unsigned int nextDimension=((m_depth+1)%(objectif.size()-1))+1;
    if (m_leftChild && min(fabs(objectif[nextDimension]-m_leftChild->getMinAxis()),fabs(objectif[nextDimension]-m_leftChild->getMaxAxis()))<=distanceCible)
        m_leftChild->getNeighboursByDistanceRecursive(objectif, neighbours, distanceCible);
    if (m_rightChild && min(fabs(objectif[nextDimension]-m_rightChild->getMinAxis()),fabs(objectif[nextDimension]-m_rightChild->getMaxAxis()))<=distanceCible)
        m_rightChild->getNeighboursByDistanceRecursive(objectif, neighbours, distanceCible);
}


float kdtree::distance(vector<float> & objectif)
{
    float dist=0;
    for (unsigned int i=1; i<objectif.size(); i++)
        dist+=pow(objectif[i]-m_location[i], 2);
    return dist;
}

float kdtree::maxDistance3DPoints(vector<float> &objectif)
{
    float dist=0;
    for (unsigned int i=1; i<objectif.size(); i+=3)
    {
        float tempDist=0;
        for (int j=0; j<3; j++)
             tempDist+=pow(objectif[j+i]-m_location[j+i], 2);
        dist=max(dist, tempDist);
    }
    return sqrt(dist);
}

void kdtree::addElement(vector<float> &element)
{
    this->addElementRecursive(element, 0);
}

void kdtree::addElementRecursive(vector<float> &element, unsigned int depth)
{
    unsigned int dimension=element.size();
    unsigned int actualDimension=depth%(dimension-1);
    if (element[actualDimension+1]<m_location[actualDimension+1])
    {
        if (m_leftChild)
            m_leftChild->addElementRecursive(element, depth+1);
        else
        {
            vector<vector<float> > elements;
            elements.push_back(element);
            m_leftChild=new kdtree(elements , depth+1);
        }
    }
    else
    {
        if (m_rightChild)
            m_rightChild->addElementRecursive(element, depth+1);
        else
        {
            vector<vector<float> > elements;
            elements.push_back(element);
            m_rightChild=new kdtree(elements, depth+1);
        }
    }
}

unsigned int kdtree::count()
{
    unsigned int total=0;
    total+=this->countRec(total);
    return total;
}

unsigned int kdtree::countRec(unsigned int actualCount)
{
    actualCount++;
    if (m_rightChild) actualCount=m_rightChild->countRec(actualCount);
    if (m_leftChild) actualCount=m_leftChild->countRec(actualCount);
    return actualCount;
}

void kdtree::getElementAtDepth(vector<int> & element, unsigned int depth)
{
    getElementAtDepthRecursive(element, depth, 0);
}

void kdtree::getElementAtDepthRecursive(vector<int> & element, unsigned int depth, unsigned int currentDepth)
{
    if (currentDepth==depth)
        element.push_back(m_location[0]);
    else
    {
        if (m_rightChild) m_rightChild->getElementAtDepthRecursive(element, depth, currentDepth+1);
        if (m_leftChild) m_leftChild->getElementAtDepthRecursive(element, depth, currentDepth+1);
    }
}
