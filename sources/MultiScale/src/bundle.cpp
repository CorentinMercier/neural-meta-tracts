#include "../include/bundle.h"
#include "../../interface/widgetopengl.h"

typedef struct {
    double r,g,b;
} COLOUR;

COLOUR GetColour(double v,double vmin,double vmax)
{
   COLOUR c = {1.0,1.0,1.0}; // white
   double dv;

   if (v < vmin)
      v = vmin;
   if (v > vmax)
      v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) {
      c.r = 0;
      c.g = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5 * dv)) {
      c.r = 0;
      c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   } else if (v < (vmin + 0.75 * dv)) {
      c.r = 4 * (v - vmin - 0.5 * dv) / dv;
      c.b = 0;
   } else {
      c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
      c.b = 0;
   }

   return(c);
}

void advanceBarCout(float percent)
{
    int nbBar = int(percent/5);
    if(percent < 100)
        std::cout << std::setfill('0') << std::setw(2) << int(percent) << " % |";
    else
    {
        std::cout << "Done |";
        return;
    }
    for(int i = 0; i < nbBar; ++i)
        std::cout << "–";
    for(int i = nbBar; i < 20; ++i)
        std::cout << " ";
    std::cout << "|\r";
    std::cout.flush();
}

Bundle::Bundle(string fileFuzzySet, string fileEndpoints1, string fileEndpoints2)
{
    if (!fileFuzzySet.empty())
    {
        cout << " Load fuzzy sets from file " << fileFuzzySet << " ... " << endl;
        m_fuzzySet = new Fuzzy(fileFuzzySet);
        cout << "   [OK] " << endl;
        cout << " Load mask 1 from file " << fileEndpoints1 << " ... " << endl;
        m_maskEndpoints1 = new Fuzzy(fileEndpoints1);
        m_maskEndpoints1->computeDistanceField();
        cout << "   [OK] " << endl;
        cout << " Load mask 2 from file " << fileEndpoints2 << " ... " << endl;
        m_maskEndpoints2 = new Fuzzy(fileEndpoints2);
        m_maskEndpoints2->computeDistanceField();
        cout << "   [OK] " << endl;
    }
}
Bundle::Bundle(){}

Bundle::~Bundle()
{
    m_fibers.clear();
    m_changes.clear();
    delete(m_fuzzySet);
    delete(m_maskEndpoints1);
    delete(m_maskEndpoints2);
    m_paires = priority_queue<paire, vector<paire>, comparisonPaires >();
    m_pairesVector.clear();
}

bool hasEnding(const string & path, const string & extension)
{
    if (path.length()<extension.length()) return false;
    return (path.compare(path.length()-extension.length(), extension.length(), extension) == 0);
}

void Bundle::loadFibers(string filename, bool reduceNbFibers)
{
    srand(time(NULL));
     if (hasEnding(filename, ".tck"))
        loadTckFibers(filename, reduceNbFibers);
     else if (hasEnding(filename, ".vtk"))
         loadVtkFibers(filename, reduceNbFibers);
     else
     {
         cerr << "Error, file format is not supported" << endl;
         exit(EXIT_FAILURE);
     }
}

void Bundle::loadVtkFibers(const string & filename, bool reduceNbFibers)
{
#ifdef USE_VTK
    // Get all data from the file
    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // All of the standard data types can be checked and obtained like this:
    if(!reader->IsFilePolyData())
    {
        cout << "output is not a polydata" << endl;
    }
    cout << "output is a polydata" << endl;
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData = reader->GetPolyDataOutput();

    unsigned int numberOfPoints = polyData->GetNumberOfPoints();
    unsigned int numberOfFibers = polyData->GetNumberOfLines();
    cout << "This bundle has " << numberOfPoints << " points and " << numberOfFibers << " fibers." << endl;

    vtkSmartPointer<vtkCellArray> Lines = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType *indices;
    vtkIdType length;
    VectorXi NumberPointsPerFiber;
    MatrixXf Points;

    Lines = polyData->GetLines();

    NumberPointsPerFiber.resize(numberOfFibers);
    int lineCount = 0;
    for (Lines->InitTraversal(); Lines->GetNextCell(length, indices); lineCount++)
    {
        NumberPointsPerFiber(lineCount)=length;
    }
    if( (unsigned int)NumberPointsPerFiber.sum() != numberOfPoints )
    {
        cout << NumberPointsPerFiber.sum() << " != " << numberOfPoints << endl;
        throw runtime_error("Total Number of Points Mismatched!");
    }

    Points.resize(numberOfPoints, 3);
    for (unsigned int i = 0; i < numberOfPoints; i++)
    {
        double p[3];
        polyData->GetPoint(i, p);
        for (int dim = 0; dim < 3; dim++)
        {
            float pf = (float)(p[dim]);
            Points(i, dim) = pf;
        }
    }

//    Bundle removedFibers;

    // List Points per Fiber
    if (reduceNbFibers)
    {
        unsigned int divisionValue=5;
        numberOfFibers/=divisionValue;
        unsigned int temp = 0;

        for (unsigned int lineCount = 0; lineCount<numberOfFibers*divisionValue; lineCount++)
        {
            if (lineCount%divisionValue==0)
            {
                Fiber newFib(Points.block(temp,0,NumberPointsPerFiber(lineCount),3), NumberPointsPerFiber(lineCount,0));
//                if (newFib.length()<120)
//                {
                    newFib.resample(100);
                    this->addFiber(newFib);
//                }
//                else
//                {
//                    removedFibers.addFiber(newFib);
//                }
            }
            temp+=NumberPointsPerFiber(lineCount);
        }
    }
    else
    {
        unsigned int temp = 0;

        for (unsigned int lineCount = 0; lineCount<numberOfFibers; lineCount++)
        {
            Fiber newFib(Points.block(temp,0,NumberPointsPerFiber(lineCount),3), NumberPointsPerFiber(lineCount,0));
//            if (!(newFib.length()<40))
//            {
                newFib.resample(100);
                this->addFiber(newFib);
//            }
//            else
//            {
//                removedFibers.addFiber(newFib);
//            }
            temp+=NumberPointsPerFiber(lineCount);
        }
    }
    //removedFibers.saveCurrentFibersAsVTK("Removed.vtk");
    //cout << "Removed fibers saved" << endl;

    //    string filepath = filename;
//    while(filepath.back()!='/')
//        filepath.pop_back();
//    filepath = filepath+"FibersWithLengthGreaterThan40mm.vtk";
//    this->saveCurrentFibersAsVTK(filepath);
//    cout << "File " << filepath << " saved" << endl;

    for (unsigned int i=0; i<this->size(); i++)
        m_fibers[i].addToHierarchy(i);
#endif
}

bool isnan(const Vector3f v)
{
    return ((isnan(v.x()) || isnan(v.y())) || isnan(v.z()));
}

bool isinf(const Vector3f v)
{
    return ((isinf(v.x()) || isinf(v.y())) || isinf(v.z()));
}

///
/// \brief loadTckFibers Function to load fibers from a tck file
/// \param filename Path of the file to load
/// \param fibers Table in which the loaded bundle is stored
///

void Bundle::loadTckFibers(const string & filename, bool reduceNbFibers)
{
    unsigned int nbOfFibers;
    if (filename.find(".tck")==string::npos)
    {
        cerr << "File must be in the tck format" << endl;
        exit(EXIT_FAILURE);
    }
    fstream file;
    file.open(filename, ios_base::in | ios_base::binary);
    if (!file.is_open())
    {
        cerr << "Impossible to find file " << filename << endl;
        exit(EXIT_FAILURE);
    }
    string message;
    file >> message;
    if (message!="mrtrix")
    {
        cerr << "Tck file not correct" << endl;
        exit(EXIT_FAILURE);
    }
    file >> message;
    if (message!="tracks")
    {
        cerr << "Tck file not correct" << endl;
        exit(EXIT_FAILURE);
    }
    unsigned offset=0;
    bool count = false, float32le = false;
    while (message!="END")
    {
        file >> message;
        if (message=="count:")
        {
            file >> message;
            cout << stoi(message) << " fibers found ..." << endl;
            nbOfFibers = stoi(message);
            count = true;
        }
        if (message.find("datatype")!=string::npos)
        {
            file >> message;
            if (message!="Float32LE")
            {
                cerr << "Format not supported: " << message << endl;
                exit(EXIT_FAILURE);
            }
            else
                float32le = true;
        }
        if (message.find("file")!=string::npos)
        {
            file >> message;//"."
            file >> message;
            offset = stoi(message);
        }
    }
    if (offset==0)
    {
        cerr << "file needs to be specified in tck file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
    if (!(float32le && count))
    {
        cerr << "Error with file header, count or datatype not specified" << endl;
        exit(EXIT_FAILURE);
    }
    file.seekg(offset, file.beg);
    Vector3f newPoint;
    newPoint << 0,0,0;
    std::vector<Vector3f> newFiber;
    newFiber.clear();
    unsigned int divisionValue=5;
    nbOfFibers/=divisionValue;
    unsigned int lineCount = 0;
    while (!isinf(newPoint))
    {
        file.read(reinterpret_cast<char *>(&newPoint), sizeof (Vector3f));
        if (isnan(newPoint) || isinf(newPoint))
        {
            if (!reduceNbFibers || lineCount%divisionValue==0)
            {
                MatrixXf Points;
                Points.resize(newFiber.size(), 3);
                for (unsigned int i=0; i<newFiber.size(); i++)
                    Points.row(i) = newFiber[i];
                Fiber newFib(Points, newFiber.size());
                newFib.resample(100);
                this->addFiber(newFib);
            }
            lineCount++;
            newFiber.clear();
            if (isinf(newPoint))
            {
                break;
            }
            file.read(reinterpret_cast<char *>(&newPoint), sizeof (Vector3f));
            if (!isnan(newPoint) && !isinf(newPoint))
                newFiber.push_back(newPoint);
        }
        else
        {
            newFiber.push_back(newPoint);
        }
    }
    file.close();
}

vector<vector<Vector3f> > Bundle::getFibers()
{
    vector<vector<Vector3f> > fibres(m_fibers.size());
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {fibres[i].resize(m_fibers[i].size());
        for (unsigned int j=0; j<m_fibers[i].size(); j++)
            fibres[i][j]=m_fibers[i][j];
    }
    return fibres;
}

Fiber& Bundle::operator[](unsigned int element)
{
    if (element>=m_fibers.size())
    {
        cerr << "Erreur : accès à un élément inexistant" << endl;
        exit(EXIT_FAILURE);
    }
    return m_fibers[element];
}

void Bundle::addFiber(Fiber &new_fiber)
{
#if (TEST_FIBERS_VALIDITY)
    {
        float ecart=0.00001f;
        bool same=false;
        unsigned int i=0;
        while (!same && i<m_fibers.size())
        {
            new_fiber.PutInSameDirection(m_fibers[i]);
            if (new_fiber[0]==m_fibers[i][0] && new_fiber.size()==m_fibers[i].size())
            {
                same=true;
                for (unsigned int j=0; j<new_fiber.size(); j++)
                    same = same && (new_fiber[j]==m_fibers[i][j]);
            }
            if (!same && (new_fiber[0]==m_fibers[i][0] || new_fiber[new_fiber.size()-1]==m_fibers[i][m_fibers[i].size()-1]))
            {
                cerr << "Extrémité identique" << endl;
                new_fiber.moveExtremities(Vector3f(ecart, ecart, ecart));
            }
            i++;
        }
        if (same)
            cerr << "Fibre en double trouvée" << endl;
        else
        {
            m_nbOfPoints+=new_fiber.size();
            m_fibers.push_back(new_fiber);
        }
    }
#else
    {
        //new_fiber.resample(100);
        m_nbOfPoints+=new_fiber.size();
//        float value=0;
//        for (unsigned int j=0; j<new_fiber.size()-1; j++)
//        {
//            float ptValue = m_fuzzySet->getValueFromSegment(new_fiber[j], new_fiber[j+1]);
//            value += (new_fiber[j+1]-new_fiber[j]).norm()/new_fiber.length() * ptValue;
//        }
//        float extremity = exp(-min(m_maskEndpoints1->getDistanceToMask(new_fiber[0])+m_maskEndpoints2->getDistanceToMask(new_fiber[new_fiber.size()-1]),
//                m_maskEndpoints1->getDistanceToMask(new_fiber[new_fiber.size()-1])+m_maskEndpoints2->getDistanceToMask(new_fiber[0]))/m_gaussianDistance);
//        value*=extremity;
//        new_fiber.setValue((float)value/255);

//        value=0;
//        float ptValue = m_fuzzySet->getValueFromVolume(new_fiber.getProfileTransform(0).center, new_fiber.getProfileTransform(0), new_fiber.getProfileTransform(1), new_fiber.getProfileTransform(2).center);
//        value += (new_fiber.getProfileTransform(1).center -new_fiber.getProfileTransform(0).center).norm()/new_fiber.cylinderLength() * ptValue;
//        for (unsigned int j=1; j<new_fiber.size()-2; j++)
//        {
//            ptValue = m_fuzzySet->getValueFromVolume(new_fiber.getProfileTransform(j-1).center, new_fiber.getProfileTransform(j), new_fiber.getProfileTransform(j+1), new_fiber.getProfileTransform(j+2).center);
//            value += (new_fiber.getProfileTransform(j+1).center-new_fiber.getProfileTransform(j).center).norm()/new_fiber.cylinderLength() * ptValue;
//        }
//        extremity = exp(-min(m_maskEndpoints1->getDistanceToMaskFromSurface(new_fiber.getProfileTransform(0), new_fiber.getProfileTransform(1).center-new_fiber.getProfileTransform(0).center)+
//                                   m_maskEndpoints2->getDistanceToMaskFromSurface(new_fiber.getProfileTransform(new_fiber.size()-1), new_fiber.getProfileTransform(new_fiber.size()-2).center-new_fiber.getProfileTransform(new_fiber.size()-1).center),
//                                   m_maskEndpoints1->getDistanceToMaskFromSurface(new_fiber.getProfileTransform(new_fiber.size()-1), new_fiber.getProfileTransform(new_fiber.size()-2).center-new_fiber.getProfileTransform(new_fiber.size()-1).center)+
//                                   m_maskEndpoints2->getDistanceToMaskFromSurface(new_fiber.getProfileTransform(0), new_fiber.getProfileTransform(1).center-new_fiber.getProfileTransform(0).center))/m_gaussianDistance);
//        value*=extremity;
//        ptValue = m_fuzzySet->getValueFromVolume(new_fiber.getProfileTransform(new_fiber.size()-3).center, new_fiber.getProfileTransform(new_fiber.size()-2), new_fiber.getProfileTransform(new_fiber.size()-1), new_fiber.getProfileTransform(new_fiber.size()-1).center);
//        value += (new_fiber.getProfileTransform(new_fiber.size()-1).center-new_fiber.getProfileTransform(new_fiber.size()-2).center).norm()/new_fiber.cylinderLength() * ptValue;
//        new_fiber.setCylinderValue((float)value/255);
        //cout << new_fiber.getValue() << " " << new_fiber.getValueOfCylinder() << endl;

        m_fibers.push_back(new_fiber);
    }
#endif
}

void Bundle::computeFiberFuzzy(unsigned int i)
{
    float value=0;
    float cylinderValue=0;
    float ptValue = m_fuzzySet->getValueFromVolume(m_fibers[i].getProfileTransform(0).center, m_fibers[i].getProfileTransform(0), m_fibers[i].getProfileTransform(1), m_fibers[i].getProfileTransform(2).center);
    cylinderValue += (m_fibers[i].getProfileTransform(1).center -m_fibers[i].getProfileTransform(0).center).norm()/m_fibers[i].cylinderLength() * ptValue;
    ptValue = m_fuzzySet->getValueFromSegment(m_fibers[i][0], m_fibers[i][1]);
    value += (m_fibers[i][1]-m_fibers[i][0]).norm()/m_fibers[i].length() * ptValue;
#pragma omp parallel for reduction (+:value, cylinderValue)
    for (unsigned int j=1; j<m_fibers[i].size()-2; j++)
    {
        float ptValueFiber = m_fuzzySet->getValueFromSegment(m_fibers[i][j], m_fibers[i][j+1]);
        value += (m_fibers[i][j+1]-m_fibers[i][j]).norm()/m_fibers[i].length() * ptValueFiber;
        float ptValueCylinder = m_fuzzySet->getValueFromVolume(m_fibers[i].getProfileTransform(j-1).center, m_fibers[i].getProfileTransform(j), m_fibers[i].getProfileTransform(j+1), m_fibers[i].getProfileTransform(j+2).center);
        cylinderValue += (m_fibers[i].getProfileTransform(j+1).center-m_fibers[i].getProfileTransform(j).center).norm()/m_fibers[i].cylinderLength() * ptValueCylinder;
    }
    ptValue = m_fuzzySet->getValueFromSegment(m_fibers[i][m_fibers[i].size()-2], m_fibers[i][m_fibers[i].size()-1]);
    value += (m_fibers[i][m_fibers[i].size()-1]-m_fibers[i][m_fibers[i].size()-2]).norm()/m_fibers[i].length() * ptValue;
    float extremityFiber = exp(-min(m_maskEndpoints1->getDistanceToMask(m_fibers[i][0])+m_maskEndpoints2->getDistanceToMask(m_fibers[i][m_fibers[i].size()-1]),
            m_maskEndpoints1->getDistanceToMask(m_fibers[i][m_fibers[i].size()-1])+m_maskEndpoints2->getDistanceToMask(m_fibers[i][0]))/m_gaussianDistance);
    value*=extremityFiber;
    m_fibers[i].setValue((float)value/255);

    float extremityCylinder = exp(-min(m_maskEndpoints1->getDistanceToMaskFromSurface(m_fibers[i].getProfileTransform(0), m_fibers[i].getProfileTransform(1).center-m_fibers[i].getProfileTransform(0).center)+
                               m_maskEndpoints2->getDistanceToMaskFromSurface(m_fibers[i].getProfileTransform(m_fibers[i].size()-1), m_fibers[i].getProfileTransform(m_fibers[i].size()-2).center-m_fibers[i].getProfileTransform(m_fibers[i].size()-1).center),
                               m_maskEndpoints1->getDistanceToMaskFromSurface(m_fibers[i].getProfileTransform(m_fibers[i].size()-1), m_fibers[i].getProfileTransform(m_fibers[i].size()-2).center-m_fibers[i].getProfileTransform(m_fibers[i].size()-1).center)+
                               m_maskEndpoints2->getDistanceToMaskFromSurface(m_fibers[i].getProfileTransform(0), m_fibers[i].getProfileTransform(1).center-m_fibers[i].getProfileTransform(0).center))/m_gaussianDistance);
    cylinderValue*=extremityCylinder;
    ptValue = m_fuzzySet->getValueFromVolume(m_fibers[i].getProfileTransform(m_fibers[i].size()-3).center, m_fibers[i].getProfileTransform(m_fibers[i].size()-2), m_fibers[i].getProfileTransform(m_fibers[i].size()-1), m_fibers[i].getProfileTransform(m_fibers[i].size()-1).center);
    cylinderValue += (m_fibers[i].getProfileTransform(m_fibers[i].size()-1).center-m_fibers[i].getProfileTransform(m_fibers[i].size()-2).center).norm()/m_fibers[i].cylinderLength() * ptValue;
    m_fibers[i].setCylinderValue((float)cylinderValue/255);
}

void Bundle::computeFiberFuzzyNotParallel(unsigned int i)
{
    float value=0;
    for (unsigned int j=0; j<m_fibers[i].size()-1; j++)
    {
        float ptValue = m_fuzzySet->getValueFromSegment(m_fibers[i][j], m_fibers[i][j+1]);
        value += (m_fibers[i][j+1]-m_fibers[i][j]).norm()/m_fibers[i].length() * ptValue;
    }
    float extremity = exp(-min(m_maskEndpoints1->getDistanceToMask(m_fibers[i][0])+m_maskEndpoints2->getDistanceToMask(m_fibers[i][m_fibers[i].size()-1]),
            m_maskEndpoints1->getDistanceToMask(m_fibers[i][m_fibers[i].size()-1])+m_maskEndpoints2->getDistanceToMask(m_fibers[i][0]))/m_gaussianDistance);
    value*=extremity;
    m_fibers[i].setValue((float)value/255);

    value=0;
    float ptValue = m_fuzzySet->getValueFromVolume(m_fibers[i].getProfileTransform(0).center, m_fibers[i].getProfileTransform(0), m_fibers[i].getProfileTransform(1), m_fibers[i].getProfileTransform(2).center);
    value += (m_fibers[i].getProfileTransform(1).center -m_fibers[i].getProfileTransform(0).center).norm()/m_fibers[i].cylinderLength() * ptValue;
    for (unsigned int j=1; j<m_fibers[i].size()-2; j++)
    {
        ptValue = m_fuzzySet->getValueFromVolume(m_fibers[i].getProfileTransform(j-1).center, m_fibers[i].getProfileTransform(j), m_fibers[i].getProfileTransform(j+1), m_fibers[i].getProfileTransform(j+2).center);
        value += (m_fibers[i].getProfileTransform(j+1).center-m_fibers[i].getProfileTransform(j).center).norm()/m_fibers[i].cylinderLength() * ptValue;
    }
    extremity = exp(-min(m_maskEndpoints1->getDistanceToMaskFromSurface(m_fibers[i].getProfileTransform(0), m_fibers[i].getProfileTransform(1).center-m_fibers[i].getProfileTransform(0).center)+
                               m_maskEndpoints2->getDistanceToMaskFromSurface(m_fibers[i].getProfileTransform(m_fibers[i].size()-1), m_fibers[i].getProfileTransform(m_fibers[i].size()-2).center-m_fibers[i].getProfileTransform(m_fibers[i].size()-1).center),
                               m_maskEndpoints1->getDistanceToMaskFromSurface(m_fibers[i].getProfileTransform(m_fibers[i].size()-1), m_fibers[i].getProfileTransform(m_fibers[i].size()-2).center-m_fibers[i].getProfileTransform(m_fibers[i].size()-1).center)+
                               m_maskEndpoints2->getDistanceToMaskFromSurface(m_fibers[i].getProfileTransform(0), m_fibers[i].getProfileTransform(1).center-m_fibers[i].getProfileTransform(0).center))/m_gaussianDistance);
    value*=extremity;
    ptValue = m_fuzzySet->getValueFromVolume(m_fibers[i].getProfileTransform(m_fibers[i].size()-3).center, m_fibers[i].getProfileTransform(m_fibers[i].size()-2), m_fibers[i].getProfileTransform(m_fibers[i].size()-1), m_fibers[i].getProfileTransform(m_fibers[i].size()-1).center);
    value += (m_fibers[i].getProfileTransform(m_fibers[i].size()-1).center-m_fibers[i].getProfileTransform(m_fibers[i].size()-2).center).norm()/m_fibers[i].cylinderLength() * ptValue;
    m_fibers[i].setCylinderValue((float)value/255);
}

void Bundle::computeFuzzy()
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        computeFiberFuzzyNotParallel(i);
    }
}

unsigned int Bundle::addFiberSize(Fiber &new_fiber)
{
    addFiber(new_fiber);
    return (unsigned int)m_fibers.size();
}

void Bundle::computeDistances()
{
    unsigned int size = m_fibers.size();
    vector<priority_queue<paire, vector<paire>, oppositeComparisonPaires > > paires(omp_get_max_threads());
#pragma omp parallel for schedule(dynamic, CHUNK)
    for (unsigned int i=0; i<size; i++)
    {
        for (set<unsigned int>::iterator ite=m_fibers[i].m_neighbours.begin(); ite!=m_fibers[i].m_neighbours.end(); ++ite)
        {
            unsigned int j=*ite;
            if (i<j)//To compute each pair only once
            {
                paire a;
                a.fib1=i;
                a.dist=m_metric->metriqueTest(m_fibers[i], m_fibers[j]);
                a.fib2=j;
                paires[omp_get_thread_num()].push(a);
            }
        }
    }

    for (unsigned int i=0; i<paires.size(); i++)
        while(!paires[i].empty())
        {
            m_paires.push(paires[i].top());
            paires[i].pop();
        }
}

///
/// MultiScale
///

void Bundle::computeMultiScale()
{
    //Statistiques
#if (STATS)
    float valenceMax=std::numeric_limits<float>::min(), valenceMin=std::numeric_limits<float>::max(), valenceMoyenne=0;
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        float value=m_fibers[i].m_neighbours.size();
        valenceMoyenne+=value;
        if (value<valenceMin)
            valenceMin=value;
        if (value>valenceMax)
            valenceMax=value;
    }

    vector<float>moyennesValence;
    float maxValence=valenceMax;

    fstream file;
    file.open("Stats.txt", ios::out);
    file << valenceMin << " " << valenceMoyenne/m_fibers.size() << " " << valenceMax << endl;
    unsigned int nbValidElements=m_fibers.size();
#endif

    m_sizeWhenComputingMultiScale=m_fibers.size();
    unsigned int count=0;
    while (!m_paires.empty())
    {
        //Sélection de la paire de distance minimale
        paire top;
        float angle=-1;//cos(M_PI_2);
        while (!m_paires.empty() && angle<m_limitAngle)
        {
            while (!m_paires.empty() && (!m_fibers[m_paires.top().fib1].is_valid() || !m_fibers[m_paires.top().fib2].is_valid()))
                m_paires.pop();
            if (!m_paires.empty())
            {
                top=m_paires.top();
                float x=m_fibers[top.fib1].getSelfScalarProduct(7);
                float y=m_fibers[top.fib2].getSelfScalarProduct(7);
//                float wc=scalarProductWeightedCurrentsParallel(m_fibers[top.fib1], m_fibers[top.fib2], 6,6,8);
                float wc = fabs(top.dist);
                angle = wc/sqrt(x*y);
                //cout << wc << " " << x << " " << y << " " << angle << m_limitAngle << endl;
                if (angle<m_limitAngle)
                    m_paires.pop();
            }
        }
        if (m_paires.empty())
            break;

        Fiber newFiber=m_metric->mergeOfFibersV2(m_fibers[top.fib1], m_fibers[top.fib2]);
        //newFiber.resample(100);

        m_fibers[top.fib1].setValid(false);
        m_fibers[top.fib2].setValid(false);
        //Suppression des anciennes fibres
        //Ajout des voisins des deux fibres à la nouvelle fibre
        for (set<unsigned int>::iterator ite=m_fibers[top.fib1].m_neighbours.begin(); ite!=m_fibers[top.fib1].m_neighbours.end(); ++ite)
            if (m_fibers[*ite].is_valid())
                newFiber.m_neighbours.insert(*ite);
        for (set<unsigned int>::iterator ite=m_fibers[top.fib2].m_neighbours.begin(); ite!=m_fibers[top.fib2].m_neighbours.end(); ++ite)
            if (m_fibers[*ite].is_valid())
                newFiber.m_neighbours.insert(*ite);

        //Trying to remove previous neighbors for memory saving
        m_fibers[top.fib1].m_neighbours.clear();
        m_fibers[top.fib2].m_neighbours.clear();
//        m_fibers[top.fib1].clearCentersAndTangents();
//        m_fibers[top.fib2].clearCentersAndTangents();

        newFiber.computeCentersAndTangents();
        unsigned int nbFibres=addFiberSize(newFiber);
        //Computation of the fuzzy value of the new fiber
        if (thereIsAFuzzySet())
            computeFiberFuzzy(nbFibres-1);
        triplet t;
        t.fib1=top.fib1;
        t.fib2=top.fib2;
        t.newFib=nbFibres-1;
        m_changes.push_back(t);

        //Calcul des nouvelles paires
        unsigned int nbNeighbours=newFiber.m_neighbours.size();
        vector<vector<paire> > pairesThread(omp_get_max_threads());
        vector<unsigned int> voisins(nbNeighbours);
        unsigned int compte=0;
        for (set<unsigned int>::iterator ite=newFiber.m_neighbours.begin(); ite!=newFiber.m_neighbours.end(); ++ite)
        {
            voisins[compte]=*ite;
            compte++;
        }
        //cout << "Nb voisins: " << nbNeighbours << endl;
#pragma omp parallel for
        for (unsigned int i=0; i<nbNeighbours; i++)
        {
            unsigned int j=voisins[i];
            if (m_fibers[j].is_valid())
            {
                paire a;
                a.dist=m_metric->metriqueTest(m_fibers[nbFibres-1], m_fibers[j]);
                a.fib1=nbFibres-1;
                a.fib2=j;
#if (STATS)
                m_fibers[j].m_neighbours.erase(top.fib1);
                m_fibers[j].m_neighbours.erase(top.fib2);
#endif
                m_fibers[j].m_neighbours.insert(nbFibres-1);
                pairesThread[omp_get_thread_num()].push_back(a);
            }
            else
            {
                newFiber.m_neighbours.erase(j);
            }
        }
        for (int i=0; i<omp_get_max_threads(); i++)
            for (unsigned int j=0; j<pairesThread[i].size(); j++)
                m_paires.push(pairesThread[i][j]);

        //Statistiques
#if (STATS)
        nbValidElements--;
        vector<float> valMoyenne(omp_get_max_threads());
        vector<float> valMax(omp_get_max_threads());
        vector<float> valMin(omp_get_max_threads());
#pragma omp parallel
        {
            valMoyenne[omp_get_thread_num()]=0;
            valMax[omp_get_thread_num()]=std::numeric_limits<float>::min();;
            valMin[omp_get_thread_num()]=std::numeric_limits<float>::max();
        }
        //valenceMoyenne=0; valenceMin=std::numeric_limits<float>::max();
        //valenceMax=std::numeric_limits<float>::min();
#pragma omp parallel for
        for (unsigned int i=0; i<m_fibers.size(); i++)
        {
            if (m_fibers[i].is_valid())
            {
                float value=m_fibers[i].m_neighbours.size();
                valMoyenne[omp_get_thread_num()]+=value;
                if (value<valMin[omp_get_thread_num()])
                    valMin[omp_get_thread_num()]=value;
                if (value>valMax[omp_get_thread_num()])
                    valMax[omp_get_thread_num()]=value;
            }
        }
        valenceMoyenne=0; valenceMin=std::numeric_limits<float>::max(); valenceMax=std::numeric_limits<float>::min();
        for (unsigned int i=0; i<omp_get_max_threads(); i++)
        {
            valenceMoyenne+=valMoyenne[i];
            if (valMin[i]<valenceMin) valenceMin=valMin[i];
            if (valMax[i]>valenceMax) valenceMax=valMax[i];
        }
        moyennesValence.push_back(valenceMoyenne/nbValidElements);
        if (valenceMax>maxValence)
            maxValence=valenceMax;
        file << valenceMin << " " << valenceMoyenne/nbValidElements << " " << valenceMax << endl;
#endif

        //Récupération des paires de fibres minimisant la distance
        if(count%101==0)
            advanceBarCout((float)(m_fibers.size()-m_sizeWhenComputingMultiScale)/m_sizeWhenComputingMultiScale*100);
        count++;
    }
#if (STATS)
    file.close();
    //    file2.close();
    //    file.close();
    float moy=0;
    for (unsigned int i=0; i<moyennesValence.size(); i++)
        moy+=moyennesValence[i];
    moy/=count;
    float esp=0;
    for (unsigned int i=0; i<moyennesValence.size(); i++)
        esp+=pow(moyennesValence[i]-moy,2);
    esp=sqrt(esp/count);

    cout << "Valence max : " << maxValence << endl;
    cout << "Valence moyenne : " << moy  << "+-" << esp << endl;
#endif

    cout <<endl << "Nb de contractions : " << count << endl;
    m_positionFibreHierarchie=m_changes.size()-1;
    m_fin=m_sizeWhenComputingMultiScale;
    m_percentage=100;
    setMultiScale(0.f);

    //Record of the calculation that has just be done
    if (!recordBundle())
    {
        cerr << "------------------ No record of the bundle has been made -----------------" << endl;
    }
}

void Bundle::reOrienteFibers(unsigned int reference)
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
        m_fibers[i].PutInSameDirection(m_fibers[reference]);
}

//MAJ de la validité ou non des fibres du bundle
void Bundle::setMultiScale(int diffNbFibres)
{
    int debut=m_fin;
    m_fin=debut+diffNbFibres;
    if (m_fin<0) m_fin=0;
    if (m_fin>(int)m_changes.size()) m_fin=m_changes.size();
    changeDisplayMultiScale(debut, m_fin);
    m_percentage=m_fin*100.f/float(m_changes.size());
}
//Le pourcentage est à 0 lorsque les fibres de départ sont affichées et à 100 lorsque le niveau de multi-résolution est au plus bas
void Bundle::setMultiScale(float percentage)
{
    if (percentage==m_percentage)
        return;
    else if (percentage<0 || percentage >100)
    {
        cerr << "Erreur, le pourcentage " << percentage << " n'est pas valide" << endl;
        return;
    }
    int debut, fin;
    debut=m_percentage*float(m_changes.size())/100.f;
    fin=percentage*float(m_changes.size())/100.f;
    changeDisplayMultiScale(debut, fin);
    m_percentage=percentage;
}

void Bundle::changeDisplayMultiScale(int debut, int fin)
{
    if (debut<fin)
    {
        for (int i=debut; i<fin; i++)
        {
            m_fibers[m_changes[i].fib1].setValid(false);
            m_fibers[m_changes[i].fib2].setValid(false);
            m_fibers[m_changes[i].newFib].setValid(true);
        }
        m_lastChange=fin-1;
    }
    else
    {
        for (int i=debut-1; i>=fin; i--)
        {
            m_fibers[m_changes[i].fib1].setValid(true);
            m_fibers[m_changes[i].fib2].setValid(true);
            m_fibers[m_changes[i].newFib].setValid(false);
        }
        if (fin!=(int)m_changes.size())
            m_lastChange=fin;
    }
    m_fin=fin;
}

//Fonction permettant de marquer les fibres à l'origine de la première fibre parmi celles visibles à l'écran
void Bundle::getOriginOfHierarchy()
{
    //Remonter la hiérarchie
    triplet lastPaire=m_changes[m_positionFibreHierarchie];
    unsigned int fib1=lastPaire.fib1;
    unsigned int fib2=lastPaire.fib2;
    forward_list<int> fibresIntermediaires;
    int taille=0;
    if (fib1>=m_sizeWhenComputingMultiScale) {fibresIntermediaires.push_front(fib1); taille++;}
    if (fib2>=m_sizeWhenComputingMultiScale) {fibresIntermediaires.push_front(fib2); taille++;}
    int position=m_positionFibreHierarchie-1;
    while(!fibresIntermediaires.empty() && position>=0)
    {
        triplet currentPaire=m_changes[position];
        auto iter=fibresIntermediaires.begin();
        for (iter=fibresIntermediaires.begin(); iter!=fibresIntermediaires.end(); ++iter)
        {
            if (currentPaire.newFib==(unsigned int)*iter)
            {
                break;
            }
        }
        if (iter!=fibresIntermediaires.end() && currentPaire.newFib==(unsigned int)*iter)
        {
            fibresIntermediaires.remove(*iter);
            taille--;
            if (currentPaire.fib1>=m_sizeWhenComputingMultiScale)
            {
                fibresIntermediaires.push_front(currentPaire.fib1);
                taille++;
            }
            else
            {
                m_fibers[currentPaire.fib1].setValid(true);
                m_fibers[currentPaire.fib1].setValue(-1);//Fibre de départ signalée
            }
            if (currentPaire.fib2>=m_sizeWhenComputingMultiScale)
            {
                fibresIntermediaires.push_front(currentPaire.fib2);
                taille++;
            }
            else
            {
                m_fibers[currentPaire.fib1].setValid(true);
                m_fibers[currentPaire.fib1].setValue(-1);//Fibre de départ signalée
            }
        }
        position--;
    }
}

void Bundle::updateNeighbours()
{
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        for (set<unsigned int>::iterator ite=m_fibers[i].m_neighbours.begin(); ite!=m_fibers[i].m_neighbours.end(); ++ite)
        {
            unsigned int j=*ite;
            m_fibers[j].m_neighbours.insert(i);
        }
    }
}

void Bundle::colorSubBundles()
{
    int marqueur=0;
    vector<bool> marqued(m_fibers.size(), false);
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        if (marqued[i])
            continue;
        queue<int> nodeQueue;
        nodeQueue.push(i);
        marqueur++;
        marqued[i]=true;
        m_fibers[i].setValue(marqueur);
        while(!nodeQueue.empty())
        {
            int s=nodeQueue.front();
            nodeQueue.pop();
            for (set<unsigned int>::iterator ite=m_fibers[s].m_neighbours.begin(); ite!=m_fibers[s].m_neighbours.end(); ++ite)
            {
                unsigned int voisin = *ite;
                if (!marqued[voisin])
                {
                    nodeQueue.push(voisin);
                    marqued[voisin]=true;
                    m_fibers[voisin].setValue(marqueur);
                }
            }
        }
    }
}


void Bundle::resample(unsigned int resampleValue)
{
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        //cout << i << " " << m_fibers[i].size() <<  endl;
        //if (m_fibers[i].size()<resampleValue)
        m_fibers[i].resample(resampleValue);
    }
}

void Bundle::computeDelaunay()
{
    tetgenio test;
    test.initialize();
    //test.numberofpoints=m_fibers.size()*2;
    test.numberofpoints=(int)getCurrentNberOfFibers()*2;
    test.pointlist=new REAL[test.numberofpoints*3];
    unsigned int counter=0;
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        if (m_fibers[i].is_valid())
        {
            int positionFibre=0;
            for (unsigned int j=0; j<3; j++)
                test.pointlist[counter*6+j]=m_fibers[i].getPoint(positionFibre,j);
            positionFibre=m_fibers[i].size()-1;
            for (unsigned int j=0; j<3; j++)
                test.pointlist[counter*6+j+3]=m_fibers[i].getPoint(positionFibre,j);
            counter++;
        }
    }
    tetgenio output;
    output.initialize();
    char paramTetgen[] = "RNIB";
    tetrahedralize(paramTetgen, &test, &output, NULL, NULL);
    cout << "Number of tetrahedron " << output.numberoftetrahedra << endl;
    for (unsigned int i=0; i<(unsigned int)output.numberoftetrahedra; i++)
    {
        vector<unsigned int> tetra(4);
        for (unsigned int j=0; j<4; j++)
            tetra[j]=output.tetrahedronlist[i*4+j]/2;//%m_fibers.size();
        for (unsigned int j=0; j<4; j++)
            for (unsigned int k=0; k<3; k++)
                m_fibers[tetra[j]].m_neighbours.insert(tetra[(j+k)%4]);
    }
    //Tetrahedralization on the middle points
//    tetgenio test2;
//    test2.initialize();
//    test2.numberofpoints=(int)getCurrentNberOfFibers();
//    test2.pointlist=new REAL[test2.numberofpoints*3];
//    counter=0;
//    for (unsigned int i=0; i<m_fibers.size(); i++)
//    {
//        if (m_fibers[i].is_valid())
//        {
//            int positionFibre=m_fibers[i].size()/2;
//            for (unsigned int j=0; j<3; j++)
//                test2.pointlist[counter*3+j]=m_fibers[i].getPoint(positionFibre,j);
//            counter++;
//        }
//    }
//    tetgenio output2;
//    output2.initialize();
//    tetrahedralize(paramTetgen, &test2, &output2, NULL, NULL);
//    cout << "Number of tetrahedron " << output2.numberoftetrahedra << endl;
//    for (unsigned int i=0; i<(unsigned int)output2.numberoftetrahedra; i++)
//    {
//        vector<unsigned int> tetra(4);
//        for (unsigned int j=0; j<4; j++)
//            tetra[j]=output2.tetrahedronlist[i*4+j]/2;//%m_fibers.size();
//        for (unsigned int j=0; j<4; j++)
//            for (unsigned int k=0; k<3; k++)
//                m_fibers[tetra[j]].m_neighbours.insert(tetra[(j+k)%4]);
//    }

    for (unsigned int i=0; i<m_fibers.size(); i++)
        m_fibers[i].m_neighbours.erase(i);
}

void Bundle::makeAllFibersNeighbors()
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
        for (unsigned j=0; j<m_fibers.size(); j++)
            if (j != i)
                m_fibers[i].m_neighbours.insert(j);
}

vector<Vector3f> Bundle::drawNeighbors()
{
    vector<Vector3f> edges;
    for (unsigned int i=0; i<m_fibers.size(); i++)
        if (m_fibers[i].is_valid())
            for (set<unsigned int>::iterator ite=m_fibers[i].m_neighbours.begin(); ite!=m_fibers[i].m_neighbours.end(); ++ite)
                if (m_fibers[*ite].is_valid())
                {
                    edges.push_back(m_fibers[i][0]);
                    edges.push_back(m_fibers[*ite][0]);
                }
    return edges;
}

vector<unsigned int> Bundle::drawNeighborsFibers()
{
    vector<unsigned int> neighbors;
    unsigned int i=(float)rand()/RAND_MAX*(float)getCurrentNberOfFibers();
    unsigned int j=0;
    unsigned int k=0;
    while (j!=i)
    {
        if (m_fibers[k].is_valid())
            j++;
        k++;
    }
    i=j;
    neighbors.push_back(i);
    for (set<unsigned int>::iterator ite=m_fibers[i].m_neighbours.begin(); ite!=m_fibers[i].m_neighbours.end(); ++ite)
        if (m_fibers[*ite].is_valid())
        {
            neighbors.push_back(*ite);
        }
    return neighbors;
}

vector<Vector3f> Bundle::drawOneNeighbour()
{
    vector<Vector3f> edges;
    unsigned int i=(float)rand()/RAND_MAX*(float)getCurrentNberOfFibers();
    unsigned int j=0;
    unsigned int k=0;
    while (j!=i)
    {
        if (m_fibers[k].is_valid())
            j++;
        k++;
    }
    i=j;
    for (set<unsigned int>::iterator ite=m_fibers[i].m_neighbours.begin(); ite!=m_fibers[i].m_neighbours.end(); ++ite)
        if (m_fibers[*ite].is_valid())
        {
            edges.push_back(m_fibers[i][0]);
            edges.push_back(m_fibers[*ite][0]);
        }
    return edges;
}

void Bundle::printNeighboursInFile()
{
    fstream file;
    file.open("Neighbours.txt", ios_base::out);
    for (unsigned int i=0; i<m_fibers.size(); i++)
        file << i << " " << m_fibers[i].m_neighbours.size() << endl;
    file.close();
}

bool Bundle::recordBundle()
{
    string filename = m_pwidget->getFilename();
    if (filename.empty())
        return false;
    fstream file;
    file.open(filename,  ios_base::out | ios_base::binary);
    file << this->size() << endl;
    file << m_sizeWhenComputingMultiScale << endl;
    Vector3f point;
    Vector3f color;
    Transfo transfo;
    unsigned int nbPoints;
    unsigned int hierarchy;
    unsigned int nbNeighbours;
    unsigned int neighbour;
    for (unsigned int i=0; i<this->size(); i++)
    {
        nbPoints=m_fibers[i].size();
        file.write(reinterpret_cast<const char *>(&nbPoints), sizeof(unsigned int));
        color=m_fibers[i].getColor();
        file.write(reinterpret_cast<const char *>(&color), sizeof(color));
        hierarchy=m_fibers[i].getNbFibHierarchy();
        file.write(reinterpret_cast<const char *>(&hierarchy), sizeof(unsigned int));
        for (unsigned int j=0; j<m_fibers[i].size(); j++)
        {
            point=m_fibers[i][j];
            file.write(reinterpret_cast<const char *>(&point), sizeof(point));
            transfo=m_fibers[i].getProfileTransform(j);
            file.write(reinterpret_cast<const char *>(&transfo), sizeof(Transfo));
        }
        //Save of the neighbours
        nbNeighbours=m_fibers[i].m_neighbours.size();
        file.write(reinterpret_cast<const char *>(&nbNeighbours), sizeof(unsigned int));
        for (set<unsigned int>::iterator ite=m_fibers[i].m_neighbours.begin(); ite!=m_fibers[i].m_neighbours.end(); ++ite)
        {
            neighbour=*ite;
            file.write(reinterpret_cast<const char *>(&neighbour), sizeof(unsigned int));
        }
    }
    //Save of changes -> the order of the merges
    file << m_changes.size() << endl;
    for (unsigned int i=0; i<m_changes.size(); i++)
    {
        triplet change=m_changes[i];
        file.write(reinterpret_cast<const char *>(&change), sizeof(triplet));
    }
    file.close();
    cout << "Bundle saved : " << filename << endl;
    return true;
}

bool Bundle::loadBundle(string filename)
{
    fstream file;
    file.open(filename, ios_base::in | ios_base::binary);
    unsigned int nbFibers;
    file >> nbFibers;
    cout << nbFibers << " fibers found in the file " << filename << endl;
    file >> m_sizeWhenComputingMultiScale;
    char machin='a';
    while (machin!='\n'){file.get(machin);}
    m_fibers.clear();
    unsigned int nbPoints;
    unsigned int hierarchy;
    unsigned int nbNeighbours;
    unsigned int neighbour;
    Vector3f point;
    Vector3f color;
    MatrixXf Points;
    vector<Transfo> listTransfos;
    Transfo transfo;
    for (unsigned int i=0; i<nbFibers; i++)
    {
        if (!file.read(reinterpret_cast<char *>(&nbPoints), sizeof(unsigned int)))
            cerr << "Error while reading file" << endl;
        Points.resize(nbPoints, 3);
        listTransfos.resize(nbPoints);
        if (!file.read(reinterpret_cast<char *>(&color), sizeof(color)))
            cerr << "Error while reading file" << endl;
        if (!file.read(reinterpret_cast<char *>(&hierarchy), sizeof(unsigned int)))
            cerr << "Error while reading file" << endl;
        for (unsigned int j=0; j<nbPoints; j++)
        {
            if (!file.read(reinterpret_cast<char *>(&point), sizeof(point)))
                cerr << "Error while reading file" << endl;
            Points.row(j)=point;
            if (!file.read(reinterpret_cast<char *>(&transfo), sizeof(Transfo)))
                cerr << "Error while reading file" << endl;
            if (transfo.a<MIN_RADIUS)
                transfo.a=MIN_RADIUS;
            if (transfo.b<MIN_RADIUS)
                transfo.b=MIN_RADIUS;
            listTransfos[j]=transfo;
        }
        Fiber newFib(Points, nbPoints);
        if (!file.read(reinterpret_cast<char *>(&nbNeighbours), sizeof(unsigned int)))
            cerr << "Error while reading file" << endl;
        for (unsigned int j=0; j<nbNeighbours; j++)
        {
            if (!file.read(reinterpret_cast<char *>(&neighbour), sizeof(unsigned int)))
                cerr << "Error while reading file" << endl;
            newFib.m_neighbours.insert(neighbour);
        }
        color*=255;
        newFib.setColor(color(0), color(1), color(2));
        newFib.addFibHierarchy(hierarchy-1);
        newFib.setProfileTransform(listTransfos);
//        if(i<nbFibers-1)
//            newFib.setValid(false);
//        else
//            newFib.setValid(true);
        this->addFiber(newFib);
    }
    //Get back order of changes
    unsigned int nbOfChanges;
    file >> nbOfChanges;
    machin='a';
    while (machin!='\n'){file.get(machin);}
    m_changes.resize(nbOfChanges);
    triplet change;
    for (unsigned int i=0; i<nbOfChanges; i++)
    {
        if (!file.read(reinterpret_cast<char *>(&change), sizeof(triplet)))
            cerr << "Error while reading file" << endl;
        m_changes[i]=change;
    }
    file.close();
    cout << "End of file" << endl;
    for (unsigned int i=0; i<this->size(); i++)
        m_fibers[i].addToHierarchy(i);
    setMultiScale(0.f);
    return true;
}


bool Bundle::saveCurrentFibersAsVTK(string filename)
{
    unsigned int nbPoints = getCurrentNberOfPoints();
    unsigned int nbFibers = getCurrentNberOfFibers();

    fstream file;
    file.open(filename, ios_base::out);

    file << "# vtk DataFile Version 7.1" << endl;
    file << "Segmentation obtained using Fiberreductor, made by Corentin Mercier" << endl;
    file << "ASCII" << endl;
    file << "DATASET POLYDATA" << endl;
    file << "POINTS " << nbPoints <<  " float" << endl;

    for (unsigned int i=0; i<size(); i++)
    {
        if (m_fibers[i].is_valid())
        {
            for (unsigned int j=0; j<m_fibers[i].size(); j++)
                file << m_fibers[i][j].x() << " " << m_fibers[i][j].y() << " " << m_fibers[i][j].z() << endl;
        }
    }

    file << "LINES " << nbFibers << " " << nbPoints+nbFibers << endl;
    unsigned long int counter = 0;
    for (unsigned int i=0; i<size(); i++)
    {
        if (m_fibers[i].is_valid())
        {
            file << m_fibers[i].size() << " ";
            for (unsigned int j=0; j<m_fibers[i].size(); j++)
            {
                file << counter << " ";
                counter++;
            }
            file << endl;
        }
    }
    file << endl;
    file << "POINT_DATA " << counter << endl;
    file << "SCALARS acs float 1" << endl;
    file << "LOOKUP_TABLE default" << endl;
    for (unsigned int i=0; i<size(); i++)
    {
        if (m_fibers[i].is_valid())
        {
            for (unsigned int j=0; j<m_fibers[i].size(); j++)
                file << m_fibers[i].getValue() << endl;
        }
    }
    file.close();
    return true;
}


void Bundle::set_widget(WidgetOpenGL* widget_param)
{
    m_pwidget = widget_param;
}

void Bundle::geometricSmoothing()
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        m_fibers[i].geometricSmoothing();
    }
}

void Bundle::pointsSmoothing()
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        m_fibers[i].pointsSmoothing();
    }
}

unsigned int Bundle::getCurrentNberOfPoints()
{
    unsigned int value=0;
#pragma omp parallel for reduction (+:value)
    for (unsigned int i=0; i<m_fibers.size(); i++)
        if (m_fibers[i].is_valid())
            value=value+m_fibers[i].size();
    return value;
}

unsigned int Bundle::getCurrentNberOfFibers()
{
    unsigned int value=0;
#pragma omp parallel for reduction (+:value)
    for (unsigned int i=0; i<m_fibers.size(); i++)
        if (m_fibers[i].is_valid())
            value++;
    return value;
}

//void Bundle::cutFibers()
//{
//    //unsigned int toDelete = 0;
//    vector<Fiber> fibers;
//    fibers.reserve(size());
//    m_nbOfPoints = 0;
//    for (unsigned int i=0; i<size(); i++)
//    {
//        Fiber newFib = m_fibers[i].cut(m_surface);
//        if (newFib.is_valid() && newFib.size() > 2)
//        {
//            fibers.push_back(newFib);
//            m_nbOfPoints+=newFib.size();
//            //cout << i << endl;
//        }
//        //m_fibers[i].setValid(m_fibers[i].cut(m_surface));
//        //if (!m_fibers[i].is_valid()) toDelete++;
//    }
//    //Remove fibers not used
////    unsigned int position = 0;
////    for (unsigned int i=0; i<size(); i++)
////    {
////        if (m_fibers[i].is_valid())
////        {
////            fibers.push_back(m_fibers[i]);
////            m_nbOfPoints+=m_fibers[i].size();
////            //m_fibers.erase(m_fibers.begin()+i);

////            //m_fibers[position] = m_fibers[i];
////            //position++;
////        }
////        //else
////        //    m_nbOfPoints+=m_fibers[i].size();
////    }
//    m_fibers.clear();
//    m_fibers.reserve(fibers.size());
//    m_fibers.assign(fibers.begin(), fibers.end());
////    cout << m_fibers.size() << endl;
////    for (unsigned int i=0; i<m_fibers.size(); i++)
////    {
////        cout << i << endl;
////        cout << m_fibers[i].size() << endl;
////        cout << m_fibers[i][0] << endl;
////        cout << endl;
////    }
//}

void Bundle::setFuzzyLimit(float limit)
{
#pragma omp parallel for
    for (unsigned int i=0; i<size(); i++)
        m_fibers[i].setFuzzyLimit(limit);
}

void Bundle::setToFiberMode(bool mode)
{
#pragma omp parallel for
    for (unsigned int i=0; i<size(); i++)
        m_fibers[i].setToFiberMode(mode);
}

void Bundle::normalizeFuzzyValues()
{
    vector<float> maxValue(omp_get_max_threads());
    vector<float> maxValueCylinder(omp_get_max_threads());
    for (int i=0; i<omp_get_max_threads(); i++)
    {
        maxValue[i]=0.0f;
        maxValueCylinder[i]=0.0f;
    }
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        if (m_fibers[i].getValue()>maxValue[omp_get_thread_num()])
            maxValue[omp_get_thread_num()] = m_fibers[i].getValue();
        if (m_fibers[i].getValueOfCylinder()>maxValueCylinder[omp_get_thread_num()])
            maxValueCylinder[omp_get_thread_num()] = m_fibers[i].getValueOfCylinder();
    }
    float maxV=maxValue[0];
    float maxCV=maxValueCylinder[0];
    for (int i=1; i<omp_get_max_threads(); i++)
    {
        if (maxV<maxValue[i]) maxV=maxValue[i];
        if (maxCV<maxValueCylinder[i]) maxCV=maxValueCylinder[i];
    }
    if (maxV<0.0001f) maxV=1;
    if (maxCV<0.0001f) maxCV=1;
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        m_fibers[i].setValue(m_fibers[i].getValue()/maxV);
        m_fibers[i].setCylinderValue(m_fibers[i].getValueOfCylinder()/maxCV);
        COLOUR c = GetColour(m_fibers[i].getValue(), 0, 1);
        m_fibers[i].setColor(c.r*255, c.g*255, c.b*255);
        c=GetColour(m_fibers[i].getValueOfCylinder(), 0, 1);
        m_fibers[i].setCylinderColor(c.r*255, c.g*255, c.b*255);
    }
}


void Bundle::removeFromMultiRes(float fuzzyValueLimit)
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        if (m_fibers[i].getValueOfCylinder()<fuzzyValueLimit)
            m_fibers[i].setValid(false);
    }
}

void Bundle::colorFibers()
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        COLOUR c = GetColour(m_fibers[i].getValue(), 0, 1.0);
        m_fibers[i].setColor(c.r*255, c.g*255, c.b*255);
        c=GetColour(m_fibers[i].getValueOfCylinder(), 0, 1.0);
        m_fibers[i].setCylinderColor(c.r*255, c.g*255, c.b*255);
    }
}

void Bundle::randomColorFibers()
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
    {
        float r1 = (float)rand()/RAND_MAX;
        float r2 = (float)rand()/RAND_MAX;
        float r3 = (float)rand()/RAND_MAX;
        m_fibers[i].setColor(r1*255, r2*255, r3*255);
        m_fibers[i].setCylinderColor(r1*255, r2*255, r3*255);
    }
}

void Bundle::findOutliers()
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
        m_fibers[i].setOutlier(true);
#pragma omp parallel for
    for (unsigned int i=0; i<m_changes.size(); i++)
    {
        m_fibers[m_changes[i].fib1].setOutlier(false);
        m_fibers[m_changes[i].fib2].setOutlier(false);
        m_fibers[m_changes[i].newFib].setOutlier(false);
    }
}

void Bundle::displayOutliers(bool d)
{
#pragma omp parallel for
    for (unsigned int i=0; i<m_fibers.size(); i++)
        if (m_fibers[i].is_outlier())
            m_fibers[i].setValid(d);
}
