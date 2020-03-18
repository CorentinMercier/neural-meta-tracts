#ifndef BUNDLE_H
#define BUNDLE_H

#include "limits.h"
#include <forward_list>
#include <fstream>
#include <iostream>
#include "fiber.h"
#include "metrique.h"
#include "priorityqueue.h"
#include "kdtree.h"
#include "fuzzy.h"
#include <QFileDialog>
#include <iomanip>
#include "libraries.h"

#ifdef USE_VTK
#include <vtkGenericDataObjectReader.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#endif
#include "../tetgen1.5.0/tetgen.h"

#define CHUNK 50
#define TEST_FIBERS_VALIDITY false
#define STATS false

class WidgetOpenGL;

using namespace std;

struct triplet
{
	 unsigned int fib1;
	 unsigned int fib2;
	 unsigned int newFib;
};

//Class containing all fibers and doing computations on all fibers
class Bundle
{
public:
	Bundle(string fileFuzzySet, string fileEndpoints1, string fileEndpoints2);
	Bundle();
	~Bundle();

	/*Loading of the set of fibers named filename, boolean reduces the quantity of fibers by taking only some of them*/
	void loadFibers(string filename, bool reduceNbFibers=false);
	/*Return a copy of the vector of all fibers*/
	vector<vector<glm::vec3> > getFibers();

	/*Return a reference to the fiber number element*/
	Fiber& operator[](unsigned int element);

	/*Add a fiber to the bundle*/
	void addFiber(Fiber &new_fiber);
	/*Add a fiber to the bundle and return bundle size*/
	unsigned int addFiberSize(Fiber &new_fiber);

	/*Compute distances between neighboring fibers*/
	void computeDistances();
	/*Merge fibers according to the distances previously computed*/
	void computeMultiScale();

	/*Orient fibers according to fiber number reference (not usually relevant)*/
	void reOrienteFibers(unsigned int reference=0);

	/*Return the number of fibers in the bundle*/
	unsigned int size() {return m_fibers.size();}

	/*Validate and invalidate fibers such that the current resolution changes from diffNbFibres fibers*/
	void setMultiScale(int diffNbFibres);
	/*Validate and invalidate fibers such that the current resolution changes to the specified percentage of fibers, 100% being the original full resolution*/
	void setMultiScale(float percentage);
	/*Return the number of resolutions computed*/
	unsigned int getNbChanges() {return m_changes.size();}
	/*Return the last merge that occured at the current resolution*/
	triplet getLastChange() {return m_changes[m_lastChange];}

	/*Old function to display the hierarchy of the fiber number m_positionFibreHierarchie*/
	void getOriginOfHierarchy();

	/*Change the fiber from which to display the hierarchy*/
	void fibrePrecedenteHierarchie(){if (m_positionFibreHierarchie!=0) m_positionFibreHierarchie--;}
	void fibreSuivanteHierarchie(){m_positionFibreHierarchie++; if (m_positionFibreHierarchie>m_changes.size()-1) m_positionFibreHierarchie=m_changes.size()-1;}

	/*Find clusters -- to test*/
	void colorSubBundles();

	/*Resampling of the fibers to obtain resampleValue points per fiber*/
	void resample(unsigned int resampleValue);

	/*Delaunay computation to find neighborhood*/
	void computeDelaunay();

	/*Instead of Delaunay, allows for all fibers comparisons*/
	void makeAllFibersNeighbors();

	/*Print fibers numbers and their respective numbers of neighbors in a file*/
	void printNeighboursInFile();

	/*return edges between fibers to draw the Delaunay tetrahedralization*/
	vector<glm::vec3> drawNeighbors();
	/*return edges between fibers similar to one fiber*/
	vector<glm::vec3> drawOneNeighbour();
	/*return the numbers of fibers similar to one fiber*/
	vector<unsigned int> drawNeighborsFibers();

	/*return the number of fibers in the original bundle*/
	unsigned int getOriginalSize() {return m_sizeWhenComputingMultiScale;}

	/*return the number of fibers at the current resolution*/
	unsigned int getNbOfValidFibers() {return m_sizeWhenComputingMultiScale-m_lastChange;}
	/*return the total number of points in the bundle*/
	unsigned int getNbOfPoints() {return m_nbOfPoints;}

	/*save the bundle with all its resolutions in a .fred file*/
	bool recordBundle();
	/*load a bundle with all its resolutions from a .fred file*/
	bool loadBundle(string filename);
	/*Save the current fibers in a vtk file named filename*/
	bool saveCurrentFibersAsVTK(string filename);

	/*define the widget in which the bundle is displayed*/
	void set_widget(WidgetOpenGL* widget_param);

	/*Smooth the geometric representation*/
	void geometricSmoothing();
	/*Smooth the points -- function not to use */
	void pointsSmoothing();

	/* Defines the surface object to be able to cut fibers */
	void setSurface(Surface *surface){m_surface = surface;}

	/* Cut the fibers to only keep fibers inside the gray matter (so inside the surface) */
	void cutFibers();

	/*Set the fuzzy limit for all fibers*/
	void setFuzzyLimit(float limit);

	/*Change from fiber to cylinder mode for segmentation*/
	void setToFiberMode(bool mode);

	/*Normalize fuzzy values*/
	void normalizeFuzzyValues();

	/*Compute the fuzzy value of one fiber*/
	void computeFiberFuzzy(unsigned int i);

	/*Compute the fuzzy values for all fibers*/
	void computeFuzzy();

	/*Remove all fibers with a fuzzy value inferior to the threshold from the multi-resolution*/
	void removeFromMultiRes(float fuzzyValueLimit);

	/*Color fibers according to their fuzzy values*/
	void colorFibers();
	void randomColorFibers();

	/*Find outliers, marking fibers as being outliers or not*/
	void findOutliers();
	/*Change the m_valid parameters of outlier fibers*/
	void displayOutliers(bool d);

	unsigned int getCurrentNberOfPoints();
	unsigned int getCurrentNberOfFibers();

	/*Function to check if there is a fuzzy set*/
	bool thereIsAFuzzySet() {return (m_fuzzySet);}

	void setMetric(Metrique *metric) {m_metric = metric;}

private:
	void loadTckFibers(const string & filename, bool reduceNbFibers);
	void loadVtkFibers(const string & filename, bool reduceNbFibers);

	//Metric
	Metrique *m_metric;

	//float m_limitAngle=cos(90*M_PI/180);
	float m_limitAngle=0;
	//Access to the parent object
	WidgetOpenGL* m_pwidget;

	vector<Fiber> m_fibers;
	kdtree* m_tree;
	priority_queue<paire, vector<paire>, comparisonPaires > m_paires;
	vector<paire> m_pairesVector;
	vector<triplet> m_changes;
	float m_percentage=0;
	int m_fin;
	int m_lastChange;
	float m_maxLengthDifBetweenFibers=std::numeric_limits<float>::max();
	float m_maxDistBetweenEndPoints=std::numeric_limits<float>::max();
	unsigned int m_sizeWhenComputingMultiScale=0;
	unsigned int m_nbOfPoints=0;

	void changeDisplayMultiScale(int debut, int fin);
	void updateNeighbours();

	unsigned int m_positionFibreHierarchie;
	unsigned int m_depth=0;

	void computeFiberFuzzyNotParallel(unsigned int i);

	//Surface of the WM/GM interface
	Surface *m_surface;

	Fuzzy *m_fuzzySet = nullptr;
	Fuzzy *m_maskEndpoints1 = nullptr;
	Fuzzy *m_maskEndpoints2 = nullptr;
	float m_gaussianDistance = 32;
};

#endif // BUNDLE_H
