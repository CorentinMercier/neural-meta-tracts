#pragma once

#ifndef SCENE_H
#define SCENE_H

#include <string>
#include <iostream>
#include <chrono>
#include "../opengl/openglutils.h"
#include "../MultiScale/include/bundle.h"
#include "../MultiScale/include/surface.h"
#include "../MultiScale/include/fuzzy.h"
#include "fiber_helper.h"
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>

#define BUFFER_OFFSET(i) ((void*)(i))

using namespace std;

class WidgetOpenGL;

class scene
{
public:
    scene(WidgetOpenGL* widget_param, string surfaceFile, string fuzzyFile, string mask1File, string mask2File, bool multipleTractograms=false);
    ~scene();

    /**  Method called only once at the beginning (load off files ...) */
    bool load_scene(string filename);

    /**  Method called at every frame */
    void draw_scene();

    /** Set the pointer to the parent Widget */
    void set_widget(WidgetOpenGL* widget_param);

    /* Record surface file name */
    void set_surface_name(string surface){m_surfaceFile = surface;}
    /* Record fuzzy set file name */
    void set_fuzzy_name(string fuzzy){m_fuzzyFile = fuzzy;}
    /* Record masks set file names */
    void set_mask_name(string mask1, string mask2){m_mask1File = mask1; m_mask2File = mask2;}

    /** Start/stop the animation */
    void start_stop();

    bool fiber_draw = true;
    bool cylinder_draw = false;
    bool wireframe_mode = false;
    bool merge_draw = false;
    bool only_merge_draw = false;
    bool outliers = true;
    bool corticalSurface_draw = false;
    string slider_set_value(int value);
    void slider_geometry_set_value(float value);
    string slider_fuzzy_change_value(float value);

    string changeMultiscale(int value);
    void drawNeighbors();
    void drawOneNeighbour();

    void changeInner(float change){m_inner+=change;}
    void changeOuter(float change){m_outer+=change;}

    void valueChanged(){m_valueChange=true;}

    void smoothGeometry();
    void smoothPoints();

    void changeColorState(bool random) {/*m_randomColor=random;*/
        if (random)
            m_brain->colorFibers();
        else
            m_brain->randomColorFibers();
        m_valueChange=true;
    }

    void outlierAction(bool showOutliers);

    /*Save current displayed fibers into a vtk file*/
    void saveVTK();

private:
    bool m_multipleTractograms=false;
    //Access to the parent object
    WidgetOpenGL* m_pwidget;

    //GL programs
    GLuint m_basicProgram;
    GLuint m_advancedProgram;
    GLuint m_advancedProgramWireframe;
    GLuint m_singleColorProgram;
    GLuint m_surfaceProgram;

    //Metric
    Metrique *m_metric;

    //Fibers
    Bundle *m_brain;
    GLuint m_eboSurface;
    GLuint m_vboSurface;

    //Surface
    Surface* m_surface;
    void loadSurfaceBuffers();
    bool m_thereIsASurface = false;
    string m_surfaceFile;

    //Fuzzy grid
    float m_fuzzySlider = 0.0f;
    string m_fuzzyFile;
    string m_mask1File;
    string m_mask2File;

    //Storage for the VBO associated to the fibers seen as lines
    GLuint m_vboFibers;
    GLuint m_vao;
    GLuint m_ebo;
    int m_numberOfFiberLines;
    int m_numberOfPoints;
    GLuint* m_vertexIndices;
    unsigned int m_nbOfIndices=0;
    unsigned int m_currentNbOfPoints=0;
    GLint* m_count;
    GLint* m_first;

    //Storage for the VBO associated to the fibers seen as cylinders (tesselated)
    GLuint m_vboCylinders;
    float m_inner=4;
    float m_outer=4;

    //Storage for neighbours
    vector<unsigned int> m_validNeighborsFibers;
    GLuint m_vboNeighbours;
    unsigned int m_nbOfNeighbours=0;
    bool m_drawNeighbours=false;
    bool m_drawAllNeighbours=false;
    bool m_changedOnce=false;

    //Color Mode
    bool m_randomColor=true;

    //Slider
    int m_value_of_slider=100;
    bool m_valueChange=false;
    float m_sliderGeometry=1.0;
};

#endif // SCENE_H
