#include "scene.h"
#include "../interface/widgetopengl.h"


//typedef struct {
//    double r,g,b;
//} COLOUR;

//COLOUR GetColour(double v,double vmin,double vmax)
//{
//   COLOUR c = {1.0,1.0,1.0}; // white
//   double dv;

//   if (v < vmin)
//      v = vmin;
//   if (v > vmax)
//      v = vmax;
//   dv = vmax - vmin;

//   if (v < (vmin + 0.25 * dv)) {
//      c.r = 0;
//      c.g = 4 * (v - vmin) / dv;
//   } else if (v < (vmin + 0.5 * dv)) {
//      c.r = 0;
//      c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
//   } else if (v < (vmin + 0.75 * dv)) {
//      c.r = 4 * (v - vmin - 0.5 * dv) / dv;
//      c.b = 0;
//   } else {
//      c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
//      c.b = 0;
//   }

//   return(c);
//}

scene::scene(WidgetOpenGL *widget_param, string surfaceFile, string fuzzyFile, string mask1File, string mask2File, bool multipleTractograms):
    m_multipleTractograms(multipleTractograms), m_pwidget(widget_param), m_surfaceFile(surfaceFile)
{
    m_brain = new Bundle(fuzzyFile, mask1File, mask2File);
    if (m_brain->thereIsAFuzzySet())
        m_metric = new Metrique(WEIGHTEDCURRENTSFUZZY);
    else
        m_metric = new Metrique(WEIGHTEDCURRENTSSIMILARITY);
    m_brain->setMetric(m_metric);
}

scene::~scene()
{
    delete(m_brain);
    if (!m_multipleTractograms)
        delete(m_surface);
}

bool scene::load_scene(string filename)
{
    //*****************************************//
    // Preload default structure               //
    //*****************************************//
    if (!m_multipleTractograms)
    {
        ShaderInfo  shaders[] =
        {
            { GL_VERTEX_SHADER, "../sources/opengl/shaders/shader.vert", 0 },
            { GL_FRAGMENT_SHADER, "../sources/opengl/shaders/shader.frag", 0 },
            { GL_NONE, nullptr, 0 }
        };
        m_basicProgram=loadShaders(shaders);
        if (m_basicProgram==0)
            cerr << "Error due to basic shaders" << endl;

        ShaderInfo advancedShaders[] =
        {
            { GL_VERTEX_SHADER, "../sources/opengl/shaders/advancedShader.vert", 0 },
            { GL_TESS_CONTROL_SHADER, "../sources/opengl/shaders/advancedShader.cont", 0},
            { GL_TESS_EVALUATION_SHADER, "../sources/opengl/shaders/advancedShader.eval", 0},
            { GL_FRAGMENT_SHADER, "../sources/opengl/shaders/advancedShader.frag", 0 },
            { GL_NONE, nullptr, 0 }
        };
        m_advancedProgram=loadShaders(advancedShaders);
        if (m_advancedProgram==0)
            cerr << "Error due to advanced shaders" << endl;

        ShaderInfo advancedShadersWireframe[] =
        {
            { GL_VERTEX_SHADER, "../sources/opengl/shaders/advancedShader.vert", 0 },
            { GL_TESS_CONTROL_SHADER, "../sources/opengl/shaders/advancedShader.cont", 0},
            { GL_TESS_EVALUATION_SHADER, "../sources/opengl/shaders/advancedShader.eval", 0},
            { GL_FRAGMENT_SHADER, "../sources/opengl/shaders/advancedShaderWireframe.frag", 0 },
            { GL_NONE, nullptr, 0 }
        };
        m_advancedProgramWireframe=loadShaders(advancedShadersWireframe);
        if (m_advancedProgramWireframe==0)
            cerr << "Error due to advanced shaders wireframe" << endl;

        ShaderInfo singleColorShaders[] =
        {
            { GL_VERTEX_SHADER, "../sources/opengl/shaders/shaderSingleColor.vert", 0 },
            { GL_FRAGMENT_SHADER, "../sources/opengl/shaders/shaderSingleColor.frag", 0 },
            { GL_NONE, nullptr, 0 }
        };
        m_singleColorProgram=loadShaders(singleColorShaders);
        if (m_singleColorProgram==0)
            cerr << "Error due to single color shaders" << endl;

        ShaderInfo surfaceShaders[] =
        {
            { GL_VERTEX_SHADER, "../sources/opengl/shaders/shaderSurface.vert", 0 },
            { GL_FRAGMENT_SHADER, "../sources/opengl/shaders/shaderSurface.frag", 0 },
            { GL_NONE, nullptr, 0 }
        };
        m_surfaceProgram=loadShaders(surfaceShaders);
        if (m_surfaceProgram==0)
            cerr << "Error due to surface shaders" << endl;
    }

    //*****************************************//
    // Load fibers                             //
    //*****************************************//

    auto start = chrono::steady_clock::now();

    if (!m_multipleTractograms && !m_surfaceFile.empty())
    {
        cout << " Load surface from file " << m_surfaceFile << " ... " << endl;
        m_surface = new Surface(m_surfaceFile);
        loadSurfaceBuffers();
        cout << "   [OK] " << endl;
    }

    // Load fiber data
    std::size_t found = filename.find(".fred");
    if (found==std::string::npos)//Not a fred file, so it is a vtk file
    {
        cout << "  Load data from file " << filename << " ... " << endl;
        m_brain->loadFibers(filename, false);
        cout << "   [OK] " << endl;

        if (!m_multipleTractograms)
        {
            //Cutting of the fibers
            //cout << " Cutting of the fibers ..." << endl;
            cout << " Loading cortical surface ..." << endl;
            m_brain->setSurface(m_surface);
            //m_brain->cutFibers();
            cout << "   [OK] " << endl;
        }

        cout << m_brain->size() << " " << m_brain->getNbOfPoints() << endl;

        if (m_brain->thereIsAFuzzySet())
        {
            cout << "  Fuzzy values computation ..." << endl;
            //m_brain->resample(100);
            m_brain->computeFuzzy();
            //m_brain->removeFromMultiRes(0.1f);
            cout << "   [OK] " << endl;
        }

        cout << "  Delaunay Tetrahedralization ... " << endl;
        m_brain->computeDelaunay();
        //m_brain->makeAllFibersNeighbors();
        cout << "   [OK] " << endl;

        cout << "  Distances computation ..." << endl;
        m_brain->computeDistances();
        cout << "   [OK] " << endl;

        cout << "  Creation of the multi-resolution model ... " << endl;
        m_brain->computeMultiScale();
        cout << "   [OK] " << endl;
    }
    else //fred file
    {
        cout << " Loading of the file " << filename << " ..." << endl;
        m_brain->loadBundle(filename);
        cout << "   [OK] " << endl;

        if (m_brain->thereIsAFuzzySet())
        {
            cout << "  Fuzzy values computation ..." << endl;
            m_brain->computeFuzzy();
            cout << "   [OK] " << endl;
        }
    }

    if (m_brain->thereIsAFuzzySet())
    {
        m_brain->setFuzzyLimit(0.0f);
        cout << "Normalization of the values" << endl;
        m_brain->normalizeFuzzyValues();
        cout << "   [OK] " << endl;
        m_brain->colorFibers();
    }

    auto end = chrono::steady_clock::now();
    auto diff = end-start;
    cout << "Multi-resolution computational time : " << chrono::duration <double, milli> (diff).count() << " ms" << endl;
    //return false;

    //vector<vector<Vector3f> > fibers = m_brain->getFibers();
    //m_brain->printNeighboursInFile();

    cout << "   [OK] " << endl;

    if (m_multipleTractograms)
        return true;

    //*****************************************//
    // Create OPENGL VBO for the fibers        //
    //*****************************************//
    {
        cout << "  Fill all fibers as VBO ... " << endl;

        // We start with the minimum displayed to not fill the memory
        m_brain->setMultiScale(100.f);
        // Store all vertices consecutively for OpenGL drawing as GL_LINES
        vector<Vector3f> fibers;
        vector<Vector3f> colors;
        m_numberOfFiberLines=0;
        m_first=new GLint[m_brain->getNbOfValidFibers()];
        m_count=new GLint[m_brain->getNbOfValidFibers()];
        m_vertexIndices=new GLuint[3];
        int counter=0;
        for (unsigned int i=0; i<m_brain->size(); i++)
        {
            if ((*m_brain)[i].is_valid())
            {
                for (unsigned int n=0; n<(*m_brain)[i].size(); n++)
                {
                    fibers.push_back((*m_brain)[i][n]);
                    colors.push_back((*m_brain)[i].getColor());
                }
                m_first[m_numberOfFiberLines]=counter;
                m_count[m_numberOfFiberLines]=static_cast<int>((*m_brain)[i].size());
                counter+=(*m_brain)[i].size();
                m_numberOfFiberLines++;
            }
        }

        glUseProgram(m_basicProgram);

        glGenVertexArrays(1, &m_vao);
        glBindVertexArray(m_vao);

        glCreateBuffers(1, &m_vboNeighbours);
        glCreateBuffers(1, &m_vboCylinders);
        // Set up the element array buffer
        glCreateBuffers(1, &m_ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo);
        // Create the VBO and fill the data
        glCreateBuffers(1, &m_vboFibers); PRINT_OPENGL_ERROR();
        glBindBuffer(GL_ARRAY_BUFFER, m_vboFibers);                 PRINT_OPENGL_ERROR();
        glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*fibers.size()+3*sizeof(float)*colors.size(), nullptr, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
        glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*fibers.size(), fibers.data());PRINT_OPENGL_ERROR();
        glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*fibers.size(), 3*sizeof(float)*colors.size(), colors.data());PRINT_OPENGL_ERROR();

        //glNamedBufferStorage(m_vboFibers, 3*sizeof(float)*lines.size(), lines.data(), GL_DYNAMIC_STORAGE_BIT); PRINT_OPENGL_ERROR();
        //        glNamedBufferStorage(m_vboFibers, sizeof(vertex_positions), vertex_positions, GL_DYNAMIC_STORAGE_BIT); PRINT_OPENGL_ERROR();
        //Autre possibilité :
        //        glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*lines.size(), NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
        //        void * data = glMapBuffer(GL_COPY_WRITE_BUFFER, GL_WRITE_ONLY);
        //        fread(data, 1, 3*sizeof(float)*lines.size(), lines.data());
        //        glUnmapBuffer(GL_COPY_WRITE_BUFFER);

        GLint loc1=glGetAttribLocation(m_basicProgram, "position");
        glVertexAttribPointer(loc1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);PRINT_OPENGL_ERROR();
        GLint loc2=glGetAttribLocation(m_basicProgram, "color");
        glVertexAttribPointer(loc2, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*fibers.size()));PRINT_OPENGL_ERROR();

        glEnableVertexAttribArray(loc1);PRINT_OPENGL_ERROR();
        glEnableVertexAttribArray(loc2);PRINT_OPENGL_ERROR();
        m_numberOfPoints = fibers.size();

        cout << "   [OK] " << endl;
    }

    m_brain->findOutliers();
    cout << "Run drawing" << endl;

    end = chrono::steady_clock::now();
    diff = end-start;
    cout << "Total computational time : " << chrono::duration <double, milli> (diff).count() << " ms" << endl;
    return true;
}


void scene::draw_scene()
{
    //    if (m_valueChange)
    //    {
    //        m_brain->setMultiScale((float)m_value_of_slider);
    //        m_valueChange=false;
    //    }
    //    struct stat result;
    //    if (stat("../opengl/shaders/advancedShader.eval", &result)==0)
    //    {
    //        auto mod_time = result.st_mtime;
    //    }

    bool valueChanged=m_valueChange;

    if (m_thereIsASurface && corticalSurface_draw)
    {
        glUseProgram(m_surfaceProgram);
        glUniformMatrix4fv(get_uni_loc(m_surfaceProgram,"camera_modelview"),1,false, m_pwidget->cam.getModelview().data());  PRINT_OPENGL_ERROR();
        glUniformMatrix4fv(get_uni_loc(m_surfaceProgram,"camera_projection"),1,false,m_pwidget->cam.getProjection().data());   PRINT_OPENGL_ERROR();
        Vector3f gray(0.5, 0.5, 0.5);
        glUniform3fv(get_uni_loc(m_surfaceProgram,"singleColor"),1,gray.data());   PRINT_OPENGL_ERROR();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_eboSurface); PRINT_OPENGL_ERROR();
        glBindBuffer(GL_ARRAY_BUFFER, m_vboSurface); PRINT_OPENGL_ERROR();
        GLint loc1=glGetAttribLocation(m_surfaceProgram, "position"); PRINT_OPENGL_ERROR();
        glVertexAttribPointer(loc1, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
        GLint loc2=glGetAttribLocation(m_surfaceProgram, "normal");
        glVertexAttribPointer(loc2, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_surface->getNbPoints()));PRINT_OPENGL_ERROR();
        glDrawElements(GL_TRIANGLES, m_surface->getNbTriangles()*3, GL_UNSIGNED_INT, nullptr);
    }

    if (fiber_draw)
    {
        //glEnable(GL_LINE_SMOOTH);
        glEnable(GL_MULTISAMPLE);
        glUseProgram(m_basicProgram);
        glUniformMatrix4fv(get_uni_loc(m_basicProgram,"camera_modelview"),1,false, m_pwidget->cam.getModelview().data());  PRINT_OPENGL_ERROR();
        glUniformMatrix4fv(get_uni_loc(m_basicProgram,"camera_projection"),1,false,m_pwidget->cam.getProjection().data());   PRINT_OPENGL_ERROR();

        if(valueChanged)
        {
            m_brain->setToFiberMode(true);
            // Store all vertices consecutively for OpenGL drawing as GL_LINES
            vector<Vector3f> fibers;
            vector<Vector3f> colors;
            m_numberOfFiberLines=0;
            delete(m_first);
            delete(m_count);
            //            m_first=new GLint[m_brain->getNbOfValidFibers()];
            //            m_count=new GLint[m_brain->getNbOfValidFibers()];
            m_first=new GLint[m_brain->getCurrentNberOfFibers()];
            m_count=new GLint[m_brain->getCurrentNberOfFibers()];
            unsigned int counter=0;
            for (unsigned int i=0; i<m_brain->size(); i++)
            {
                if ((*m_brain)[i].is_valid())
                {
                    for (unsigned int n=0; n<(*m_brain)[i].size(); n++)
                    {
                        fibers.push_back((*m_brain)[i][n]);
                        colors.push_back((*m_brain)[i].getColor());
                        //Colors depending on inside/outside
                        //                        if (m_surface->isInside((*m_brain)[i][n]))
                        //                            colors.push_back(Vector3f(0,0,1));
                        //                        else
                        //                        {
                        //                            colors.push_back(Vector3f(1,0,0));
                        //                        }
                        //Colors depending on Fuzzy set
                        //                        int value = m_fuzzy->getValueFromPoint((*m_brain)[i][n]);
                        //                        COLOUR c = GetColour(value, 0,255);
                        //                        colors.push_back(Vector3f(c.r, c.g, c.b));
                    }
                    m_first[m_numberOfFiberLines]=counter;
                    m_count[m_numberOfFiberLines]=(*m_brain)[i].size();
                    counter+=(*m_brain)[i].size();
                    m_numberOfFiberLines++;
                }
            }

            // Create the VBO and fill the data
            glBindBuffer(GL_ARRAY_BUFFER, m_vboFibers);                 PRINT_OPENGL_ERROR();

            glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*fibers.size()+3*sizeof(float)*colors.size(), NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
            glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*fibers.size(), fibers.data());PRINT_OPENGL_ERROR();
            glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*fibers.size(), 3*sizeof(float)*colors.size(), colors.data());PRINT_OPENGL_ERROR();

            GLint loc1=glGetAttribLocation(m_basicProgram, "position");
            GLint loc2=glGetAttribLocation(m_basicProgram, "color");

            glEnableVertexAttribArray(loc1);PRINT_OPENGL_ERROR();
            glEnableVertexAttribArray(loc2);PRINT_OPENGL_ERROR();

            m_numberOfPoints=fibers.size();
            m_valueChange=false;
        }
        glBindBuffer(GL_ARRAY_BUFFER, m_vboFibers);                 PRINT_OPENGL_ERROR();
        GLint loc1=glGetAttribLocation(m_basicProgram, "position");
        glVertexAttribPointer(loc1, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
        GLint loc2=glGetAttribLocation(m_basicProgram, "color");
        glVertexAttribPointer(loc2, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_numberOfPoints));PRINT_OPENGL_ERROR();

        glEnableVertexAttribArray(loc1);PRINT_OPENGL_ERROR();
        glEnableVertexAttribArray(loc2);PRINT_OPENGL_ERROR();

        glMultiDrawArrays(GL_LINE_STRIP, m_first, m_count, m_numberOfFiberLines);
    }
    if (cylinder_draw || wireframe_mode || (!outliers && cylinder_draw))
    {
        glUseProgram(m_advancedProgram);
        glUniformMatrix4fv(get_uni_loc(m_advancedProgram,"camera_modelview"),1,false, m_pwidget->cam.getModelview().data());  PRINT_OPENGL_ERROR();
        glUniformMatrix4fv(get_uni_loc(m_advancedProgram,"camera_projection"),1,false,m_pwidget->cam.getProjection().data());   PRINT_OPENGL_ERROR();
        glUniform1f(get_uni_loc(m_advancedProgram, "Outer"), m_outer);
        glUniform1i(get_uni_loc(m_advancedProgram, "colorFromOrientation"), !m_randomColor);
        glUniform1i(get_uni_loc(m_advancedProgram, "neighbors"), m_drawNeighbours);
        glUniform1f(get_uni_loc(m_advancedProgram, "geometricSize"), m_sliderGeometry);
        if (valueChanged)
        {
            m_brain->setToFiberMode(false);
            if (!outliers)//We want to remove the outliers
            {
                m_currentNbOfPoints=0;
                m_nbOfIndices=0;
                for (unsigned int i=0; i<m_brain->size(); i++)
                {
                    if ((*m_brain)[i].is_valid() && ((*m_brain)[i].getNbFibHierarchy()>1))
                    {
                        m_currentNbOfPoints+=(*m_brain)[i].size();
                        m_nbOfIndices+=(*m_brain)[i].size()-1;
                    }
                }
            }
            else
            {
                m_currentNbOfPoints=m_brain->getCurrentNberOfPoints();
                m_nbOfIndices=m_currentNbOfPoints-m_brain->getCurrentNberOfFibers();
            }

            // Store all vertices consecutively for OpenGL
            vector<Vector3f> fibers(m_currentNbOfPoints);
            vector<Vector3f> colors(m_currentNbOfPoints);

            //Ellipse parameters
            vector<float> ab(2*m_currentNbOfPoints); //a then b, so to read as a vec2 in GLSL
            vector<Vector3f> ellipseOrientation(m_currentNbOfPoints); //Vector of direction from the axes a
            delete(m_vertexIndices);
            m_numberOfFiberLines=0;
            m_vertexIndices=new GLuint[m_nbOfIndices*4];
            unsigned int counter=0;
            unsigned int position=0;
            unsigned int counterNeighbors=0;
            Vector3f red(-1.0, 0.0, 0.0);
            bool drawNeighbor=m_drawNeighbours;
            if (drawNeighbor)
            {
                sort(m_validNeighborsFibers.begin(), m_validNeighborsFibers.end());
            }
            for (unsigned int i=0; i<m_brain->size(); i++)
            {
                if ((*m_brain)[i].is_valid() && (outliers || (*m_brain)[i].getNbFibHierarchy()>1))
                {
                    for (unsigned int n=0; n<(*m_brain)[i].size(); n++)
                    {
                        if (n<(*m_brain)[i].size()-1)
                        {
                            if (n==0)
                                m_vertexIndices[position]=counter+n;
                            else
                                m_vertexIndices[position]=counter+n-1;
                            m_vertexIndices[position+1]=counter+n;
                            m_vertexIndices[position+2]=counter+n+1;
                            if (n<(*m_brain)[i].size()-2)
                                m_vertexIndices[position+3]=counter+n+2;
                            else
                                m_vertexIndices[position+3]=counter+n+1;
                            position+=4;
                        }
                        fibers[counter+n]=(*m_brain)[i].getProfileTransform(n).center;//(*m_brain)[i][n];
                        if (drawNeighbor && m_validNeighborsFibers[counterNeighbors]==i)
                        {
                            colors[counter+n]=red;
                        }
                        else
                            colors[counter+n]=(*m_brain)[i].getCylinderColor();
                        ab[2*counter+2*n]=(*m_brain)[i].getProfileTransform(n).a;
                        ab[2*counter+2*n+1]=(*m_brain)[i].getProfileTransform(n).b;
                        ellipseOrientation[counter+n]=(*m_brain)[i].getProfileTransform(n).ellipseOrientation;
                    }
                    if (drawNeighbor && m_validNeighborsFibers[counterNeighbors]==i)
                    {
                        counterNeighbors++;
                        if (counterNeighbors>m_validNeighborsFibers.size()-1)
                            drawNeighbor=false;
                    }
                    counter+=(*m_brain)[i].size();
                    m_numberOfFiberLines++;
                }
                else if (drawNeighbor && m_validNeighborsFibers[counterNeighbors]==i)
                {
                    counterNeighbors++;
                    if (counterNeighbors>m_validNeighborsFibers.size()-1)
                        drawNeighbor=false;
                }
            }

            // Set up the element array buffer
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo); PRINT_OPENGL_ERROR();
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_nbOfIndices*4*sizeof(GLuint), m_vertexIndices, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();

            // Create the VBO and fill the data

            glBindBuffer(GL_ARRAY_BUFFER, m_vboCylinders);                 PRINT_OPENGL_ERROR();

            glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*fibers.size()+3*sizeof(float)*colors.size()+sizeof(float)*ab.size()+3*sizeof(float)*ellipseOrientation.size(), NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
            glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(float)*fibers.size(), fibers.data());PRINT_OPENGL_ERROR();
            glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*fibers.size(), 3*sizeof(float)*colors.size(), colors.data());PRINT_OPENGL_ERROR();
            glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*fibers.size()+3*sizeof(float)*colors.size(), sizeof(float)*ab.size(), ab.data());PRINT_OPENGL_ERROR();
            glBufferSubData(GL_ARRAY_BUFFER, 3*sizeof(float)*fibers.size()+3*sizeof(float)*colors.size()+sizeof(float)*ab.size(), 3*sizeof(float)*ellipseOrientation.size(), ellipseOrientation.data());PRINT_OPENGL_ERROR();

            GLint loc1=glGetAttribLocation(m_advancedProgram, "position");
            GLint loc2=glGetAttribLocation(m_advancedProgram, "color");
            GLint loc3=glGetAttribLocation(m_advancedProgram, "ab");
            GLint loc4=glGetAttribLocation(m_advancedProgram, "ellipseOrientation");

            glEnableVertexAttribArray(loc1);PRINT_OPENGL_ERROR();
            glEnableVertexAttribArray(loc2);PRINT_OPENGL_ERROR();
            glEnableVertexAttribArray(loc3);PRINT_OPENGL_ERROR();
            glEnableVertexAttribArray(loc4);PRINT_OPENGL_ERROR();

            m_numberOfPoints=fibers.size();
            m_valueChange=false;
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo); PRINT_OPENGL_ERROR();
        glBindBuffer(GL_ARRAY_BUFFER, m_vboCylinders);                 PRINT_OPENGL_ERROR();
        GLint loc1=glGetAttribLocation(m_advancedProgram, "position");
        glVertexAttribPointer(loc1, 3, GL_FLOAT, GL_FALSE, 0, NULL); PRINT_OPENGL_ERROR();
        glPatchParameteri(GL_PATCH_VERTICES, 4); PRINT_OPENGL_ERROR();
        GLint loc2=glGetAttribLocation(m_advancedProgram, "color");
        glVertexAttribPointer(loc2, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_numberOfPoints));PRINT_OPENGL_ERROR();
        GLint loc3=glGetAttribLocation(m_advancedProgram, "ab");
        glVertexAttribPointer(loc3, 2, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_numberOfPoints+3*sizeof(float)*m_numberOfPoints));PRINT_OPENGL_ERROR();
        GLint loc4=glGetAttribLocation(m_advancedProgram, "ellipseOrientation");
        glVertexAttribPointer(loc4, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)(3*sizeof(float)*m_numberOfPoints+3*sizeof(float)*m_numberOfPoints+2*sizeof(float)*m_numberOfPoints));PRINT_OPENGL_ERROR();

        if (cylinder_draw)
        {
            glDrawElements(GL_PATCHES, m_nbOfIndices*4, GL_UNSIGNED_INT, 0);
        }
        if (wireframe_mode)
        {
            glUseProgram(m_advancedProgramWireframe);
            glUniformMatrix4fv(get_uni_loc(m_advancedProgramWireframe,"camera_modelview"),1,false, m_pwidget->cam.getModelview().data());  PRINT_OPENGL_ERROR();
            glUniformMatrix4fv(get_uni_loc(m_advancedProgramWireframe,"camera_projection"),1,false,m_pwidget->cam.getProjection().data());   PRINT_OPENGL_ERROR();
            glUniform1f(get_uni_loc(m_advancedProgramWireframe, "Outer"), m_outer);
            glUniform1f(get_uni_loc(m_advancedProgramWireframe, "geometricSize"), m_sliderGeometry);
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1,1);
            glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);
            glDrawElements(GL_PATCHES, m_nbOfIndices*4, GL_UNSIGNED_INT, 0);
            glPolygonMode( GL_FRONT_AND_BACK, GL_FILL);
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
    //Merge mode
    if (merge_draw)
    {
        //triplet changes=m_brain->getLastChange();
    }
    if (m_drawAllNeighbours)
    {
        glEnable(GL_MULTISAMPLE);
        glUseProgram(m_singleColorProgram);
        glUniformMatrix4fv(get_uni_loc(m_singleColorProgram,"camera_modelview"),1,false, m_pwidget->cam.getModelview().data());  PRINT_OPENGL_ERROR();
        glUniformMatrix4fv(get_uni_loc(m_singleColorProgram,"camera_projection"),1,false,m_pwidget->cam.getProjection().data());   PRINT_OPENGL_ERROR();
        Vector3f singleColor=Vector3f(1.0, 0.0, 0.0);
        glUniform3fv(get_uni_loc(m_singleColorProgram,"singleColor"),1,singleColor.data());   PRINT_OPENGL_ERROR();
        glBindBuffer(GL_ARRAY_BUFFER, m_vboNeighbours);                 PRINT_OPENGL_ERROR();
        GLint loc1=glGetAttribLocation(m_singleColorProgram, "position");
        glVertexAttribPointer(loc1, 3, GL_FLOAT, GL_FALSE, 0, NULL);PRINT_OPENGL_ERROR();
        glDrawArrays(GL_LINES, 0, m_nbOfNeighbours);
    }
}

void scene::loadSurfaceBuffers()
{
    //    vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    //    reader->SetFileName(m_surface.c_str());
    //    reader->Update();
    //    vtkSmartPointer<vtkPolyData> mesh  = vtkSmartPointer<vtkPolyData>::New();
    //    mesh=reader->GetPolyDataOutput();
    //    unsigned int nbVertices = mesh->GetNumberOfPoints();
    //    unsigned int nbTriangles = mesh->GetNumberOfPolys();

    //    cout << "This surface is composed of " << nbTriangles << " triangles and " << nbVertices << " vertices" << endl;

    //    //Acquiring the points of the surface
    //    m_surfacePoints.resize(nbVertices);
    //    m_surfaceNormals.resize(nbVertices);
    //    vector<int> nbNormals(nbVertices);
    //    for (unsigned int i=0; i<nbVertices; i++)
    //    {
    //        double p[3];
    //        mesh->GetPoint(i,p);
    //        m_surfacePoints[i] = Vector3f(p[0], p[1], p[2]);//p[1]-19.5, p[2]+36);
    //        m_surfaceNormals[i] = Vector3f(0);
    //        nbNormals[i]=0;
    //    }
    //    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    //    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    //    triangles = mesh->GetPolys();
    //    m_surfaceTriangles.resize(nbTriangles);
    //    for (unsigned int i=0; i<nbTriangles; i++)
    //    {
    //        triangles->GetNextCell(idList);
    //        m_surfaceTriangles[i] = Vector3i(idList->GetId(0), idList->GetId(1), idList->GetId(2));
    //    }
    //    //Creation of the normal per vertex
    //    for (unsigned int i=0; i<nbTriangles; i++)
    //    {
    //        Vector3f currentNormal = (m_surfacePoints[m_surfaceTriangles[i](1)]-m_surfacePoints[m_surfaceTriangles[i](0)]).cross(m_surfacePoints[m_surfaceTriangles[i](2)]-m_surfacePoints[m_surfaceTriangles[i](1)]);
    //        for (unsigned int j=0; j<3; j++)
    //        {
    //            nbNormals[m_surfaceTriangles[i](j)]++;
    //            m_surfaceNormals[m_surfaceTriangles[i](j)]+=currentNormal;
    //        }
    //    }
    //    for (unsigned int i=0; i<nbVertices; i++)
    //    {
    //        m_surfaceNormals[i]/=nbNormals[i];
    //        m_surfaceNormals[i].normalized();
    //    }

    //    glCreateBuffers(1, &m_vboSurface);
    //    glCreateBuffers(1, &m_eboSurface);
    //    glBindBuffer(GL_ARRAY_BUFFER, m_vboSurface); PRINT_OPENGL_ERROR();
    //    glBufferData(GL_ARRAY_BUFFER, m_surfacePoints.size()*3*sizeof(float)*2, NULL, GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
    //    glBufferSubData(GL_ARRAY_BUFFER, 0, m_surfacePoints.size()*3*sizeof(float), m_surfacePoints.data()); PRINT_OPENGL_ERROR();
    //    glBufferSubData(GL_ARRAY_BUFFER, m_surfacePoints.size()*3*sizeof(float), m_surfacePoints.size()*3*sizeof(float), m_surfaceNormals.data()); PRINT_OPENGL_ERROR();
    //    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_eboSurface); PRINT_OPENGL_ERROR();
    //    glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_surfaceTriangles.size()*3*sizeof(int), m_surfaceTriangles.data(), GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
    m_vboSurface=m_surface->getVboSurface();
    m_eboSurface=m_surface->getEboSurface();
    m_thereIsASurface = true;
}

void scene::set_widget(WidgetOpenGL *widget_param)
{
    m_pwidget=widget_param;
    m_brain->set_widget(widget_param);
}

void scene::start_stop()
{

}

string scene::slider_set_value(int value)
{
    m_value_of_slider = value;
    m_brain->setMultiScale((float)m_value_of_slider);
    m_valueChange=true;
    //draw_scene();
    unsigned int nbOfFibersDrawn=0;
    for (unsigned int i=0; i<m_brain->size(); i++)
        if ((*m_brain)[i].is_valid())
            nbOfFibersDrawn++;
    string affichage;
    affichage = "Resolution : " + to_string(100-value) + (string)"%\n"
            + to_string(nbOfFibersDrawn) + " / " + to_string(m_brain->getOriginalSize()) + " fibers\n" +
            "Last distance : \n" + to_string(m_metric->metriqueTest((*m_brain)[m_brain->getLastChange().fib1], (*m_brain)[m_brain->getLastChange().fib2])) + "\n" +
            "Segmentation : " + to_string((int)(m_fuzzySlider*255));
    //    cout << "Test" << endl;
    return affichage;
}

void scene::slider_geometry_set_value(float value)
{
    if (value>1)
        value = 1;
    if (value<0)
        value = 0;
    m_sliderGeometry = value;
}

string scene::slider_fuzzy_change_value(float value)
{
    if (value>1)
        value = 1;
    if (value<0)
        value = 0;
    m_fuzzySlider = value;
    m_brain->setFuzzyLimit(m_fuzzySlider);
    m_valueChange = true;
    unsigned int nbOfFibersDrawn=0;
    for (unsigned int i=0; i<m_brain->size(); i++)
        if ((*m_brain)[i].is_valid())
            nbOfFibersDrawn++;
    string affichage;
    affichage = "Resolution : " + to_string(100-m_value_of_slider) + (string)"%\n"
            + to_string(nbOfFibersDrawn) + " / " + to_string(m_brain->getOriginalSize()) + " fibers\n" +
//            "Last distance : \n" + to_string(metriqueTest((*m_brain)[m_brain->getLastChange().fib1], (*m_brain)[m_brain->getLastChange().fib2])) + "\n" +
            "Segmentation : " + to_string((m_fuzzySlider));
    return affichage;
}

string scene::changeMultiscale(int value)
{
    m_brain->setMultiScale(value);
    m_valueChange=true;
    unsigned int nbOfFibersDrawn=0;
    for (unsigned int i=0; i<m_brain->size(); i++)
        if ((*m_brain)[i].is_valid())
            nbOfFibersDrawn++;
    string affichage;
    affichage = "Resolution : " + to_string((int)((float)nbOfFibersDrawn/m_brain->getOriginalSize()*100)) + (string)"%\n"
            + to_string(nbOfFibersDrawn) + " / " + to_string(m_brain->getOriginalSize()) + " fibers\n" +
            "Last distance : \n" + to_string(m_metric->metriqueTest((*m_brain)[m_brain->getLastChange().fib1], (*m_brain)[m_brain->getLastChange().fib2]));
    return affichage;
}

void scene::drawNeighbors()
{
    if (m_drawAllNeighbours)
    {
        m_drawAllNeighbours=false;
        return;
    }
    //    m_validNeighborsFibers=m_brain->drawNeighborsFibers();
    //    //cout << m_validNeighborsFibers.size() << endl;
    //    m_drawNeighbours=true;
    //    m_valueChange=true;


    // Store all vertices consecutively for OpenGL drawing as GL_LINES
    vector<Vector3f> neighbours=m_brain->drawNeighbors();
    m_nbOfNeighbours=neighbours.size();
    // Create the VBO and fill the data
    glBindBuffer(GL_ARRAY_BUFFER, m_vboNeighbours);                 PRINT_OPENGL_ERROR();
    glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_nbOfNeighbours, &neighbours[0], GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
    GLint loc1=glGetAttribLocation(m_singleColorProgram, "position");
    glEnableVertexAttribArray(loc1);PRINT_OPENGL_ERROR();
    m_drawAllNeighbours=true;

}

void scene::drawOneNeighbour()
{
    if (m_drawNeighbours)
    {
        if (m_changedOnce)
        {
            m_drawNeighbours=false;
            m_changedOnce=false;
            return;
        }
        else
            m_changedOnce=true;
    }
    m_validNeighborsFibers=m_brain->drawNeighborsFibers();
    //cout << m_validNeighborsFibers.size() << endl;
    m_drawNeighbours=true;
    m_valueChange=true;


    //    if (m_drawNeighbours)
    //    {
    //        m_drawNeighbours=false;
    //        m_valueChange=true;
    //        return;
    //    }
    //    // Store all vertices consecutively for OpenGL drawing as GL_LINES
    //    vector<Vector3f> neighbours=m_brain->drawOneNeighbour();
    //    m_nbOfNeighbours=neighbours.size();
    //    // Create the VBO and fill the data
    //    glBindBuffer(GL_ARRAY_BUFFER, m_vboNeighbours);                 PRINT_OPENGL_ERROR();
    //    glBufferData(GL_ARRAY_BUFFER, 3*sizeof(float)*m_nbOfNeighbours, &neighbours[0], GL_STATIC_DRAW); PRINT_OPENGL_ERROR();
    //    GLint loc1=glGetAttribLocation(m_singleColorProgram, "position");
    //    glEnableVertexAttribArray(loc1);PRINT_OPENGL_ERROR();
    //    m_drawNeighbours=true;
}

void scene::smoothGeometry()
{
    m_brain->geometricSmoothing();
    m_valueChange=true;
}


void scene::smoothPoints()
{
    m_brain->pointsSmoothing();
    m_valueChange=true;
}

void scene::saveVTK()
{
    string filename = m_pwidget->askForVTKFilename();
    if (filename.empty())
    {
        cerr << "Error, filename not valid, not recording" << endl;
        return;
    }
    cout << "Recording of file " << filename << endl;
    if (m_brain->saveCurrentFibersAsVTK(filename))
        cout << "   [OK]" << endl;
    else
        cerr << "----------- Error, recording failed ---------------";
}

void scene::outlierAction(bool showOutliers)
{
    m_brain->displayOutliers(showOutliers);
}
