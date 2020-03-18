#include "widgetopengl.h"
#include "libraries.h"

WidgetOpenGL::WidgetOpenGL(const QGLFormat& format, QGLWidget *parent) :
	QGLWidget(format, parent)
{
#ifdef USE_VTK
	m_filename = QFileDialog::getOpenFileName(this, "Open a file containing fibers", "/media/Donnees/IsbiFiles", "vtk, tck or fred files (*.vtk *.tck *.fred)").toStdString();
#else
	m_filename = QFileDialog::getOpenFileName(this, "Open a file containing fibers", "/media/Donnees/IsbiFiles", "tck or fred files (*.tck *.fred)").toStdString();
#endif
	if (m_filename.empty())
	{
		cerr << "Error, a fiber tract needs to be specified" << endl;
		exit(EXIT_FAILURE);
	}
	std::size_t found = m_filename.find(".fred");
	if (found==std::string::npos)
		m_recordName = QFileDialog::getSaveFileName(this, "Recording of the computed bundle: ", "/media/Donnees/IsbiFiles", "fred files (*.fred)").toStdString();
#ifdef USE_VTK
	m_surface = QFileDialog::getOpenFileName(this, "Open a file containing the cortical surface", "/media/Donnees/Owncloud-Data/Data/CorticalSurface", "vtk (*.vtk)").toStdString();
#endif
	m_fuzzy = QFileDialog::getOpenFileName(this, "Open a file containing a fuzzy set", "/media/Donnees/IsbiFiles", "txt (*.txt)").toStdString();
	if (!m_fuzzy.empty())
	{
		m_setPoint1 = QFileDialog::getOpenFileName(this, "Open a file containing a first mask of endpoints", "/media/Donnees/IsbiFiles", "txt (*.txt)").toStdString();
		m_setPoint2 = QFileDialog::getOpenFileName(this, "Open a file containing a second mask of endpoints", "/media/Donnees/IsbiFiles", "txt (*.txt)").toStdString();
	}
	QWidget::setFocusPolicy(Qt::WheelFocus);
	m_scene = new scene(this, m_surface, m_fuzzy, m_setPoint1, m_setPoint2);
	cout << "Scene created" << endl;
	startTimer(25); //start timer every 25ms
}

//WidgetOpenGL::WidgetOpenGL(const QGLFormat& format, QGLWidget *parent) :
//    QGLWidget(format, parent)
//{
//    string recapFile = QFileDialog::getOpenFileName(this, "Select a file containing paths", "/media/Donnees/IsbiFiles", "txt files (*.txt)").toStdString();
//    fstream pathFile;
//    pathFile.open(recapFile, ios_base::in);
//    string bundle;
//    int nbOfSubjects;
//    pathFile >> nbOfSubjects;
//    for (unsigned int i=0; i<nbOfSubjects; i++)
//    {
//        pathFile >> m_filename;//Path of the vtk file
//        pathFile >> bundle;//Name of the bundle
//        pathFile >> m_fuzzy;//Path of the fuzzy set file
//        pathFile >> m_setPoint1;//Path of the first endpoint mask file
//        pathFile >> m_setPoint2;//Path of the second endpoint mask file
//        m_recordName = m_filename;
//        while(m_recordName.back()!='.')
//            m_recordName.pop_back();
//        m_recordName.pop_back();
//        m_recordName = m_recordName+bundle+"-V3.fred";
//        m_surface="";
//        cout << m_filename << " " << m_recordName << " " << m_surface << " " << m_fuzzy << " " << m_setPoint1 << " " << m_setPoint2 << endl;
//        m_scene = new scene(this, m_surface, m_fuzzy, m_setPoint1, m_setPoint2, true);
//        m_scene->set_widget(this);
//        m_scene->set_surface_name(m_surface);
//        m_scene->set_fuzzy_name(m_fuzzy);
//        m_scene->set_mask_name(m_setPoint1, m_setPoint2);
//        m_scene->load_scene(m_filename);
//        delete(m_scene);
//    }
//    cout << "Finished all subjects" << endl;
//    exit(EXIT_SUCCESS);
//}

WidgetOpenGL::~WidgetOpenGL()
{

}

void WidgetOpenGL::initializeGL()
{
	cout << "Opengl setup" << endl;
	//Init OpenGL
	setup_opengl();
	cout << "Done" << endl;

	//Init Camera
	cam.setupCamera();

	//Init Scene 3D
	m_scene->set_widget(this);
	m_scene->set_surface_name(m_surface);
	m_scene->set_fuzzy_name(m_fuzzy);
	m_scene->set_mask_name(m_setPoint1, m_setPoint2);
	if (!m_scene->load_scene(m_filename))
	{
		this->window()->close();
		exit(0);
	}
	//Activate depth buffer
	glEnable(GL_DEPTH_TEST); PRINT_OPENGL_ERROR();
}

void WidgetOpenGL::paintGL()
{
	//Compute current cameras
	cam.setupCamera();

	//clear screen
	glViewport (0, 0, cam.getScreenSizeX(), cam.getScreenSizeY()); PRINT_OPENGL_ERROR();
	glClearColor (1.0f, 1.0f, 1.0f, 1.0f);                      PRINT_OPENGL_ERROR();
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);        PRINT_OPENGL_ERROR();

	m_scene->draw_scene();
}


void WidgetOpenGL::resizeGL(int const width, int const height)
{
	cam.setScreenSize(width, height);
	glViewport(0,0, width, height); PRINT_OPENGL_ERROR();
}

void WidgetOpenGL::setup_opengl()
{
	setup_glad();
	print_current_opengl_context();
}

void WidgetOpenGL::setup_glad()
{
	if(!gladLoadGL())
	{
		std::cerr<<"Error initializing GLAD\n";
		exit(EXIT_FAILURE);
	}
	std::cout << "GLAD initialized\n";
}

void WidgetOpenGL::print_current_opengl_context() const
{
	std::cout << "OpenGl informations: VENDOR:       " << glGetString(GL_VENDOR)<<std::endl;
	std::cout << "                     RENDERDER:    " << glGetString(GL_RENDERER)<<std::endl;
	std::cout << "                     VERSION:      " << glGetString(GL_VERSION)<<std::endl;
	std::cout << "                     GLSL VERSION: " << glGetString(GL_SHADING_LANGUAGE_VERSION)<<std::endl;
	std::cout<<"Current OpenGL context: "<< context()->format().majorVersion() << "." << context()->format().minorVersion()<<std::endl;
}

scene& WidgetOpenGL::get_scene()
{
	return *m_scene;
}

void WidgetOpenGL::keyPressEvent(QKeyEvent *event)
{
	int current=event->key();
	Qt::KeyboardModifiers mod=event->modifiers();

	// We can quit the scene with 'Q'
	if( (mod&Qt::ShiftModifier)!=0 && (current==Qt::Key_Q) )
	{
		std::cout<<"\n[EXIT OK]\n\n"<<std::endl;
		this->window()->close();
	}
	//Camera current position is recorded to a file called Camera.txt
	if (current==Qt::Key_S)
	{
		fstream file;
		file.open("Camera.txt" , ios_base::out);
		glm::fquat q = cam.getQuat();
		file << q.x << " " << q.y << " " << q.z << " " << q.w << endl;
		glm::vec3 t=cam.getTranslation();
		file << t[0] << " " << t[1] << " " << t[2] << endl;
		file << cam.getDist() << endl;
		file.close();
		cout << "Position de la caméra enregistrée" << endl;
	}
	//Camera position is loaded from file Camera.txt
	if (current==Qt::Key_C)
	{
		fstream file;
		file.open("Camera.txt", ios_base::in);
		float x, y, z, w;
		file >> x; file >> y; file >> z; file >> w;
		glm::fquat q(w, x, y, z);
		cam.setQuaternion(q);
		file >> x; file >> y; file >> z;
		cam.setTranslation(glm::vec3(x, y, z));
		file >> x;
		cam.setDist(x);
		m_scene->draw_scene();
		updateGL(); PRINT_OPENGL_ERROR();
		cout << "Position de la caméra chargée" << endl;
	}
//	//Camera current position is recorded to a file called Camera.txt
//	if (current==Qt::Key_Q)
//	{
//		fstream file;
//		file.open("Camera2.txt" , ios_base::out);
//		Quaternion<float> q = cam.getQuat();
//		file << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << endl;
//		Vector3f t=cam.getTranslation();
//		file << t[0] << " " << t[1] << " " << t[2] << endl;
//		file << cam.getDist() << endl;
//		file.close();
//		cout << "Position de la caméra enregistrée" << endl;
//	}
//	//Camera position is loaded from file Camera.txt
//	if (current==Qt::Key_X)
//	{
//		fstream file;
//		file.open("Camera2.txt", ios_base::in);
//		float x, y, z, w;
//		file >> x; file >> y; file >> z; file >> w;
//		Quaternion<float> q(w, x, y, z);
//		cam.setQuaternion(q);
//		file >> x; file >> y; file >> z;
//		cam.setTranslation(Vector3f(x, y, z));
//		file >> x;
//		cam.setDist(x);
//		m_scene->draw_scene();
//		updateGL(); PRINT_OPENGL_ERROR();
//		cout << "Position de la caméra chargée" << endl;
//	}
	if (current==Qt::Key_X)
	{
		cam.alignX();
		updateGL(); PRINT_OPENGL_ERROR();
	}
	if (current==Qt::Key_Y)
	{
		cam.alignY();
		updateGL(); PRINT_OPENGL_ERROR();
	}
	if (current==Qt::Key_Z)
	{
		cam.alignZ();
		updateGL(); PRINT_OPENGL_ERROR();
	}
	//Draw the Delaunay Tetrahedralization, creating lines between fibers extremities (one extremity each time)
	if (current==Qt::Key_V)
	{
		m_scene->drawNeighbors();
	}
	//Draw fibers in black and white, letting one set of neighbors in red
	if (current==Qt::Key_B)
	{
		m_scene->drawOneNeighbour();
	}
	if (current==Qt::Key_P)//Smoothing of the points
	{
		m_scene->smoothPoints();
	}
	if (current==Qt::Key_M)//Smoothing of the geometry
	{
		m_scene->smoothGeometry();
	}
	//Change the tessellation level, increasing it
	if (current==Qt::Key_O)
	{
		m_scene->changeOuter(0.1);
	}
	//Change the tessellation level, decreasing it
	if (current==Qt::Key_L)
	{
		m_scene->changeOuter(-0.1);
	}
	if (current==Qt::Key_A)
	{
		m_turnCam = !m_turnCam;
	}

	QGLWidget::keyPressEvent(event);
	updateGL();
}

void WidgetOpenGL::timerEvent(QTimerEvent *event)
{
	event->accept();
	if (m_turnCam)
		cam.rotateAlongY(0.01);
	updateGL(); PRINT_OPENGL_ERROR();
}


void WidgetOpenGL::mousePressEvent(QMouseEvent *event)
{
	cam.xPrevious()=event->x();
	cam.yPrevious()=event->y();

	updateGL(); PRINT_OPENGL_ERROR();
}

void WidgetOpenGL::mouseMoveEvent(QMouseEvent *event)
{
	int const x=event->x();
	int const y=event->y();

	int const ctrl_pressed  = (event->modifiers() & Qt::ControlModifier);
	int const shift_pressed = (event->modifiers() & Qt::ShiftModifier);

	if(!ctrl_pressed && !shift_pressed && (event->buttons() & Qt::LeftButton))
		cam.rotation(x, y);
	if(!ctrl_pressed && !shift_pressed && (event->buttons() & Qt::RightButton))
	{
		cam.zoom(y);
		m_steps = y;
	}

	// Shift+Left button controls the window translation (left/right, bottom/up)
	float const dL=0.0001f*(1+10*fabs(cam.getDist()));
	if( !ctrl_pressed && shift_pressed && (event->buttons() & Qt::LeftButton) )
	{
		cam.goUp(dL*(y-cam.yPrevious()));
		cam.goRight(-dL*(x-cam.xPrevious()));
	}

	// Shift+Right button enables to translate forward/backward
	if( !ctrl_pressed && shift_pressed && (event->buttons() & Qt::RightButton) )
		cam.goForward(5.0f*dL*(y-cam.yPrevious()));

	cam.xPrevious()=x;
	cam.yPrevious()=y;

	updateGL(); PRINT_OPENGL_ERROR();
}

void WidgetOpenGL::wheelEvent(QWheelEvent * event)
{
	int numDegrees = event->delta() / 8;
	int numSteps = numDegrees;

	m_steps += numSteps;
	if (event->orientation() == Qt::Vertical)
		cam.zoomWheel(m_steps);
	event->accept();
}

string WidgetOpenGL::getFilename()
{
	return m_recordName;
}

string WidgetOpenGL::askForVTKFilename()
{
	return QFileDialog::getSaveFileName(this, "Recording of the displayed fibers as vtk: ", "/home/comercier/Bureau", "vtk files (*.vtk)").toStdString();
}
