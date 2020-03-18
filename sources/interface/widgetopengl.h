#pragma once

#ifndef WIDGETOPENGL_H
#define WIDGETOPENGL_H

#include <glad/glad.h>
#include <GL/gl.h>

#include <QtOpenGL/QGLWidget>
#include <QMouseEvent>
#include <QFileDialog>
#include <QDialogButtonBox>

#include <cstring>
#include <iostream>

#include "../scene/scene.h"
#include "../opengl/openglutils.h"
#include "../opengl/camera.h"

using namespace std;

/** Qt Widget to render OpenGL scene */
class WidgetOpenGL : public QGLWidget
{
	Q_OBJECT

public:

	WidgetOpenGL(const QGLFormat& format, QGLWidget *parent = 0);
	~WidgetOpenGL();

	scene& get_scene();
	//Camera
	camera cam;

	string getFilename();

	string askForVTKFilename();

protected:

	/** Setup the OpenGL rendering mode */
	void initializeGL();
	/** The actual rendering function */
	void paintGL();
	/** Function called when the window is resized */
	void resizeGL(const int width, const int height);

	/** Function called a button of the mouse is pressed */
	void mousePressEvent(QMouseEvent *event);
	/** Function called when the wheel of the mouse is moved **/
	void wheelEvent(QWheelEvent *event);
	/** Function called when the mouse is moved */
	void mouseMoveEvent(QMouseEvent *event);
	/** Function called in a timer loop */
	void timerEvent(QTimerEvent *event);
	/** Function called when keyboard is pressed */
	void keyPressEvent(QKeyEvent *event);

private:

	/** Init the OpenGL rendering mode once at the beginning */
	void setup_opengl();
	/** Init Glew once at the beginning */
	void setup_glad();
	/** Print on the command line the actual version of the OpenGL Context */
	void print_current_opengl_context() const;
	/** All the content of the 3D scene */
	scene * m_scene;

	string m_filename;
	string m_recordName;
	string m_surface;
	string m_fuzzy;
	string m_setPoint1;
	string m_setPoint2;

	bool m_turnCam=false;

	int m_steps = -1325;
};

#endif // WIDGETOPENGL_H
