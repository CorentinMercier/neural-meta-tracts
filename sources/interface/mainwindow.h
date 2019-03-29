#pragma once

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "widgetopengl.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    /** Quit the application */
    void action_quit();
    /** Set the Wireframe mode for the meshes */
    void action_wireframe();

    void action_fiber();
    void action_cylinder();
    void action_slider();
    void action_geometry_slider();
    void action_merge();
    void action_only_merge();
    void action_color_mode();
    void action_outliers();
    void action_corticalSurface();
    void action_save_vtk();

    void keyPressEvent(QKeyEvent *event);

private:

    /** Layout for the Window */
    Ui::MainWindow *ui;
    /** The OpenGL Widget */
    WidgetOpenGL *glWidget;
};

#endif // MAINWINDOW_H
