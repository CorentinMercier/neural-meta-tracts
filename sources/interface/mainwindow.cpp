#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent), ui(new Ui::MainWindow)
{
	//Setup window layout
	ui->setupUi(this);

	this->resize(1920*2/3,1280*2/3);

	//Create openGL context
	QGLFormat qglFormat;
	qglFormat.setVersion(4,5);
	qglFormat.setDoubleBuffer(false);
	qglFormat.setSampleBuffers(false);
	//Create OpenGL Widget renderer
	glWidget=new WidgetOpenGL(qglFormat);

	//Add the OpenGL Widget into the layout
	ui->layout_scene->addWidget(glWidget);

	//Connect slot and signals
	connect(ui->quit,SIGNAL(clicked()),this,SLOT(action_quit()));
	connect(ui->wireframe,SIGNAL(clicked()), this, SLOT(action_wireframe()));
	connect(ui->fibers,SIGNAL(clicked()), this, SLOT(action_fiber()));
	connect(ui->cylinders,SIGNAL(clicked()), this, SLOT(action_cylinder()));
	connect(ui->slider,SIGNAL(sliderMoved(int)), this, SLOT(action_slider()));
	connect(ui->geometrySlider,SIGNAL(sliderMoved(int)), this, SLOT(action_geometry_slider()));
	//connect(ui->Merge ,SIGNAL(clicked()), this, SLOT(action_merge()));
	connect(ui->outliers ,SIGNAL(clicked()), this, SLOT(action_outliers()));
	connect(ui->ColorMode,SIGNAL(clicked()), this, SLOT(action_color_mode()));
	connect(ui->corticalSurface, SIGNAL(clicked()), this, SLOT(action_corticalSurface()));
	connect(ui->SaveVTK,SIGNAL(clicked()),this,SLOT(action_save_vtk()));
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::action_quit()
{
	close();
}


void MainWindow::action_wireframe()
{
	bool const state_wireframe=ui->wireframe->isChecked();
	glWidget->get_scene().wireframe_mode=state_wireframe;
//	if (state_wireframe)
//	{
//		glWidget->get_scene().fiber_draw = false;
//		ui->fibers->setCheckState(Qt::Unchecked);
//	}
	glWidget->get_scene().valueChanged();
}

void MainWindow::action_fiber()
{
	bool const state = ui->fibers->isChecked();
	glWidget->get_scene().fiber_draw = state;
//    if (state)
//    {
//        glWidget->get_scene().cylinder_draw = false;
//        glWidget->get_scene().wireframe_mode = false;
//        ui->cylinders->setCheckState(Qt::Unchecked);
//        ui->wireframe->setCheckState(Qt::Unchecked);
//    }
	glWidget->get_scene().valueChanged();
}

void MainWindow::action_cylinder()
{
	bool const state = ui->cylinders->isChecked();
	glWidget->get_scene().cylinder_draw = state;
//    if (state)
//    {
//        glWidget->get_scene().fiber_draw = false;
//        ui->fibers->setCheckState(Qt::Unchecked);
//    }
	glWidget->get_scene().valueChanged();
}

void MainWindow::action_slider()
{
	int const value = ui->slider->sliderPosition();// value();
	string textToPrint = glWidget->get_scene().slider_set_value(value);
	ui->informations->setText(QString::fromStdString(textToPrint));
}

void MainWindow::action_geometry_slider()
{
	int const value = ui->geometrySlider->sliderPosition();// value();
	float val = (float)value/2550;
	string textToPrint = glWidget->get_scene().slider_fuzzy_change_value(val);
	ui->informations->setText(QString::fromStdString(textToPrint));
}

void MainWindow::action_corticalSurface()
{
	bool const state = ui->corticalSurface->isChecked();
	glWidget->get_scene().corticalSurface_draw = state;
}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
	int current=event->key();
	if (current==Qt::Key_Right)
	{
		string textToPrint=glWidget->get_scene().changeMultiscale(1);
		ui->informations->setText(QString::fromStdString(textToPrint));
	}
	if (current==Qt::Key_Left)
	{
		string textToPrint=glWidget->get_scene().changeMultiscale(-1);
		ui->informations->setText(QString::fromStdString(textToPrint));
	}
}

void MainWindow::action_merge()
{
	//bool const state = ui->Merge->isChecked();
	//glWidget->get_scene().merge_draw = state;
}

void MainWindow::action_outliers()
{
	bool const state = ui->outliers->isChecked();
	//glWidget->get_scene().outliers = state;
	glWidget->get_scene().outlierAction(state);
	glWidget->get_scene().valueChanged();
}

void MainWindow::action_only_merge()
{
	//bool const state = ui->onlyMerged->isChecked();
}

void MainWindow::action_color_mode()
{
	bool const state = ui->ColorMode->isChecked();
	glWidget->get_scene().changeColorState(state);
}

void MainWindow::action_save_vtk()
{
	glWidget->get_scene().saveVTK();
}
