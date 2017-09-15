/*********************************************
**
**	Created on: 	28/12/2013 2013
**	Author: 	matt - Matt Edmunds
**	File:		main.cpp
**  This is the English version of this programme
**
**********************************************/

#include <iostream>
#include <QApplication>
#include <QDir>
#include <QDebug>

#include "DebugOutput.h"
#include "MainWindow.h"
#include "OptimisationModel.h"

//Programme entry point.
int main(int argc, char *argv[])
{
    //Initialises the QT library application event loops.
    QApplication app(argc, argv);
    app.setWindowIcon( QIcon(":/images/AerOpt.png") );
    //app.setAttribute(Qt::AA_DontUseNativeMenuBar);

    //Application main interaction classes.
    DebugOutput::Instance();
    qInfo() << " *** Welcome to AerOpt ***";
    qInfo() << " ***     Have a nice day     ***";

    // Setup optimisation model and selection model
    OptimisationModel* optimisationModel = new OptimisationModel();
    QItemSelectionModel* optimisationSelectionModel = new QItemSelectionModel(optimisationModel);
    optimisationModel->setSelectionModel(optimisationSelectionModel);

    //Main window setup and show.
    MainWindow w;
    w.setWindowTitle("AerOpt");
    w.setOptimisationModel(optimisationModel);
    w.show();

    return app.exec();
}
