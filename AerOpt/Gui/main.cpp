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
#include "AppController.h"
#include "Canvas.h"
#include "OptimisationRun.h"
#include "PlotterDialog.h"

//Programme entry point.
int main(int argc, char *argv[])
{
	//Initialises the QT library application event loops.
    QApplication app(argc, argv);
    app.setWindowIcon( QIcon(":/images/AerOpt.png") );
    //app.setAttribute(Qt::AA_DontUseNativeMenuBar);

	//Application main interaction classes.
    DebugOutput& debugOutput = DebugOutput::Instance();
    qInfo() << " *** Welcome to AerOpt ***";
    qInfo() << " ***     Have a nice day     ***";

	//Main window setup and show.
    MainWindow w;
    w.setWindowTitle("AerOpt");
    w.show();

    return app.exec();
}
