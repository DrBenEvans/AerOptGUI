/*********************************************
**
**	Created on: 	28/12/2013 2013
**	Author: 	matt - Matt Edmunds
**	File:		main.cpp
**  This is the English version of this programme
**
**********************************************/

#include "MainWindow.h"
#include <iostream>
#include <QApplication>
#include <QDir>
#include <QDebug>

#include "DebugOutput.h"
#include "TreeView.h"
#include "Canvas.h"
#include "ProjectData.h"

//Programme entry point.
int main(int argc, char *argv[])
{
	//Initialises the QT library application event loops.
	QApplication a(argc, argv);
	a.setWindowIcon( QIcon(":/images/AerOpt.png") );

	//Application main interaction classes.
	DebugOutput& debugOutput = DebugOutput::Instance();
	ProjectData data;
	Canvas canvas(data);
	TreeView treeView(data, canvas);

	//Main window setup and show.
	MainWindow w(debugOutput, treeView, canvas);
	w.setWindowTitle("AerOpt");
	w.show();

	return a.exec();
}
