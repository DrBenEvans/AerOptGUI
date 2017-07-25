/*********************************************
**
**	Created on: 	03/01/2014 2014
**	Author: 	matt - Matt Edmunds
**	File:		MainWindow.cpp
**
**********************************************/

#include "MainWindow.h"
#include "DebugOutput.h"
#include "TreeView.h"
#include "Canvas.h"
#include <iostream>

MainWindow::MainWindow(DebugOutput& debug, TreeView& tree, Canvas& canvas)
	: mDebugOutput(debug), mTreeView(tree), mCanvas(canvas)
{
	setupUi(this);

	this->resize(1280, 800);
	horizontalLayout = new QHBoxLayout(centralwidget);

	splitterHorizontal = new QSplitter(centralwidget);
	splitterHorizontal->setOrientation(Qt::Horizontal);

	splitterHorizontal->addWidget(&mTreeView);

	widgetSplitter = new QWidget(splitterHorizontal);

	verticalLayout = new QVBoxLayout(widgetSplitter);
	verticalLayout->setSpacing(0);
	verticalLayout->setContentsMargins(0, 0, 0, 0);

	splitterVertical = new QSplitter(widgetSplitter);
	splitterVertical->setOrientation(Qt::Vertical);

	mCanvas.setMinimumWidth(800);
	mCanvas.setMinimumHeight(480);

	splitterVertical->addWidget(&mCanvas);
	splitterVertical->addWidget(&mDebugOutput);

	verticalLayout->addWidget(splitterVertical);
	splitterHorizontal->addWidget(widgetSplitter);
	horizontalLayout->addWidget(splitterHorizontal);

	splitterHorizontal->setStretchFactor(0, 0);
	splitterHorizontal->setStretchFactor(1, 16777215);

	splitterVertical->setStretchFactor(0, 16777215);
	splitterVertical->setStretchFactor(1, 0);
}

MainWindow::~MainWindow()
{
	layout()->removeWidget( &mTreeView );
	mTreeView.setParent(nullptr);

	layout()->removeWidget( &mCanvas );
	mCanvas.setParent(nullptr);

	layout()->removeWidget( &mDebugOutput );
	mDebugOutput.setParent(nullptr);
}

void MainWindow::on_actionStop_Mesher_triggered()
{
	mTreeView.stopMesher();
}

void MainWindow::on_actionStop_AerOpt_triggered()
{
	mTreeView.stopAerOpt();
}
