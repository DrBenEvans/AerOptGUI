/*********************************************
**
**	Created on: 	03/01/2014 2014
**	Author: 	matt - Matt Edmunds
**	File:		MainWindow.h
**
**********************************************/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "ui_MainWindow.h"
#include <QHBoxLayout>
#include <QSplitter>

class DebugOutput;
class TreeView;
class Canvas;

class MainWindow : public QMainWindow, private Ui::MainWindow
{
	Q_OBJECT
	
public:
	explicit MainWindow(DebugOutput& debug, TreeView& tree, Canvas& canvas);
	MainWindow() = delete;
	~MainWindow();

	QSplitter *splitterVertical;
	QVBoxLayout *verticalLayout;
	QWidget *widgetSplitter;
	QSplitter *splitterHorizontal;
	QHBoxLayout *horizontalLayout;


	DebugOutput& mDebugOutput;
	TreeView& mTreeView;
	Canvas& mCanvas;
private slots:
	void on_actionStop_Mesher_triggered();
	void on_actionStop_AerOpt_triggered();
};

#endif // MAINWINDOW_H
