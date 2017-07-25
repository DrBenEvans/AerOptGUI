/*********************************************
**
**	Created on: 	08/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		TreeView.cpp
**
**********************************************/

#include <QDebug>
#include <QFileDialog>
#include <QStringList>
#include <QString>
#include <QTimer>
#include <QEventLoop>

#include <fstream>
#include <map>
#include <utility>
#include <stdexcept>
#include <algorithm>
#include <tuple>
#include <cmath>

#include "TreeView.h"
#include "Canvas.h"
#include "Menu.h"
#include "ProjectData.h"
#include "Enumerations.h"
#include "ObjectiveDialog.h"
#include "BoundaryDialog.h"
#include "OptimiserDialog.h"
#include "MeshDialog.h"
#include "PlotterDialog.h"

TreeView::TreeView(ProjectData& data, Canvas& canvas) : mCanvas(canvas), mData(data)
{
	setupUi(this);

	myMeshProcess.setParent(this);
	connect (&myMeshProcess, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(meshingFinished(int,QProcess::ExitStatus)));
	connect (&myMeshProcess, SIGNAL(started()), this, SLOT(meshingStarted()));
//	connect (myMeshProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(meshOutput()));
	connect (&myMeshProcess, SIGNAL(readyReadStandardError()), this, SLOT(meshError()));
	mProjectDirectory = "";

	myOptProcess.setParent(this);
	connect (&myOptProcess, SIGNAL(finished(int,QProcess::ExitStatus)),
			 this, SLOT(optimiserFinished(int,QProcess::ExitStatus)));
	connect (&myOptProcess, SIGNAL(started()), this, SLOT(optimiserStarted()));
//	connect (myOptProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(processOutput()));
//	connect (myOptProcess, SIGNAL(readyReadStandardError()), this, SLOT(processError()));

	myDirWatcher.setParent(this);
//	QObject::connect(&myDirWatcher, SIGNAL(fileChanged(const QString&)), this, SLOT(readFitness(const QString&)));
	QObject::connect(&myDirWatcher, SIGNAL(directoryChanged(const QString&)), this, SLOT(readDirectory(const QString&)));

	treeWidget->setIndentation(20);
	treeWidget->setRootIsDecorated(true);
	treeWidget->setUniformRowHeights(true);
	treeWidget->setItemsExpandable(true);
	treeWidget->setAnimated(true);
	treeWidget->setExpandsOnDoubleClick(true);
	treeWidget->setColumnCount(3);
	treeWidget->setContextMenuPolicy(Qt::CustomContextMenu);

	QStringList headers;
    headers << "AerOpt Settings";
    headers << "Set";
	headers << " ";
	treeWidget->setHeaderLabels(headers);
	treeWidget->setColumnWidth(0, 205);
	treeWidget->setColumnWidth(1, 33);

	mRoot = new QTreeWidgetItem(treeWidget, Enum::TreeType::ROOT);
	mRoot->setText(0, "Aerofoil Project");
    mRoot->setText(1, "No");
	treeWidget->expandItem(mRoot);
	treeWidget->setCurrentItem(mRoot);
	mCurrentNode = mRoot;

	mMenus.clear();
	mMenus.insert( std::pair<int, Menu*>(Enum::TreeType::DEFAULT,   new Menu(mData, mCanvas, Enum::TreeType::DEFAULT,   this)) );
	mMenus.insert( std::pair<int, Menu*>(Enum::TreeType::ROOT,      new Menu(mData, mCanvas, Enum::TreeType::ROOT,      this)) );
	mMenus.insert( std::pair<int, Menu*>(Enum::TreeType::DATA,      new Menu(mData, mCanvas, Enum::TreeType::DATA,      this)) );
	mMenus.insert( std::pair<int, Menu*>(Enum::TreeType::MESH,      new Menu(mData, mCanvas, Enum::TreeType::MESH,      this)) );
	mMenus.insert( std::pair<int, Menu*>(Enum::TreeType::FUNCTION,  new Menu(mData, mCanvas, Enum::TreeType::FUNCTION,  this)) );
	mMenus.insert( std::pair<int, Menu*>(Enum::TreeType::BOUNDARY,  new Menu(mData, mCanvas, Enum::TreeType::BOUNDARY,  this)) );
	mMenus.insert( std::pair<int, Menu*>(Enum::TreeType::OPTIMISER, new Menu(mData, mCanvas, Enum::TreeType::OPTIMISER, this)) );
	mMenus.insert( std::pair<int, Menu*>(Enum::TreeType::RUNTIME,   new Menu(mData, mCanvas, Enum::TreeType::RUNTIME,   this)) );

	QObject::connect(treeWidget, SIGNAL( customContextMenuRequested(const QPoint &)), this, SLOT(showContextMenu(const QPoint &)));
	mMenusSet = true;

	//Determin application directory and binary locations.
	int index = 0;
	QString appPath;
#ifdef Q_OS_UNIX
	// do fancy unix stuff here
	QString mesherExe = "/AerOpt/FLITE/Mesher/MeshGenerator";
	QString aeroptExe = "/AerOpt/FLITE/AerOpt";
#endif
#ifdef Q_OS_WIN32
	// do windows stuff here
	QString mesherExe = "/AerOpt/FLITE/Mesher/MeshGenerator.exe";
	QString aeroptExe = "/AerOpt/FLITE/AerOpt.exe";
#endif


	appPath = QCoreApplication::applicationFilePath();
//	qDebug() << appPath;
	appPath = QDir::fromNativeSeparators(QCoreApplication::applicationFilePath());

	index = appPath.lastIndexOf("/AerOpt/");

	if (index != -1)
	{
		//Determin app path
		appPath.remove(index, appPath.size() - index);
		appPath = QDir::toNativeSeparators(appPath);
//		qDebug() << appPath;
		mAppPath = appPath;

		//Determin mesher path
		mesherExe = QDir::toNativeSeparators(mesherExe);
		mesherExe.prepend( appPath );
//		qDebug() << mesherExe;
		mMesherPath = mesherExe;

		//Determin flite path << lol
		aeroptExe = QDir::toNativeSeparators(aeroptExe);
		aeroptExe.prepend( appPath );
//		qDebug() << aeroptExe;
		mAerOptPath = aeroptExe;

        QFileInfo checkMesher(mMesherPath);
		QFileInfo checkAerOpt(mAerOptPath);

        bool c = true;
        c &= checkMesher.exists();
		c &= checkMesher.isFile();
        if (!c)
        {
            qCritical() << "Mesher executables not found at: " << mMesherPath;
        }

        c = true;
		c &= checkAerOpt.exists();
		c &= checkAerOpt.isFile();
        if (!c)
        {
            qCritical() << "AerOpt executables not found at: " << mAerOptPath;
        }


	}
	else
	{
		qCritical() << "Application root directory could not be established!";
	}

	//Get source folders (input and output)
	QDir path(mAppPath);
	QDir appPathOut = QDir(path.absolutePath() + QDir::separator() + "AerOpt/FLITE/Output_Data").absolutePath();
	QDir appPathIn = QDir(path.absolutePath() + QDir::separator() + "AerOpt/FLITE/Input_Data").absolutePath();
	appPathOut = QDir::toNativeSeparators(appPathOut.path());
	appPathIn = QDir::toNativeSeparators(appPathIn.path());

	if (!appPathOut.exists())
	{
		appPathOut.mkpath(appPathOut.path());
	}

	if (!appPathIn.exists())
	{
		appPathIn.mkpath(appPathIn.path());
	}

	sGenNo = 0;

	//Set Defaults.
	addProfileObject();
	addMeshObject();
	addFunctionObject();
	addBoundaryObject();
	addOptimiserObject();
	addRuntimeObject();
    qInfo() << " *** Welcome to AerOpt ***";
    qInfo() << " ***     Have a nice day     ***";

	mPlotter = new PlotterDialog(this);
	mPlotter->hide();
}

TreeView::~TreeView()
{
	myMeshProcess.kill();
	myOptProcess.kill();
}

//Public slots

//Main menus

void TreeView::addProfileObject()
{
	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();
	if (selectedList.size() > 0) mCurrentNode = selectedList.at(0);

	for (QTreeWidgetItem* node : selectedList)
	{
		bool b = false;
		int qty = node->childCount();

		for (int i = 0; i < qty; ++i)
		{
			b |= node->child(i)->type() == Enum::TreeType::DATA;
		}

		if ( b )
		{
            qDebug() << "Item 'Profile' already exists: ";
		}
		else
		{
			QTreeWidgetItem* item = new QTreeWidgetItem(node, Enum::TreeType::DATA);
            QString name = "Profile";
			item->setText(0, name);
            item->setText(1, "No");
			mData.setProfile(false);
			treeWidget->expandItem(item);
		}
	}
}

void TreeView::addMeshObject()
{
	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();
	if (selectedList.size() > 0) mCurrentNode = selectedList.at(0);

	for (QTreeWidgetItem* node : selectedList)
	{
		bool b = false;
		int qty = node->childCount();

		for (int i = 0; i < qty; ++i)
		{
			b |= node->child(i)->type() == Enum::TreeType::MESH;
		}

		if ( b )
		{
            qDebug() << "Item 'Mesher' already exists: ";
		}
		else
		{
			QTreeWidgetItem* item = new QTreeWidgetItem(node, Enum::TreeType::MESH);
            QString name = "Mesher";
			item->setText(0, name);
            item->setText(1, "No");
			mData.setMesh(false);
			treeWidget->expandItem(item);
		}
	}
}

void TreeView::addFunctionObject()
{
	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();
	if (selectedList.size() > 0) mCurrentNode = selectedList.at(0);

	for (QTreeWidgetItem* node : selectedList)
	{
		bool b = false;
		int qty = node->childCount();

		for (int i = 0; i < qty; ++i)
		{
			b |= node->child(i)->type() == Enum::TreeType::FUNCTION;
		}

		if ( b )
		{
            qInfo() << "Item 'Fitness Function' already exists";
		}
		else
		{
			QTreeWidgetItem* item = new QTreeWidgetItem(node, Enum::TreeType::FUNCTION);
            QString name = "Fitness Function";
			item->setText(0, name);
            item->setText(1, "No");
			mData.setFunction(false);
			treeWidget->expandItem(item);
		}
	}
}

void TreeView::addBoundaryObject()
{
	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();
	if (selectedList.size() > 0) mCurrentNode = selectedList.at(0);

	for (QTreeWidgetItem* node : selectedList)
	{
		bool b = false;
		int qty = node->childCount();

		for (int i = 0; i < qty; ++i)
		{
			b |= node->child(i)->type() == Enum::TreeType::BOUNDARY;
		}

		if ( b )
		{
            qInfo() << "Item 'Flow Conditions' already exists";
		}
		else
		{
			QTreeWidgetItem* item = new QTreeWidgetItem(node, Enum::TreeType::BOUNDARY);
            QString name = "Flow Conditions";
			item->setText(0, name);
            item->setText(1, "No");
			mData.setBoundary(false);
			treeWidget->expandItem(item);
		}
	}
}

void TreeView::addOptimiserObject()
{
	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();
	if (selectedList.size() > 0) mCurrentNode = selectedList.at(0);

	for (QTreeWidgetItem* node : selectedList)
	{
		bool b = false;
		int qty = node->childCount();

		for (int i = 0; i < qty; ++i)
		{
			b |= node->child(i)->type() == Enum::TreeType::OPTIMISER;
		}

		if ( b )
		{
            qInfo() << "Item 'Optimiser Parameters' already exists";
		}
		else
		{
			QTreeWidgetItem* item = new QTreeWidgetItem(node, Enum::TreeType::OPTIMISER);
			QString name = "Optimiser Parameters";
			item->setText(0, name);
            item->setText(1, "No");
			mData.setOptimiser(false);
			treeWidget->expandItem(item);
		}
	}
}

void TreeView::addRuntimeObject()
{
	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();
	if (selectedList.size() > 0) mCurrentNode = selectedList.at(0);

	for (QTreeWidgetItem* node : selectedList)
	{
		bool b = false;
		int qty = node->childCount();

		for (int i = 0; i < qty; ++i)
		{
			b |= node->child(i)->type() == Enum::TreeType::RUNTIME;
		}

		if ( b )
		{
            qInfo() << "Item 'Runtime Monitoring' already exists";
		}
		else
		{
			QTreeWidgetItem* item = new QTreeWidgetItem(node, Enum::TreeType::RUNTIME);
			QString name = "Runtime Monitoring";
			item->setText(0, name);
            item->setText(1, "No");
			mData.setRunTime(false);
			treeWidget->expandItem(item);
		}
	}
}

void TreeView::deleteObject()
{
	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();

	for (QTreeWidgetItem* node : selectedList)
	{
        qInfo() << "Deleted " << node->text(0);
		delete node;
	}
}

void TreeView::clearProject()
{
    qInfo() << "Project data cleared!";
    mRoot->setText(1, "No");

    //All sub items set to 'No'
	for(int i = 0; i < mRoot->childCount(); ++i)
	{
        mRoot->child(i)->setText(1, "No");
	}

	mData.clearProject();
	mCanvas.update();
}

void TreeView::loadProject()
{
	QDir projPath;

#ifdef Q_OS_UNIX
	// do fancy unix stuff
	projPath = QFileDialog::getExistingDirectory(this, "Select Project Directory", QDir::homePath()+"/Documents/Projects/AerOptProject/");
#endif
#ifdef Q_OS_WIN32
	// do windows stuff here
	projPath = QFileDialog::getExistingDirectory(this, "Select Project Directory", QDir::homePath());
#endif

	mProjectDirectory = projPath.path();
	mData.setProjectPathSet(true);
    mRoot->setText(1, "Yes");

    qInfo() << "Project directory set: " << mProjectDirectory;
}

//Sub menus

void TreeView::importProfile()
{
	if (!mData.projectPathSet())
	{
		qWarning() << "The Project directory is not set! Please set the Project directory.";
		return;
	}

	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();
	if (selectedList.size() > 0) mCurrentNode = selectedList.at(0);

	bool success = true;
	QString fileName;
	QStringList fileNames;


#ifdef Q_OS_UNIX
	// do fancy unix stuff
	fileNames.append(
        QFileDialog::getOpenFileName(this, "Select Profile File", QDir::homePath()+"/Documents/Projects/AerOptProject/", "Profile Files (*.prf)")
					);
#endif
#ifdef Q_OS_WIN32
	// do windows stuff here
	fileNames.append(
        QFileDialog::getOpenFileName(this, "Select Profile File", QDir::homePath(), "Profile Files (*.prf)")
					);
#endif

	if (fileNames.size() > 0)
	{
		for (const QString &f: fileNames)
		{
			fileName = f;
		}

        qInfo() << "File selected: " << fileName;

		success &= loadProfile(fileName.toStdString(), mData);

		if (success)
		{
			//Set text OK here!
			for (QTreeWidgetItem* node : selectedList)
			{
                node->setText(1, "Yes");
				mData.setProfile(true);
				mData.setRenderProfile(true);
				mData.setRenderMesh(false);
			}
            qInfo() << "File successfully loaded.";
		}
		else
		{
			//Set text Not OK here!
			for (QTreeWidgetItem* node : selectedList)
			{
                node->setText(1, "No");
				mData.setProfile(false);
				mData.setRenderProfile(false);
				mData.setRenderMesh(true);
			}
			qWarning() << "File failed to load correctly.";
		}
	}
	else
	{
		//Set text Not OK here!
		for (QTreeWidgetItem* node : selectedList)
		{
            node->setText(1, "No");
			mData.setProfile(false);
			mData.setRenderProfile(false);
			mData.setRenderMesh(true);
		}
		mData.clearProfile();
		qWarning() << "Profile data not imported!";
	}
	mCanvas.update();
}

void TreeView::runMesher()
{
	sGenNo = 0;
	//Clear mesh data first
    //Set menu to 'No' here
	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();
	if (selectedList.size() > 0) mCurrentNode = selectedList.at(0);
	for (QTreeWidgetItem* node : selectedList)
	{
        node->setText(1, "No");
		mData.setMesh(false);
		mData.setRenderMesh(false);
		mData.setRenderProfile(true);
	}
	mData.clearMesh();

	if (!mData.profile())
	{
		qWarning() << "Profile not loaded! Please load a valid profile.";
		mData.setRenderProfile(false);
		mCanvas.update();
		return;
	}

	MeshDialog diag(mData, this);

	diag.show();
	diag.exec();

	bool r = true;
	QDir meshPath(mProjectDirectory);
	QString meshInFile;
	QString meshBacFile;
	QString meshGeoFile;
	QString meshDatFile;

	QString meshWorkDir = mAppPath;
	meshWorkDir += "/AerOpt/FLITE/Mesher/";
	meshWorkDir = QDir::toNativeSeparators(meshWorkDir);

	meshInFile = QDir(meshPath.absolutePath() + QDir::separator() + meshPath.dirName() + ".in").absolutePath();
	meshInFile = QDir::toNativeSeparators(meshInFile);

	meshBacFile = QDir(meshPath.absolutePath() + QDir::separator() + meshPath.dirName() + ".bac").absolutePath();
	meshBacFile = QDir::toNativeSeparators(meshBacFile);

	meshGeoFile = QDir(meshPath.absolutePath() + QDir::separator() + meshPath.dirName() + ".geo").absolutePath();
	meshGeoFile = QDir::toNativeSeparators(meshGeoFile);

	meshDatFile = QDir(meshPath.absolutePath() + QDir::separator() + meshPath.dirName() + ".dat").absolutePath();
	meshDatFile = QDir::toNativeSeparators(meshDatFile);

	r &= createInputFile(meshInFile.toStdString(),
						 meshBacFile.toStdString(),
						 meshGeoFile.toStdString(),
						 meshDatFile.toStdString());
	r &= QFile::exists(meshInFile);

	r &= createBacFile(meshBacFile.toStdString());
	r &= QFile::exists(meshBacFile);

	r &= createGeoFile(meshGeoFile.toStdString(), mData);
	r &= QFile::exists(meshGeoFile);

	if (r)
	{
		mMenusSet = false;

		myMeshProcess.setWorkingDirectory( meshWorkDir );
		myMeshProcess.setStandardInputFile( meshInFile );
		myMeshProcess.start( mMesherPath );
	}
}

void TreeView::setObjective()
{
	bool r = true;
	if (!mData.mesh())
	{
		qWarning() << "Please first generate an initial mesh!";
		return;
	}
	ObjectiveDialog diag(mData, this);

	diag.show();
	diag.exec();

	r &= mData.function();
    if (r) mCurrentNode->setText(1, "Yes");
}

void TreeView::setBoundary()
{
	bool r = true;
	if (!mData.function())
	{
		qWarning() << "Please first set the objective function!";
		return;
	}
	BoundaryDialog diag(mData, this);

	diag.show();
	diag.exec();

	r &= mData.boundary();
    if (r) mCurrentNode->setText(1, "Yes");
}

void TreeView::setOptimiser()
{
	bool r = true;
	if (!mData.boundary())
	{
        qWarning() << "Please first set the flow conditions!";
		return;
	}
	OptimiserDialog diag(mData, this);

	diag.show();
	diag.exec();

	r &= mData.optimiser();
    if (r) mCurrentNode->setText(1, "Yes");
}

void TreeView::runAerOpt()
{
	bool r = true;

	if (!mData.optimiser())
	{
		qWarning() << "Please first set the optimiser parameters!";
		return;
	}
	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();
	if (selectedList.size() > 0) mCurrentNode = selectedList.at(0);

	QString AerOpt = mAerOptPath;

	QString AerOptWorkDir = mAppPath;
	AerOptWorkDir += "/AerOpt/FLITE/";
	AerOptWorkDir = QDir::toNativeSeparators(AerOptWorkDir);

	QString AerOptInFile = AerOptWorkDir;
	AerOptInFile += "Input_Data/AerOpt_InputParameters.txt";
	AerOptInFile = QDir::toNativeSeparators(AerOptInFile);

	QString AerOptNodeFile = AerOptWorkDir;
	AerOptNodeFile += "Input_Data/Control_Nodes.txt";
	AerOptNodeFile = QDir::toNativeSeparators(AerOptNodeFile);

	QString path = mAppPath + "/AerOpt/FLITE/Output_Data/";
	path = QDir::toNativeSeparators(path);

	QDir meshPath(mProjectDirectory);
	QString meshDatFile;
	meshDatFile = QDir(meshPath.absolutePath() + QDir::separator() + meshPath.dirName() + ".dat").absolutePath();
	meshDatFile = QDir::toNativeSeparators(meshDatFile);

	//Empty the input and output directories
	QString inFolder = mAppPath;
	inFolder += "/AerOpt/FLITE/Input_Data";
	inFolder = QDir::toNativeSeparators(inFolder);
	r &= emptyFolder(inFolder);

	QString outFolder = mAppPath;
	outFolder += "/AerOpt/FLITE/Output_Data";
	outFolder = QDir::toNativeSeparators(outFolder);
	r &= emptyFolder(outFolder);

	//Copy mesh file as input to AerOpt from project directory to input directory
	QString dest = mAppPath;
	dest += "/AerOpt/FLITE/Input_Data/Mesh.dat";
	dest = QDir::toNativeSeparators(dest);
	r &= copyFile(meshDatFile, dest);

	//Set input files in input directory
	r &= createAerOptInFile(AerOptInFile.toStdString(), mData);
	r &= createAerOptNodeFile(AerOptNodeFile.toStdString(), mData);

	sGenNo = 0;
	mData.resetBoundary();

	//then run aeropt
	if (r)
	{
		mMenusSet = false;
		myDirWatcher.addPath(path);
		myOptProcess.setWorkingDirectory( AerOptWorkDir );
		myOptProcess.start( AerOpt );
	}
	else
	{
		qWarning() << "Process Aborted!";
	}
}

void TreeView::stopMesher()
{
	myMeshProcess.kill();
    qInfo() << "Any running mesh jobs have been stopped!";

	mMenusSet = true;
}

void TreeView::stopAerOpt()
{
	myOptProcess.kill();
    qInfo() << "Any running AerOpt jobs have been stopped!";

	mMenusSet = true;
}

void TreeView::showGraph()
{
	mPlotter->show();
	mPlotter->adjustSize();
}

//Private slots

void TreeView::showContextMenu(const QPoint &point)
{

	if (!mMenusSet) return;

	QList<QTreeWidgetItem*> selectedList = treeWidget->selectedItems();

	if (treeWidget->indexAt(point).isValid())
	if (selectedList.size() > 0)
	{
		mCurrentNode = selectedList.at(0);
		try
		{
			mMenus.at( selectedList.at(0)->type() )->exec( treeWidget->mapToGlobal(point) );
		}
		catch (const std::out_of_range& oor)
		{
			qCritical() << "Out of Range error: " << oor.what()
						<< " Location: 'void TreeView::showContextMenu(const QPoint &point)'";
		}
	}
}

void TreeView::meshOutput()
{
	QByteArray byteArray = myMeshProcess.readAllStandardOutput();
	QStringList strLines = QString(byteArray).split("\n");

	foreach (QString line, strLines){
        qInfo() << line;
	}
}

void TreeView::meshError()
{
	QByteArray byteArray = myMeshProcess.readAllStandardError();
	QStringList strLines = QString(byteArray).split("\n");

	foreach (QString line, strLines){
		qCritical() << line;
	}
}

void TreeView::processOutput()
{
	QByteArray byteArray = myOptProcess.readAllStandardOutput();
	QStringList strLines = QString(byteArray).split("\n");

	foreach (QString line, strLines){
        qInfo() << line;
	}
}

void TreeView::processError()
{
	QByteArray byteArray = myOptProcess.readAllStandardError();
	QStringList strLines = QString(byteArray).split("\n");

	foreach (QString line, strLines){
		qCritical() << line;
	}
}

void TreeView::meshingStarted()
{
    qInfo() << "Process started normally";
}

void TreeView::meshingFinished(int exitCode, QProcess::ExitStatus exitStatus)
{
	if (exitStatus == QProcess::NormalExit && exitCode == 0)
	{
		bool r = true;

        qInfo() << "Process finished normally";

		//Set and/or check existance of output data file

		QDir meshPath(mProjectDirectory);
		QString meshDatFile;
		meshDatFile = QDir(meshPath.absolutePath() + QDir::separator() + meshPath.dirName() + ".dat").absolutePath();
		meshDatFile = QDir::toNativeSeparators(meshDatFile);

		r &= myMeshProcess.exitCode() == 0;
		r &= QFile::exists(meshDatFile);

		//Now load points into data object, ready for canvas to render.
		if (r)
		{
			r &= loadMeshProfile(sGenNo, meshDatFile.toStdString(), mData);
			r &= loadMesh(meshDatFile.toStdString(), mData);
		}

        //When done set menu to 'Yes'
		if (r)
		{
            mCurrentNode->setText(1, "Yes");
			mData.setMesh(true);
			mData.setRenderMesh(true);
			mData.setRenderProfile(false);
            qInfo() << "Mesh successfully created.";
		}
		else
		{
            mCurrentNode->setText(1, "No");
			mData.setMesh(false);
			mData.setRenderMesh(false);
			mData.setRenderProfile(true);
			qWarning() << "Mesh not created!";
		}
		mCanvas.update();
	}
	else
	{
		qCritical() << "Process finished abnormally with exit code: " << exitCode;
		//Definately set Y/N in treeview here!!
	}
	mMenusSet = true;
}

void TreeView::optimiserStarted()
{
    qInfo() << "Process started normally";
	mPlotter->clearData();
	mPlotter->show();
}

void TreeView::readDirectory(const QString& path)
{
	//Read initial fitness
	if (sGenNo == 0)
	{
		bool b = true;
		QString fit = path;
		fit += "Fitness_";
		fit += QString::number(sGenNo);
		fit += ".txt";
		QFileInfo checkFit(fit);
		b &= checkFit.exists();
		b &= checkFit.isFile();
		if (b)
		{
			readFitness(fit);
			sGenNo++;
		}
	}

	//When meshfiles exists, load it.
	//get list of .dat files and trim "Geometry_" and ".dat"
	//then compare number with previous iteration.
	bool c = true;

	//Fitness file
	QString fit = path;
	fit += "Fitness_";
	fit += QString::number(sGenNo);
	fit += ".txt";

	//Mesh file
	QString mesh = path;
	mesh += "Geometry_";
	mesh += QString::number(sGenNo);
	mesh += ".dat";

	//Results file
	QString results = path;
	results += "Geometry_";
	results += QString::number(sGenNo);
	results += ".resp";

	//Load Mesh + Results
	QFileInfo checkFit(fit);
	QFileInfo checkMesh(mesh);
	QFileInfo checkResults(results);

	c &= checkFit.exists();
	c &= checkFit.isFile();
	c &= checkMesh.exists();
	c &= checkMesh.isFile();
	c &= checkResults.exists();
	c &= checkResults.isFile();

	if (c)
	{
		loadMeshProfile(sGenNo, mesh.toStdString(), mData);
		loadMesh(mesh.toStdString(), mData);
		loadResults(results.toStdString(), mData);
		readFitness(fit);

		//Output .prf file here
		QString prf = path;
		prf += "Profile_";
		prf += QString::number(sGenNo);
		prf += ".prf";
		saveCurrentProfile(prf, mData);

		sGenNo++;
	}
	mCanvas.update();
}

void TreeView::readFitness(const QString& path)
{
	bool r = true;
	std::string line = "";
	QString output = "";
	QString plotData = "";

	std::ifstream infile(path.toStdString(), std::ifstream::in);
	r &= infile.is_open();

	if (r)
	{
		std::getline(infile, line);

		output = QString::fromStdString(line);
		output = output.trimmed();
		plotData = output;
		output.prepend( " >  " );
		output.prepend( "Generation: " );
        qInfo() << output;


		//Extract from plotData gen no. and then nest values
		QStringList list1 = plotData.split(" ", QString::SkipEmptyParts);
		QVector<double> nests;
		int genNo = 0;
		bool ok;
		for (int i = 0; i < list1.size(); ++i)
		{
			if (i == 0) genNo = list1.at(i).toInt(&ok);
			else nests.push_back( list1.at(i).toDouble(&ok) );

			if (genNo == 0 && i == 1)
			{
				int cols = mData.noAgents() -1;
				for (int j = 0; j < cols; ++j)
				{
					nests.push_back( nests.first() );
				}
			}
		}

		//set plotter with extracted data
		if (ok)
		{
			genNo += 1;
			mPlotter->setData(genNo, nests);
		}

		//Also need to clear it someware too!
		//On start Optimiser
		//And on clear project

	}
	infile.close();
}

void TreeView::optimiserFinished(int exitCode, QProcess::ExitStatus exitStatus)
{

	//Windows bug!! re read last file set if exists.
	QEventLoop loop;
	QTimer::singleShot(1000, &loop, SLOT(quit()));
	loop.exec();

	QDir lastPath(mAppPath);
	lastPath = QDir(lastPath.absolutePath() + QDir::separator() + "AerOpt/FLITE/Output_Data").absolutePath();
	lastPath = QDir::toNativeSeparators(lastPath.path());
	readDirectory(lastPath.path());

	QEventLoop loop2;
	QTimer::singleShot(2000, &loop2, SLOT(quit()));
	loop2.exec();

	if (exitStatus == QProcess::NormalExit && exitCode == 0)
	{
        qInfo() << "Process finished normally";

		bool r = true;

		r &= myOptProcess.exitCode() == 0;

        //When done set menu to 'Yes'
        if (r) mCurrentNode->setText(1, "Yes");
        else mCurrentNode->setText(1, "No");
	}
	else
	{
		qCritical() << "Process finished abnormally with exit code: " << exitCode;
		//Definately set Y/N in treeview here!!
        mCurrentNode->setText(1, "No");
	}

	//Move files to project folder here!
	bool r = true;

	//Get source folders (input and output)
	QDir appPath(mAppPath);
	QDir appPathOut = QDir(appPath.absolutePath() + QDir::separator() + "AerOpt/FLITE/Output_Data").absolutePath();
	QDir appPathIn = QDir(appPath.absolutePath() + QDir::separator() + "AerOpt/FLITE/Input_Data").absolutePath();
	appPathOut = QDir::toNativeSeparators(appPathOut.path());
	appPathIn = QDir::toNativeSeparators(appPathIn.path());

	//Get destination folders (input and output)
	QDir projPath(mProjectDirectory);
	QDir projPathOut = QDir(projPath.absolutePath() + QDir::separator() + "Output_Data").absolutePath();
	QDir projPathIn = QDir(projPath.absolutePath() + QDir::separator() + "Input_Data").absolutePath();
	projPathOut = QDir::toNativeSeparators(projPathOut.path());
	projPathIn = QDir::toNativeSeparators(projPathIn.path());

	//Copy Output folders
	r &= copyFolder(appPathOut.path(), projPathOut.path());

	//Copy Input folders
	r &= copyFolder(appPathIn.path(), projPathIn.path());

	//delete cruft(2D*, Delaunay*, FileCreateDir*)
	QString cruft = QDir(appPath.path() + QDir::separator() + "AerOpt/FLITE/").absolutePath();
	cruft = QDir::toNativeSeparators(cruft);

	QFile::remove(QDir::toNativeSeparators(cruft + QDir::separator() + "Delaunay_elements.txt"));
	QFile::remove(QDir::toNativeSeparators(cruft + QDir::separator() + "Delaunay_nodes.txt"));
	QFile::remove(QDir::toNativeSeparators(cruft + QDir::separator() + "FileCreateDir.bat"));
	QFile::remove(QDir::toNativeSeparators(cruft + QDir::separator() + "FileCreateDir.scr"));

	QStringList filters;
	filters << "2D*";
	QStringList folders = QDir(cruft).entryList(filters, QDir::Dirs);

	foreach(QString dir, folders)
	{
		QString temp = cruft + QDir::separator() + dir;
		temp = QDir::toNativeSeparators(temp);
		r &= removeFolder( temp );
	}

	if (!r) qWarning() << "Something went wrong with file cleanup!";

	mCanvas.update();
	myDirWatcher.removePaths( myDirWatcher.directories() );
	mMenusSet = true;
}


//Private functions

bool TreeView::loadProfile(const std::string& filePath, ProjectData& data)
{
	bool r = true;
	float x, y;

	std::ifstream infile(filePath, std::ifstream::in);
	r &= infile.is_open();

	if (r)
	{
		data.clearProfile();
		while (infile >> x >> y)
		{
			data.addPoint(x, y);
		}
	}
	infile.close();

	r &= data.checkProfileIntegrity();

	if (!r) data.clearProfile();

	return r;
}

bool TreeView::createInputFile(const std::string& meshInFile,
							   const std::string& meshBacFile,
							   const std::string& meshGeoFile,
							   const std::string& meshDatFile)
{
	bool r = true;

	std::ofstream outfile(meshInFile, std::ofstream::out);
	r &= outfile.is_open();

	if (r)
	{
		outfile << meshGeoFile << std::endl;
		outfile << meshBacFile << std::endl;
		outfile << meshDatFile << std::endl;
		outfile << "0" << std::endl;
		outfile << "1" << std::endl;
        outfile << "N" << std::endl;
		outfile << "1" << std::endl;
		outfile << "1" << std::endl;
		outfile << "1" << std::endl;
		outfile << "1" << std::endl;
		outfile << "3" << std::endl;
		outfile << "5" << std::endl;
	}
	outfile.close();

	return r;
}

bool TreeView::createBacFile(const std::string& meshBacFile)
{
	bool r = true;

	std::ofstream outfile(meshBacFile, std::ofstream::out);
	r &= outfile.is_open();

	if (r)
	{
		outfile << "title" << std::endl;
		outfile << "   4    2" << std::endl;
		outfile << "   1 -0.600081E+02 -0.600000E+02  0.100000E+01  0.100000E+01  0.100000E+01" << std::endl;
		outfile << "      0.100000E+01  2.500000E+01  2.500000E+01  0.100000E+01  0.100000E+01" << std::endl;
		outfile << "   2  0.600081E+02 -0.600000E+02  0.100000E+01  0.100000E+01  0.100000E+01" << std::endl;
		outfile << "      0.100000E+01  2.500000E+01  2.500000E+01  0.100000E+01  0.100000E+01" << std::endl;
		outfile << "   3  0.600081E+02  0.600000E+02  0.100000E+01  0.100000E+01  0.100000E+01" << std::endl;
		outfile << "      0.100000E+01  2.500000E+01  2.500000E+01  0.100000E+01  0.100000E+01" << std::endl;
		outfile << "   4 -0.600081E+02  0.600000E+02  0.100000E+01  0.100000E+01  0.100000E+01" << std::endl;
		outfile << "      0.100000E+01  2.500000E+01  2.500000E+01  0.100000E+01  0.100000E+01" << std::endl;
		outfile << "  1    1   2   3" << std::endl;
		outfile << "  2    1   3   4" << std::endl;
		outfile << "  sources" << std::endl;
		outfile << "  0  1	 0" << std::endl;
		outfile << "  point" << std::endl;
		outfile << "  line" << std::endl;
		outfile << "InnerRadius" << std::endl;
		switch (mData.meshDensity())
		{
			case Enum::Mesh::COURSE :
                outfile << "    0.0  0.0   0.15  0.4  0.6" << std::endl;//<< course
                outfile << "    1.0  0.0   0.15  0.4  0.6" << std::endl;//<< course
				break;
			case Enum::Mesh::MEDIUM :
                outfile << "    0.0  0.0   0.07  0.4  0.6" << std::endl;//<< medium
                outfile << "    1.0  0.0   0.07  0.4  0.6" << std::endl;//<< medium
				break;
			case Enum::Mesh::FINE :
                outfile << "    0.0  0.0   0.02  0.3  0.7" << std::endl;//<< fine
                outfile << "    1.0  0.0   0.02  0.3  0.7" << std::endl;//<< fine
				break;
			default :
                outfile << "    0.0  0.0   0.15  0.4  0.6" << std::endl;//<< course
                outfile << "    1.0  0.0   0.15  0.4  0.6" << std::endl;//<< course
		}
		outfile << " 0 " << std::endl;
		outfile << " 0 " << std::endl;
	}
	outfile.close();

	return r;
}

bool TreeView::createGeoFile(const std::string& meshGeoFile, ProjectData& data)
{
	bool r = true;
	int noPoints = 0;
	int noSegments = 6;
	int noDomainPoints = 4;

	std::ofstream outfile(meshGeoFile, std::ofstream::out);
	r &= outfile.is_open();

	noPoints = data.getProfile().size();

	if (r)
	{
		outfile << "npoin nseg nvseg nlayer hmin" << std::endl;
		outfile << int(noPoints + noDomainPoints) << "	  " << noSegments <<  "	  0	   0	    0.000" << std::endl;

		outfile.precision(7);
//		outfile << std::scientific;
		outfile << std::fixed;

		int i = 1;
		for (auto &p : data.getProfile())
		{
			outfile << i << "	  " << p.first << "	   " << p.second << std::endl;
			++i;
		}

		outfile << i << " 	 -50.000000	 -50.000000" << std::endl;
		++i;
		outfile << i << " 	  50.000000	 -50.000000" << std::endl;
		++i;
		outfile << i << " 	  50.000000	  50.000000" << std::endl;
		++i;
		outfile << i << " 	 -50.000000	  50.000000" << std::endl;


		int seg = 1;
		int prof1 = std::round(noPoints*0.5);
		int prof2 = noPoints - prof1 + 2;

		outfile << seg << " " << prof1 << " " << 1 << std::endl;
		for (i = 1; i <= prof1; ++i)
		{
			outfile << i << std::endl;
		}
		++seg;


		outfile << seg << " " << prof2 << " " << 1 << std::endl;
		for (i = prof1; i <= noPoints; ++i)
		{
			outfile << i << std::endl;
		}
		outfile << 1 << std::endl;
		++seg;


		outfile << seg << " " << 2 << " " << 3 << std::endl;
		outfile << i << " " << i+1 << std::endl;
		++seg;
		++i;

		outfile << seg << " " << 2 << " " << 3 << std::endl;
		outfile << i << " " << i+1 << std::endl;
		++seg;
		++i;

		outfile << seg << " " << 2 << " " << 3 << std::endl;
		outfile << i << " " << i+1 << std::endl;
		++seg;
		++i;

		outfile << seg << " " << 2 << " " << 3 << std::endl;
		outfile << i << " " << i-3 << std::endl;
	}
	outfile.close();

	return r;
}

bool TreeView::loadMeshProfile(const uint genNo, const std::string &filePath, ProjectData &data)
{
	bool r = true;
	int type = 1;
	std::string word1 = "";

	std::ifstream infile(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		while (infile >> word1)
		{
			if (word1 == "clean")
			{
				type = 2;
				break;
			}
		}

		infile.close();

		if (type == 1)
		{
			r &= loadMeshProfileType1(genNo, filePath, data);
		}
		else if (type == 2)
		{
			r &= loadMeshProfileType2(genNo, filePath, data);
		}
	}
	else
	{
		infile.close();
	}

	return r;
}

bool TreeView::loadMeshProfileType1(const uint genNo, const std::string& filePath, ProjectData& data)
{
	bool r = true;

	std::list<int> iBounds;
	std::list<std::pair<uint,uint>> bConnectivity;
	std::list<std::pair<float,float>> pPoints;

	//Read boundary indicies and connectivity
	std::ifstream infile(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		std::string word1("");
		std::string word2("");
		while (infile >> word1)
		{
			if ( word1 == "boundary" )
			{
				infile >> word2;
				break;
			}
		}

		if ( word1 == "boundary" &&  word2 == "sides")
		{
			//Do >> get boundary point indicies
			int i,j,b,t1,t2;
			while ( infile >> i >> j >> t1 >> b >> t2)
			{
				if (b == 1)
				{
					iBounds.push_back(i);
					bConnectivity.emplace_back(i,j);
				}
			}
			iBounds.sort();
		}
	}
	infile.close();

	//Read boundary points from indecies
	infile.open(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		std::string word3("");
		while (infile >> word3)
		{
			if ( word3 == "coordinates" )
			{
				break;
			}
		}

		if ( word3 == "coordinates" )
		{
			//Do >> get coordinates
			int i;
			float x,y;
			float t1, t2;

			while (infile >> i >> x >> y >> t1 >> t2)
			{
				if (i == iBounds.front())
				{
					pPoints.emplace_back(x,y);
					iBounds.pop_front();
				}
			}

			for (auto &pair : pPoints)
			{
				data.addBoundaryPoint(genNo, pair.first, pair.second);
			}

			for (auto& pair : bConnectivity)
			{
				data.addBConnectivity(genNo, pair.first, pair.second);
			}

		}
	}
	infile.close();

	return r;
}

bool TreeView::loadMeshProfileType2(const uint genNo, const std::string& filePath, ProjectData& data)
{
	bool r = true;

	std::list<int> iBounds;
	std::list<std::pair<uint,uint>> bConnectivity;
	std::list<std::pair<float,float>> pPoints;

	//Read boundary indicies and connectivity
	std::ifstream infile(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		std::string word1("");
		std::string word2("");
		while (infile >> word1)
		{
			if ( word1 == "Boundary" )
			{
				infile >> word2;
				break;
			}
		}

		if ( word1 == "Boundary" &&  word2 == "Faces")
		{
			//Do >> get boundary point indicies
			int i,j,b;
			while ( infile >> i >> j >> b )
			{
				if (b == 1)
				{
					iBounds.push_back(i);
					bConnectivity.emplace_back(i,j);
				}
			}
			iBounds.sort();
		}
	}
	infile.close();

	//Read boundary points from indecies
	infile.open(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		std::string word3("");
		while (infile >> word3)
		{
			if ( word3 == "Coordinates" )
			{
				break;
			}
		}

		if ( word3 == "Coordinates" )
		{
			//Do >> get coordinates
			int i;
			float x,y;

			while (infile >> i >> x >> y)
			{
				if (i == iBounds.front())
				{
					pPoints.emplace_back(x,y);
					iBounds.pop_front();
				}
			}

			for (auto &pair : pPoints)
			{
				data.addBoundaryPoint(genNo, pair.first, pair.second);
			}

			for (auto& pair : bConnectivity)
			{
				data.addBConnectivity(genNo, pair.first, pair.second);
			}
		}
	}
	infile.close();

	return r;
}

bool TreeView::loadMesh(const std::string& filePath, ProjectData& data)
{
	bool r = true;
	int type = 1;
	std::string word1 = "";

	std::ifstream infile(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		while (infile >> word1)
		{
			if (word1 == "clean")
			{
				type = 2;
				break;
			}
		}

		infile.close();

		if (type == 1)
		{
			r &= loadMeshType1(filePath, data);
		}
		else if (type == 2)
		{
			r &= loadMeshType2(filePath, data);
		}
	}
	else
	{
		infile.close();
	}

	return r;
}

bool TreeView::loadMeshType1(const std::string& filePath, ProjectData& data)
{
	bool r = true;

	//Read mesh connectivity
	std::vector<std::tuple<uint,uint,uint>> mConnectivity;
	std::ifstream infile(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		std::string word1("");
		std::string word2("");
		while (infile >> word1)
		{
			if ( word1 == "connectivities" ) break;
		}

		if ( word1 == "connectivities")
		{
			//Do >> get triangle connectivities
			uint i,j,k,t1;
			while ( infile >> word2 >> i >> j >> k >> t1)
			{
				if (word2 != "coordinates")
				{
					mConnectivity.emplace_back(i,j,k);
				}
				else break;
			}

			data.clearMeshConnectivities();
			for (auto& tuple : mConnectivity)
			{
				data.addMeshConnectivity(tuple);
			}
		}
	}
	infile.close();

	//Read mesh points
	std::vector<std::pair<float,float>> mPoints;
	infile.open(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		std::string word1("");
		std::string word2("");
		while (infile >> word1)
		{
			if ( word1 == "coordinates" ) break;
		}

		if ( word1 == "coordinates" )
		{
			float x,y,t1,t2;
			while ( infile >> word2 >> x >> y >> t1 >> t2)
			{
				if (word2 != "unknown")
				{
					mPoints.emplace_back(x,y);
				}
				else break;
			}

			data.clearMeshPoints();
			for (auto &pair : mPoints)
			{
				data.addMeshPoint(pair);
			}
		}
	}
	infile.close();

	return r;
}

bool TreeView::loadMeshType2(const std::string& filePath, ProjectData& data)
{
	bool r = true;

	//Read mesh connectivity
	std::vector<std::tuple<uint,uint,uint>> mConnectivity;
	std::ifstream infile(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		std::string word1("");
		std::string word2("");
		while (infile >> word1)
		{
			if ( word1 == "Connectivities" ) break;
		}

		if ( word1 == "Connectivities")
		{
			//Do >> get triangle connectivities
			uint i,j,k;
			while ( infile >> word2 >> i >> j >> k)
			{
				if (word2 != "Coordinates")
				{
					mConnectivity.emplace_back(i,j,k);
				}
				else break;
			}

			data.clearMeshConnectivities();
			for (auto& tuple : mConnectivity)
			{
				data.addMeshConnectivity(tuple);
			}
		}
	}
	infile.close();

	//Read mesh points
	std::vector<std::pair<float,float>> mPoints;
	infile.open(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		std::string word1("");
		std::string word2("");
		while (infile >> word1)
		{
			if ( word1 == "Coordinates" ) break;
		}

		if ( word1 == "Coordinates" )
		{
			float x,y;
			while ( infile >> word2 >> x >> y)
			{
				if (word2 != "Boundary")
				{
					mPoints.emplace_back(x,y);
				}
				else break;
			}

			data.clearMeshPoints();
			for (auto &pair : mPoints)
			{
				data.addMeshPoint(pair);
			}
		}
	}
	infile.close();

	return r;
}

bool TreeView::loadResults(const std::string& filePath, ProjectData& data)
{
	bool r = true;

	float i,rho,u,v,e,p;
	std::ifstream infile(filePath, std::ifstream::in);
	r &= infile.is_open();
	if (r)
	{
		data.clearMeshData();
		while (infile >> i >> rho >> u >> v >> e >> p)
		{
			data.addMeshData( std::make_tuple(rho, u, v, e, p) );
		}
	}
	infile.close();

	return r;
}

bool TreeView::createAerOptInFile(const std::string& filePath, ProjectData& data)
{
	bool r = true;

	std::ofstream outfile(filePath, std::ofstream::out);
	r &= outfile.is_open();

	auto& cpoints = data.getControlPoints();//<< list of uints identifying
	r &= cpoints.size() > 0;

	if (r)
	{
		std::string xrange;
		std::string mxrange;
		std::string yrange;
		std::string myrange;
        std::string smoothing_str;

        for (auto& i : cpoints)
		{
			qreal x1,y1,x2,y2;
            BoundaryPoint *point = data.getControlPoint(0, i);

            point->getBoundCoords(&x1,&y1,&x2,&y2);//<< the control point and range

			xrange += " ";
			xrange += std::to_string(x2);

			yrange += " ";
			yrange += std::to_string(y2);

			mxrange += " ";
			mxrange += std::to_string(x1);

			myrange += " ";
			myrange += std::to_string(y1);

            smoothing_str += " ";
            smoothing_str += std::to_string(point->getSmoothFactor());
		}

		outfile << "&inputVariables" << std::endl;
		outfile << "IV%Ma = " << std::to_string( data.machNo() ) << std::endl;
		outfile << "IV%Tamb = " << std::to_string( data.freeTemp() ) << std::endl;// < deg kelvin
		outfile << "IV%Pamb = " << std::to_string( data.freePress() ) << std::endl;// < Pa
		outfile << "IV%R = 287" << std::endl;
		outfile << "IV%gamma = 1.4" << std::endl;
		outfile << "IV%Re = " << std::to_string( data.reNo() ) << std::endl;
		outfile << "IV%xrange =" << xrange << mxrange << std::endl;
		outfile << "IV%yrange =" << yrange << myrange << std::endl;
		outfile << "IV%zrange = 0.00" << std::endl;
        outfile << "IV%smoothfactor =" << smoothing_str << std::endl;
		outfile << "IV%angle = 0.0" << std::endl;
		outfile << "IV%Cnconnecttrans = 0" << std::endl;
		outfile << "IV%CNconnectangle = 0" << std::endl;
        outfile << "IV%Optimiser = " << std::to_string( data.getOptimisationMethod() + 1 ) << std::endl;
		outfile << "IV%engFMF = 1.0" << std::endl;
		outfile << "IV%AlphaInflowDirection = " << std::to_string( data.freeAlpha() ) << std::endl;// < angle of attack
		outfile << "IV%turbulencemodel = 0" << std::endl;
        outfile << "IV%Top2Low = " << std::to_string(float(data.getNoTop())/100) << std::endl;
		outfile << "IV%NoSnap = " << std::to_string( data.noAgents() ) << std::endl;
		outfile << "IV%NoCN = " << cpoints.size() << std::endl;// < number of control nodes
		outfile << "IV%NoDim = 2" << std::endl;
		outfile << "IV%DoF = 8" << std::endl;// < Degrees freedom
		outfile << "IV%NoG = " << std::to_string( data.noGens() ) << std::endl;// < Generations
		outfile << "IV%NoPOMod = -1" << std::endl;
		outfile << "IV%NoLeviSteps = 100" << std::endl;
		outfile << "IV%NoIter = -3" << std::endl;
		outfile << "IV%delay = 30" << std::endl;
		outfile << "IV%constrain = .True." << std::endl;
		outfile << "IV%AdaptSamp = .FALSE." << std::endl;
		outfile << "IV%waitMax = 48" << std::endl;
		outfile << "IV%maxit = 50000" << std::endl;
		outfile << "IV%Aconst = 0.01" << std::endl;
		outfile << "IV%POD = .false." << std::endl;
		outfile << "IV%NoDelBP = 4" << std::endl;
		outfile << "IV%ObjectiveFunction = " << std::to_string( data.objFunc() ) << std::endl;
		//1 - Lift/Drag
		//2 - Distortion << no point in using this.
		//3 - max Lift
		//4 - min Drag
		//5 - max Downforce
		//6 - min Lift

		outfile << "! Test Parameters for Mesh Deformation" << std::endl;
		outfile << "IV%MeshMovement = 4" << std::endl;
		//1 - 'Linear with Smoothing' << to be tested currently unavailable.
		//2 - 'FDGD with Smoothing' << in validation stage could be good.
		//3 - 'RBF' << produces poor results thus don't use.
		//4 - 'FDGD' << use this for the moment << this is prevous version 1.
		outfile << "IV%Meshtest = .false." << std::endl;

		outfile << "! TEST POD" << std::endl;
		outfile << "IV%multiquadric = .true." << std::endl;
		outfile << "IV%Pol = .true." << std::endl;
		outfile << "IV%meanP = .true." << std::endl;

		outfile << "! Strings" << std::endl;
		outfile << "IV%filename = 'Geometry'" << std::endl;
		outfile << "IV%runOnCluster = 'N'" << std::endl;

#ifdef Q_OS_UNIX
		// do fancy unix stuff
		outfile << "IV%SystemType = 'L'" << std::endl;
#endif
#ifdef Q_OS_WIN32
		// do windows stuff here
		outfile << "IV%SystemType = 'W'" << std::endl;
#endif

		outfile << "! Login Information" << std::endl;
		outfile << "IV%UserName = 'egnaumann'" << std::endl;
		outfile << "IV%Password = 'Fleur666'" << std::endl;
		outfile << "IV%version = '2.3'" << std::endl;
		outfile << "/" << std::endl;
	}
	else
	{
		qWarning() << "Number of control points is zero or AerOpt Input File could not be created.";
	}
	outfile.close();

	return r;
}

bool TreeView::createAerOptNodeFile(const std::string& filePath, ProjectData& data)
{
	bool r = true;

	std::ofstream outfile(filePath, std::ofstream::out);
	r &= outfile.is_open();

	if (r)
	{
		std::string xrange;
		std::string mxrange;
		std::string yrange;
		std::string myrange;

		auto& cpoints = data.getControlPoints();
		for (auto& i : cpoints)
		{
            outfile << std::to_string(data.getControlPoint(0, i)->x())
					<< "	"
                    << std::to_string(data.getControlPoint(0, i)->y())
					<< std::endl;
		}
	}
	outfile.close();

	return r;
}

bool TreeView::emptyFolder(const QString& path)
{
	bool r = true;

	QDir dir(path);
	dir.setNameFilters(QStringList() << "*.*");
	dir.setFilter(QDir::Files);
	foreach(QString dirFile, dir.entryList())
	{
		dir.remove(dirFile);
	}

	return r;
}

bool TreeView::copyFolder(const QString& source, const QString& dest)
{
	bool r = true;

	//Set source folder
	QDir sourcePath(source);

	//Set/create destination folders
	QDir destPath(dest);
	destPath.mkpath(dest);

	//Copy Input folders
	QStringList filesListSource = sourcePath.entryList(QDir::Files);
	QStringList filesListDest = destPath.entryList(QDir::Files);

	//Delete destination if exists
	foreach (QString filename, filesListDest)
	{
		QFileInfo f(filename);
		QString destname = f.fileName();
		QDir dest = QDir(destPath.absolutePath() + QDir::separator() + destname).absolutePath();
		destname = QDir::toNativeSeparators(dest.path());

		if (QFile::exists(destname))
		{
			QFile::remove(destname);
		}
	}

	//Copy files from source to destination
	foreach (QString filename, filesListSource)
	{
		QFileInfo f(filename);

		QString destname = f.fileName();
		QString sourcename = f.fileName();

		QDir dest = QDir(destPath.absolutePath() + QDir::separator() + destname).absolutePath();
		destname = QDir::toNativeSeparators(dest.path());

		QDir source = QDir(sourcePath.absolutePath() + QDir::separator() + sourcename).absolutePath();
		sourcename = QDir::toNativeSeparators(source.path());

		r &= QFile::copy( sourcename, destname );
	}

	return r;
}

bool TreeView::copyFile(const QString& source, const QString& dest)
{
	bool r = true;

	if (QFile::exists(dest))
	{
		QFile::remove(dest);
	}

	r &= QFile::copy(source, dest);

	return r;
}

bool TreeView::removeFolder(const QString& path)
{
	bool r = true;
	QDir dir(path);

	if (dir.exists(path))
	{
		Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst))
		{
			if (info.isDir())
			{
				r = removeFolder(info.absoluteFilePath());
			}
			else
			{
				r = QFile::remove(info.absoluteFilePath());
			}

			if (!r)
			{
				return r;
			}
		}
		r = dir.rmdir(path);
	}
	return r;
}

bool TreeView::saveCurrentProfile(const QString& path, const ProjectData& data)
{
	bool r = true;

	std::ofstream outfile(path.toStdString(), std::ofstream::out);
	r &= outfile.is_open();

	auto& currbconn = data.getBConnects().back();
	const uint size = currbconn.size();

	auto& currbpoints = data.getBoundary().back();

	if (r)
	{
		uint i = 1;
		double x, y;
		outfile.precision(15);
		outfile << std::fixed;

		for (uint j = 0; j < size; ++j)
		{
            for (auto& point : currbconn)
			{
                if (point.first == i)
				{
                    x = currbpoints.at(point.first-1)->x();
                    y = currbpoints.at(point.first-1)->y();
					outfile << x
							<< "     "
							<< y
							<< std::endl;

                    i = point.second;
					break;
				}
			}
		}
	}
	outfile.close();

	return r;
}
