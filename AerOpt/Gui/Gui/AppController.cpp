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

#include "AppController.h"
#include "Canvas.h"
#include "Menu.h"
#include "OptimisationRun.h"
#include "Enumerations.h"
#include "MeshDialog.h"
#include "PlotterDialog.h"
#include "ConfigSimulationDialog.h"

AppController::AppController(OptimisationRun& data, Canvas& canvas, QWidget *parent) : QDialog(parent), mCanvas(canvas), mData(data)
{
    mParent = parent;
    myMeshProcess.setParent(mParent);
    myOptProcess.setParent(mParent);
    myDirWatcher.setParent(mParent);

    mProjectDirectory = "";

    //Determine application directory and binary locations.
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
    qDebug() << appPath;
	appPath = QDir::fromNativeSeparators(QCoreApplication::applicationFilePath());

	index = appPath.lastIndexOf("/AerOpt/");

	if (index != -1)
	{
		//Determin app path
		appPath.remove(index, appPath.size() - index);
		appPath = QDir::toNativeSeparators(appPath);
        qDebug() << appPath;
		mAppPath = appPath;

		//Determin mesher path
		mesherExe = QDir::toNativeSeparators(mesherExe);
		mesherExe.prepend( appPath );
        qDebug() << mesherExe;
		mMesherPath = mesherExe;

		//Determin flite path << lol
		aeroptExe = QDir::toNativeSeparators(aeroptExe);
		aeroptExe.prepend( appPath );
        qDebug() << aeroptExe;
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

    qInfo() << " *** Welcome to AerOpt ***";
    qInfo() << " ***     Have a nice day     ***";
}

AppController::~AppController()
{
	myMeshProcess.kill();
	myOptProcess.kill();
}

void AppController::clearProject()
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

void AppController::new_simulation() {
    mData = OptimisationRun();
    ConfigSimulationDialog diag(mData);
    diag.show();
    diag.exec();
}

void AppController::configure_current_simulation() {
    ConfigSimulationDialog diag(mData);
    diag.show();
    diag.exec();
}

void AppController::loadProject()
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

void AppController::runMesher()
{
	sGenNo = 0;
	//Clear mesh data first
    mData.setMesh(false);
    mData.setRenderMesh(false);
    mData.setRenderProfile(true);
	mData.clearMesh();

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
		myMeshProcess.setWorkingDirectory( meshWorkDir );
		myMeshProcess.setStandardInputFile( meshInFile );
		myMeshProcess.start( mMesherPath );
	}
}

void AppController::setBoundary()
{
	bool r = true;
    if (!mData.mesh())
	{
        qWarning() << "Please first set the mesh!";
		return;
	}
    ConfigSimulationDialog diag(mData, this);

	diag.show();
	diag.exec();

	r &= mData.boundary();
    if (r) mCurrentNode->setText(1, "Yes");
}

void AppController::setOptimiser()
{
	bool r = true;
	if (!mData.boundary())
	{
        qWarning() << "Please first set the flow conditions!";
		return;
	}
    ConfigSimulationDialog diag(mData, this);

	diag.show();
	diag.exec();

	r &= mData.optimiser();
    if (r) mCurrentNode->setText(1, "Yes");
}

void AppController::runAerOpt()
{
	bool r = true;

	if (!mData.optimiser())
	{
		qWarning() << "Please first set the optimiser parameters!";
		return;
	}

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
		myDirWatcher.addPath(path);
		myOptProcess.setWorkingDirectory( AerOptWorkDir );
		myOptProcess.start( AerOpt );
	}
	else
	{
		qWarning() << "Process Aborted!";
	}
}

void AppController::stopMesher()
{
	myMeshProcess.kill();
    qInfo() << "Any running mesh jobs have been stopped!";
}

void AppController::stopAerOpt()
{
	myOptProcess.kill();
    qInfo() << "Any running AerOpt jobs have been stopped!";
}

void AppController::showGraph()
{
	mPlotter->show();
	mPlotter->adjustSize();
}

//Private slots
void AppController::meshOutput()
{
	QByteArray byteArray = myMeshProcess.readAllStandardOutput();
	QStringList strLines = QString(byteArray).split("\n");

	foreach (QString line, strLines){
        qInfo() << line;
	}
}

void AppController::meshError()
{
	QByteArray byteArray = myMeshProcess.readAllStandardError();
	QStringList strLines = QString(byteArray).split("\n");

	foreach (QString line, strLines){
		qCritical() << line;
	}
}

void AppController::processOutput()
{
	QByteArray byteArray = myOptProcess.readAllStandardOutput();
	QStringList strLines = QString(byteArray).split("\n");

	foreach (QString line, strLines){
        qInfo() << line;
	}
}

void AppController::processError()
{
	QByteArray byteArray = myOptProcess.readAllStandardError();
	QStringList strLines = QString(byteArray).split("\n");

	foreach (QString line, strLines){
		qCritical() << line;
	}
}

void AppController::meshingStarted()
{
    qInfo() << "Meshing process started normally";
}

void AppController::meshingFinished(int exitCode, QProcess::ExitStatus exitStatus)
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
}

void AppController::optimiserStarted()
{
    qInfo() << "Optimiser process started normally";
	mPlotter->clearData();
	mPlotter->show();
}

void AppController::readDirectory(const QString& path)
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

void AppController::readFitness(const QString& path)
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

void AppController::optimiserFinished(int exitCode, QProcess::ExitStatus exitStatus)
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
}


//Private functions

bool AppController::createInputFile(const std::string& meshInFile,
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

bool AppController::createBacFile(const std::string& meshBacFile)
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

bool AppController::createGeoFile(const std::string& meshGeoFile, OptimisationRun& data)
{
	bool r = true;
	int noPoints = 0;
	int noSegments = 6;
	int noDomainPoints = 4;
    int nViscSegments = 0;

	std::ofstream outfile(meshGeoFile, std::ofstream::out);
	r &= outfile.is_open();

	noPoints = data.getProfile().size();
    int nLayers = data.getNumBoundaryLayers();
    bool hasBoundaryLayer = data.getBoundaryLayerFlag();

    if(hasBoundaryLayer) {
        nViscSegments = 2;
    };

	if (r)
	{
		outfile << "npoin nseg nvseg nlayer hmin" << std::endl;
        outfile << int(noPoints + noDomainPoints) << "	  " << noSegments <<  "	  " << nViscSegments << "	  " << nLayers << "   0.000" << std::endl;

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

    if(hasBoundaryLayer) {
        qreal thickness = data.getBoundaryLayerThickness();
        qreal dist = thickness;
        for(int i=0; i<nLayers; i++) {
            dist += thickness;
            outfile << dist << std::endl;
        }
    }

	outfile.close();

	return r;
}

bool AppController::loadMeshProfile(const uint genNo, const std::string &filePath, OptimisationRun &data)
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

bool AppController::loadMeshProfileType1(const uint genNo, const std::string& filePath, OptimisationRun& data)
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

bool AppController::loadMeshProfileType2(const uint genNo, const std::string& filePath, OptimisationRun& data)
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

bool AppController::loadMesh(const std::string& filePath, OptimisationRun& data)
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

bool AppController::loadMeshType1(const std::string& filePath, OptimisationRun& data)
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

bool AppController::loadMeshType2(const std::string& filePath, OptimisationRun& data)
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

bool AppController::loadResults(const std::string& filePath, OptimisationRun& data)
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

bool AppController::createAerOptInFile(const std::string& filePath, OptimisationRun& data)
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

bool AppController::createAerOptNodeFile(const std::string& filePath, OptimisationRun& data)
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

bool AppController::emptyFolder(const QString& path)
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

bool AppController::copyFolder(const QString& source, const QString& dest)
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

bool AppController::copyFile(const QString& source, const QString& dest)
{
	bool r = true;

	if (QFile::exists(dest))
	{
		QFile::remove(dest);
	}

	r &= QFile::copy(source, dest);

	return r;
}

bool AppController::removeFolder(const QString& path)
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

bool AppController::saveCurrentProfile(const QString& path, const OptimisationRun& data)
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
