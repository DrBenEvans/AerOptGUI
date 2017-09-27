/*********************************************
**
**	Created on: 	09/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		OptimisationRun.cpp
**
**********************************************/

#include <QDebug>
#include <QPointF>
#include <QSettings>
#include <QTimer>
#include <QProcess>
#include <iostream>
#include <fstream>
#include "BoundaryPoint.h"
#include "Optimisation.h"
#include "OptimisationModel.h"
#include "FileManipulation.h"

Optimisation::Optimisation() :
    mInitMesh(new Mesh()),
    mOptimisationModel(nullptr)
{
    mObjFunc = Enum::ObjFunc::LIFTDRAG;

    mMachNo = 0.5;
    mReNo = 0.0;
    mFreeAlpha = 0.0;
    mFreePress = 101325;
    mFreeTemp = 303.15;

    mOptimisationMethod = Enum::OptMethod::MCS;
    mNoAgents = 4;
    mNoGens = 3;
    mNoTop = 75;

    mCurrentGen = 0;
    mProcess = new ProcessManager();

    auto finished = QOverload<int, QProcess::ExitStatus>::of(&ProcessManager::finished);
    auto finishedLambda = [this](int exitCode, QProcess::ExitStatus exitStatus) { optimiserFinished(exitCode, exitStatus); };
    mProcess->connect(mProcess, finished, finishedLambda);

    auto dirReadLambda = [this](const QString& path) {
        readDirectory(path);
        if(mOptimisationModel != 0) {
            mOptimisationModel->emitOptimisationDataChanged(this);
        }
    };

    auto stdOutLambda = [this](const QString line) {
        mOutputLog += line;
        if(mOptimisationModel != 0) {
            mOptimisationModel->emitOptimisationOutputChanged(this);
        }
    };

    auto stdErrLambda = [this](const QString line) {
        mOutputLog += "ERROR: "+line;
        if(mOptimisationModel != 0) {
            mOptimisationModel->emitOptimisationOutputChanged(this);
        }
    };

    mProcess->connect(mProcess, &ProcessManager::directoryChanged, dirReadLambda);
    mProcess->connect(mProcess, &ProcessManager::stdOut, stdOutLambda);
    mProcess->connect(mProcess, &ProcessManager::stdErr, stdErrLambda);
}

Optimisation::~Optimisation() {
    mProcess->cleanupProcess();
}

void Optimisation::setModel(OptimisationModel* model) {
    mOptimisationModel = model;
}

Enum::OptMethod Optimisation::getOptimisationMethod() const
{
    return mOptimisationMethod;
}

int Optimisation::getNoTop() const
{
    return mNoTop;
}

void Optimisation::setNoTop(int noTop)
{
    mNoTop = noTop;
}

//Getters and Setters
std::shared_ptr<Mesh> Optimisation::initMesh()
{
    return mInitMesh;
}

void Optimisation::setOptimisationMethod(Enum::OptMethod method)
{
    mOptimisationMethod = method;
}

int Optimisation::noGens() const
{
	return mNoGens;
}

void Optimisation::setNoGens(int noGens)
{
	mNoGens = noGens;
}

int Optimisation::noAgents() const
{
	return mNoAgents;
}

void Optimisation::setNoAgents(int noAgents)
{
	mNoAgents = noAgents;
}

float Optimisation::freeTemp() const
{
	return mFreeTemp;
}

void Optimisation::setFreeTemp(float freeTemp)
{
	mFreeTemp = freeTemp;
}

float Optimisation::freePress() const
{
	return mFreePress;
}

void Optimisation::setFreePress(float freePress)
{
	mFreePress = freePress;
}

float Optimisation::freeAlpha() const
{
	return mFreeAlpha;
}

void Optimisation::setFreeAlpha(float freeAlpha)
{
	mFreeAlpha = freeAlpha;
}

float Optimisation::reNo() const
{
	return mReNo;
}

void Optimisation::setReNo(float reNo)
{
	mReNo = reNo;
}

float Optimisation::machNo() const
{
	return mMachNo;
}

void Optimisation::setMachNo(float machNo)
{
	mMachNo = machNo;
}

Enum::ObjFunc Optimisation::objFunc() const
{
	return mObjFunc;
}

void Optimisation::setObjFunc(const Enum::ObjFunc& objFunc)
{
	mObjFunc = objFunc;
}

QString Optimisation::label() const {
    return mLabel;
}

void Optimisation::setLabel(QString label) {
    mLabel = label;
}

QString Optimisation::simulationDirectoryName() {
    QString fileName = "run_" + mLabel;
    fileName.replace(QRegExp(QString::fromUtf8("[^a-zA-Z\\d]")),"_");
    return fileName;
}

std::vector<BoundaryPoint*> Optimisation::controlPoints() {
    return mControlPoints;
}

void Optimisation::setControlPoints(std::vector<BoundaryPoint*> controlPoints) {
    mControlPoints = controlPoints;
}

int Optimisation::controlPointCount() {
    return mControlPoints.size();
}

QString Optimisation::outputDataDirectory() {
    QSettings settings;
    QString outputData = settings.value("AerOpt/workingDirectory").toString();
    outputData += "/" + simulationDirectoryName() + "/Output_Data";
    outputData = QDir::toNativeSeparators(outputData);
    return outputData;
}

bool Optimisation::run() {
    // load settings
    QSettings settings;
    QString AerOptInFile = settings.value("AerOpt/inputFile").toString();
    QString AerOptNodeFile = settings.value("AerOpt/nodeFile").toString();
    QString meshDatFile = settings.value("mesher/initMeshFile").toString();
    QString AerOpt = settings.value("AerOpt/executable").toString();
    QString AerOptWorkDir = settings.value("AerOpt/workingDirectory").toString();
    QDir scratchDir = QDir(settings.value("mesher/scratchDir").toString());
    QString inFolder = settings.value("AerOpt/inFolder").toString();

    bool r = true;

    //Empty the input directory
    r &= FileManipulation::emptyFolder(inFolder);

    //Copy mesh file as input to AerOpt from project directory to input directory
    QString dest = QDir::toNativeSeparators(inFolder + "/Mesh.dat");
    r &= FileManipulation::copyFile(meshDatFile, dest);

    //Set input files in input directory
    r &= createAerOptInFile(AerOptInFile);
    r &= createAerOptNodeFile(AerOptNodeFile);

    //then run aeropt
    if (r)
    {
        mProcess->run(AerOpt, AerOptWorkDir, outputDataDirectory());
    }
    else
    {
        qWarning() << "Process Aborted!";
    }
    return r;
}

bool Optimisation::createAerOptInFile(const QString& filePath)
{
    bool r = true;
    QSettings settings;
    QString workingDirectory = settings.value("AerOpt/workingDirectory").toString();

    std::ofstream outfile(filePath.toStdString(), std::ofstream::out);
    r &= outfile.is_open();

    int ctlPointCount = controlPointCount();
    r &= ctlPointCount > 0;

    if (r)
    {
        std::string xrange;
        std::string mxrange;
        std::string yrange;
        std::string myrange;
        std::string smoothing;

        for (auto& point : controlPoints())
        {
            QRectF rect = point->controlPointRect();
            float smoothFactor = point->getSmoothFactor();

            xrange += " " + std::to_string(rect.right());

            yrange += " " + std::to_string(rect.bottom());

            mxrange += " " + std::to_string(rect.left());

            myrange += " " + std::to_string(rect.top());

            smoothing += " " + std::to_string(smoothFactor);
        }

        Enum::OptMethod optMethod = getOptimisationMethod();
        uint optMethodIndex;
        switch(optMethod) {
            case Enum::OptMethod::MCS: optMethodIndex = 1;
            case Enum::OptMethod::DE: optMethodIndex = 2;
            case Enum::OptMethod::PSO: optMethodIndex = 3;
            default: 1;
        }

        outfile << "&inputVariables" << std::endl;
        outfile << "IV%Ma = " << std::to_string( machNo() ) << std::endl;
        outfile << "IV%Tamb = " << std::to_string( freeTemp() ) << std::endl;// < deg kelvin
        outfile << "IV%Pamb = " << std::to_string( freePress() ) << std::endl;// < Pa
        outfile << "IV%R = 287" << std::endl;
        outfile << "IV%gamma = 1.4" << std::endl;
        outfile << "IV%Re = " << std::to_string( reNo() ) << std::endl;
        outfile << "IV%xrange =" << xrange << mxrange << std::endl;
        outfile << "IV%yrange =" << yrange << myrange << std::endl;
        outfile << "IV%zrange = 0.00" << std::endl;
        outfile << "IV%smoothfactor =" << smoothing << std::endl;
        outfile << "IV%smoothconvergence = -3" << std::endl;
        outfile << "IV%angle = 0.0" << std::endl;
        outfile << "IV%Cnconnecttrans = 0" << std::endl;
        outfile << "IV%engFMF = 1.0" << std::endl;
        outfile << "IV%AlphaInflowDirection = " << std::to_string( freeAlpha() ) << std::endl;// < angle of attack
        outfile << "IV%YawInflowAngle = 0.0" << std::endl;
        outfile << "IV%turbulencemodel = 0" << std::endl;
        outfile << "IV%Low2Top = " << std::to_string(float(getNoTop())/100) << std::endl;
        outfile << "IV%NoSnap = " << std::to_string( noAgents() ) << std::endl;
        outfile << "IV%NoNests = " << std::to_string( noAgents() ) << std::endl;
        outfile << "IV%NoCN = " << ctlPointCount << std::endl;// < number of control nodes
        outfile << "IV%NoDim = 2" << std::endl;
        outfile << "IV%DoF = 8" << std::endl;// < Degrees freedom
        outfile << "IV%NoG = " << std::to_string( noGens() ) << std::endl;// < Generations
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
        outfile << "IV%ObjectiveFunction = " << std::to_string( objFunc() ) << std::endl;
        //1 - Lift/Drag
        //2 - Distortion << no point in using this.
        //3 - max Lift
        //4 - min Drag
        //5 - max Downforce
        //6 - min Lift
        outfile << "IV%Optimiser = " << std::to_string( optMethodIndex ) << std::endl;

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
        outfile << "IV%filepath = '" << workingDirectory.toStdString() << "'" << std::endl;
        outfile << "IV%SimulationName = '" << simulationDirectoryName().toStdString() << "'" << std::endl;
        outfile << "IV%filename = 'Geometry'" << std::endl;
        outfile << "IV%Meshfilename = 'Mesh'" << std::endl;
        outfile << "IV%runOnCluster = 'N'" << std::endl;

#ifdef Q_OS_UNIX
        // do fancy unix stuff
        outfile << "IV%SystemType = 'L'" << std::endl;
#endif
#ifdef Q_OS_WIN32
        // do windows stuff here
        outfile << "IV%SystemType = 'W'" << std::endl;
#endif
        outfile << "IV%NoProcessors = 1" << std::endl;
        outfile << "IV%shapeenclosed = .true." << std::endl;

        outfile << "! Login Information" << std::endl;
        outfile << "IV%UserName = 'egnaumann'" << std::endl;
        outfile << "IV%Password = 'Fleur666'" << std::endl;
        outfile << "IV%version = '3.5'" << std::endl;
        outfile << "/" << std::endl;
    }
    else
    {
        qWarning() << "Number of control points is zero or AerOpt Input File could not be created.";
    }
    outfile.close();

    return r;
}

bool Optimisation::createAerOptNodeFile(const QString& filePath)
{
    bool r = true;

    std::ofstream outfile(filePath.toStdString(), std::ofstream::out);
    r &= outfile.is_open();

    if (r)
    {
        std::string xrange;
        std::string mxrange;
        std::string yrange;
        std::string myrange;

        for (auto& point : controlPoints())
        {
            outfile << std::to_string(point->x())
                    << "	"
                    << std::to_string(point->y())
                    << std::endl;
        }
    }
    outfile.close();

    return r;
}

void Optimisation::optimiserFinished(int exitCode, QProcess::ExitStatus exitStatus)
{
    QSettings settings;

    //Windows bug!! re read last file set if exists.
    QEventLoop loop;
    QTimer::singleShot(1000, &loop, SLOT(quit()));
    loop.exec();

    QDir lastPath(settings.value("AerOpt/workingDirectory").toString());
    lastPath = QDir(lastPath.absolutePath() + QDir::separator() + "AerOpt/FLITE/Output_Data").absolutePath();
    lastPath = QDir::toNativeSeparators(lastPath.path());
    readDirectory(lastPath.path());

    QEventLoop loop2;
    QTimer::singleShot(2000, &loop2, SLOT(quit()));
    loop2.exec();

    if (exitStatus == QProcess::NormalExit && exitCode == 0)
    {
        qDebug() << "Process finished normally";

        bool r = true;

        r &= mProcess->exitCode() == 0;
    }
    else
    {
        qCritical() << "Process finished abnormally with exit code: " << exitCode;
    }

    //Move files to project folder here!
    bool r = true;

    //Get source folders (input and output)
    QDir appPath(settings.value("AerOpt/workingDirectory").toString());

    //Get destination folders (input and output)
    QDir projPath(settings.value("mesher/scratchDir").toString());
    QDir projPathOut = QDir(projPath.absolutePath() + QDir::separator() + "Output_Data").absolutePath();
    QDir projPathIn = QDir(projPath.absolutePath() + QDir::separator() + "Input_Data").absolutePath();
    projPathOut = QDir::toNativeSeparators(projPathOut.path());
    projPathIn = QDir::toNativeSeparators(projPathIn.path());

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
        r &= FileManipulation::removeFolder( temp );
    }

    if (!r) qWarning() << "Something went wrong with file cleanup!";

    //TODO
    //mCanvas.update();
}

void Optimisation::readDirectory(const QString& path)
{
    //Read initial fitness
    if (mCurrentGen == 0)
    {
        bool b = true;
        QString fit = path;
        fit += "Fitness_";
        fit += QString::number(mCurrentGen);
        fit += ".txt";
        QFileInfo checkFit(fit);
        b &= checkFit.exists();
        b &= checkFit.isFile();
        if (b)
        {
            readFitness(fit);
            mCurrentGen++;
        }
    }

    //When meshfiles exists, load it.
    bool c = true;

    //Fitness file
    QString fit = path;
    fit += "Fitness_";
    fit += QString::number(mCurrentGen);
    fit += ".txt";

    //Mesh file
    QString meshpath = path;
    meshpath += "Geometry_";
    meshpath += QString::number(mCurrentGen);
    meshpath += ".dat";

    //Results file
    QString results = path;
    results += "Geometry_";
    results += QString::number(mCurrentGen);
    results += ".resp";

    //Load Mesh + Results
    QFileInfo checkFit(fit);
    QFileInfo checkMesh(meshpath);
    QFileInfo checkResults(results);

    c &= checkFit.exists();
    c &= checkFit.isFile();
    c &= checkMesh.exists();
    c &= checkMesh.isFile();
    c &= checkResults.exists();
    c &= checkResults.isFile();

    if (c)
    {
        if(mMeshes.size() < mCurrentGen) {
            mMeshes.emplace_back(new Mesh());
        }

        Mesh* mesh = mMeshes.at(mCurrentGen-1);
        mesh->loadMesh(meshpath);
        mesh->loadResults(results);

        readFitness(fit);

        mCurrentGen++;
    }
}

void Optimisation::readFitness(const QString& path)
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
        qDebug() << output;


        //Extract from plotData gen no. and then nest values
        QStringList list1 = plotData.split(" ", QString::SkipEmptyParts);
        std::vector<double> nests;
        uint genNo = 0;
        bool ok;
        for (int i = 0; i < list1.size(); ++i)
        {
            if (i == 0)
                genNo = list1.at(i).toInt(&ok);
            else
                nests.push_back( list1.at(i).toDouble(&ok) );

            if (genNo == 0 && i == 1)
            {
                int cols = noAgents() -1;
                for (int j = 0; j < cols; ++j)
                {
                    nests.push_back( nests.at(0) );
                }
            }
        }

        //read fitness
        if (ok)
        {
            genNo += 1;
            if(mFitness.size() < genNo)
                mFitness.emplace_back(nests);
            else
                mFitness.at(genNo) = nests;
        }

    }
    infile.close();
}

std::vector<std::vector<double> > Optimisation::fitness() {
    return mFitness;
}

QString Optimisation::outputText() {
    return mOutputLog;
}
