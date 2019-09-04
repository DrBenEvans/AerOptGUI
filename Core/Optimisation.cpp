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
#include <sstream>
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

    mProcess = new ProcessManager();

    auto finished = QOverload<int, QProcess::ExitStatus>::of(&ProcessManager::finished);
    auto finishedLambda = [this](int exitCode, QProcess::ExitStatus exitStatus) { optimiserFinished(exitCode, exitStatus); };
    mProcess->connect(mProcess, finished, finishedLambda);

    auto dirReadLambda = [this](const QString& path) {
        readFitness();
        if(mOptimisationModel != nullptr) {
            mOptimisationModel->emitOptimisationFitnessChanged(this);
        }
    };

    auto stdOutLambda = [this](const QString line) {
        addToOutputLog(line);
        if(mOptimisationModel != nullptr) {
            mOptimisationModel->emitOptimisationOutputChanged(this);
        }
    };

    auto stdErrLambda = [this](const QString line) {
        addToOutputLog("Standard Error: "+line);
        if(mOptimisationModel != nullptr) {
            mOptimisationModel->emitOptimisationOutputChanged(this);
        }
    };

    mProcess->connect(mProcess, &ProcessManager::directoryChanged, dirReadLambda);
    mProcess->connect(mProcess, &ProcessManager::stdOut, stdOutLambda);
    mProcess->connect(mProcess, &ProcessManager::stdErr, stdErrLambda);


    clusterChecker = new clusterManager();
    clusterChecker->connect(clusterChecker, &clusterManager::directoryChanged, dirReadLambda);
    clusterChecker->connect(clusterChecker, &clusterManager::stdOut, stdOutLambda);
    clusterChecker->connect(clusterChecker, &clusterManager::stdErr, stdErrLambda);
}

Optimisation::~Optimisation() {
    if( !runOnCluster ) {
        mProcess->cleanupProcess();
    }
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
Mesh* Optimisation::initMesh()
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

std::vector<BoundaryPoint*> Optimisation::initialBoundaryPoints() {
    return mBoundaryPoints;
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

void Optimisation::setBoundaryPoints(std::vector<BoundaryPoint*> boundaryPoints) {
    mBoundaryPoints = boundaryPoints;
}

QString Optimisation::simulationDirectoryPath() {
    QSettings settings;
    QString dir = settings.value("AerOpt/workingDirectory").toString();
    dir += "/" + simulationDirectoryName();
    dir = QDir::toNativeSeparators(dir);
    return dir;
}

QString Optimisation::outputDataDirectory() {
    QString outputData = simulationDirectoryPath();
    outputData += "/Output_Data";
    outputData = QDir::toNativeSeparators(outputData);
    return outputData;
}

ProfilePoints Optimisation::initProfilePoints() {
    return mProfilePoints;
}

void Optimisation::setInitProfilePoints(ProfilePoints profilePoints) {
    mProfilePoints = profilePoints;
}

void Optimisation::writeProfilePointsToSimulationDir() {
    QString filePath = simulationDirectoryPath() + "/profilePoints.txt";
    QDir::toNativeSeparators(filePath);

    Mesh* mesh = initMesh();
    if(mesh) {
        std::ofstream outfile(filePath.toStdString(), std::ofstream::out);
        if(outfile.is_open()) {
            float x, y;
            for(auto& profilePoint : initProfilePoints()) {
                x = profilePoint.first;
                y = profilePoint.second;

                outfile << x << " " << y << std::endl;
            }
        }

        outfile.close();
    }
}

bool Optimisation::readProfilePointsFromSimulationDir() {
    QString filePath = simulationDirectoryPath() + "/profilePoints.txt";
    QDir::toNativeSeparators(filePath);

    Profile profile;
    bool success = profile.setFile(filePath);
    mProfilePoints = profile.getProfile();

    return success;
}

void Optimisation::copyFileToSimulationDir(QString source) {
    QDir pathdir(source);
    QString dirName = pathdir.dirName();

    QString dest = simulationDirectoryPath() + "/" + dirName;
    QDir::toNativeSeparators(dest);

    FileManipulation::copyFile(source, dest);
}


bool Optimisation::run() {
    // load settings
    QSettings settings;
    QString AerOptInFile = settings.value("AerOpt/inputFile").toString();
    QString AerOptNodeFile = settings.value("AerOpt/nodeFile").toString();
    QString AerOptBoundaryFile = settings.value("AerOpt/boundaryFile").toString();
    QString meshDatFile = settings.value("mesher/initMeshFile").toString();
    QString AerOpt = settings.value("AerOpt/executable").toString();
    QString AerOptWorkDir = settings.value("AerOpt/workingDirectory").toString();
    QString inFolder = settings.value("AerOpt/inFolder").toString();
    QString username = settings.value("Cluster/Username").toString();
    QString address = settings.value("Cluster/Address").toString();

    setInitProfilePoints(initMesh()->profilePoints());
    bool r = true;

    //Empty the input directory
    r &= FileManipulation::emptyFolder(inFolder);

    //Copy mesh file as input to AerOpt from project directory to input directory
    QString dest = QDir::toNativeSeparators(inFolder + "/Mesh.dat");
    r &= FileManipulation::copyFile(meshDatFile, dest);

    //Set input files in input directory
    r &= createAerOptInFile(AerOptInFile);
    r &= createAerOptNodeFile(AerOptNodeFile);
    r &= createAerOptBoundaryPointFile(AerOptBoundaryFile);

    //then run aeropt
    if (r)
    {
        if(runOnCluster){

            QString dirname = simulationDirectoryName();

            clusterChecker->setWorkingDirectory(dirname);
            clusterChecker->setUsername(username);
            clusterChecker->setClusterAddress(address);
            clusterChecker->setPassword(mClusterPassword);
            clusterChecker->start();

            QDir dir(simulationDirectoryPath());
            dir.mkdir(".");

        } else {

            mProcess->run(AerOpt, AerOptWorkDir, outputDataDirectory());
            qInfo() << "Copying input files to optimisation output directory";

        }
        copyFileToSimulationDir(AerOptInFile);
        copyFileToSimulationDir(AerOptNodeFile);
        copyFileToSimulationDir(AerOptBoundaryFile);
        writeProfilePointsToSimulationDir();
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
    std::string username = settings.value("Cluster/Username").toString().toStdString();
    std::string account = settings.value("Cluster/Account").toString().toStdString();

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
            xrange += " " + std::to_string(point->controlPointRect().right());
            yrange += " " + std::to_string(point->controlPointRect().bottom());
            mxrange += " " + std::to_string(point->controlPointRect().left());
            myrange += " " + std::to_string(point->controlPointRect().top());
            smoothing += " " + std::to_string(point->getSmoothFactor());
        }

        Enum::OptMethod optMethod = getOptimisationMethod();
        uint optMethodIndex;
        switch(optMethod) {
        case Enum::OptMethod::MCS:
            optMethodIndex = 1;
            break;
        case Enum::OptMethod::DE:
            optMethodIndex = 2;
            break;
        case Enum::OptMethod::PSO:
            optMethodIndex = 3;
            break;
        default:
            optMethodIndex = -1;
            break;
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
        if(runOnCluster){
            std::ifstream infile("clusterlogin", std::ifstream::in);
            std::string runname = simulationDirectoryName().toStdString();
            std::string path = "/home/" + username + "/AerOpt/" + runname;
            outfile << "IV%filepath = '" << path << "'" << std::endl;
            outfile << "IV%clusterpath = '" << path << "'" << std::endl;
            outfile << "IV%OutputDirectory = '" << runname << "'" << std::endl;
            outfile << "IV%filename = 'Geometry'" << std::endl;
            outfile << "IV%Account = '" << account << "'" << std::endl;
            outfile << "IV%Meshfilename = 'Mesh'" << std::endl;
            outfile << "IV%runOnCluster = 'Y'" << std::endl;
            outfile << "IV%SystemType = 'B'" << std::endl;
        } else{
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

        }
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
        qWarning() << "Number of control nodes is zero or AerOpt Input File could not be created.";
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
                    << "\t"
                    << std::to_string(point->y())
                    << std::endl;
        }
    }
    outfile.close();

    return r;
}

bool Optimisation::createAerOptBoundaryPointFile(const QString &filePath) {

    bool r = true;

    std::ofstream outfile(filePath.toStdString(), std::ofstream::out);
    r &= outfile.is_open();

    if (r)
    {
        std::string xcoord;
        std::string ycoord;
        std::string controlPoint;
        std::string smoothingFactor;
        std::string xmin;
        std::string xmax;
        std::string ymin;
        std::string ymax;

        for (auto& point : mBoundaryPoints)
        {
            xcoord += std::to_string(point->x()) + "\t";
            ycoord += std::to_string(point->y()) + "\t";
            controlPoint += std::to_string(point->isControlPoint()) + "\t";
            smoothingFactor += std::to_string(point->getSmoothFactor()) + "\t";
            xmin += std::to_string(point->controlPointRect().left()) + "\t";
            xmax += std::to_string(point->controlPointRect().right()) + "\t";
            ymin += std::to_string(point->controlPointRect().top()) + "\t";
            ymax += std::to_string(point->controlPointRect().bottom()) + "\t";

        }

        // Boundary Point Info
        outfile << "IV%x = " << xcoord << std::endl;
        outfile << "IV%y = " << ycoord << std::endl;
        outfile << "IV%isCP = " << controlPoint << std::endl;
        outfile << "IV%smoothingFactor = " << smoothingFactor << std::endl;
        outfile << "IV%xMin = " << xmin << std::endl;
        outfile << "IV%xMax = " << xmax << std::endl;
        outfile << "IV%yMin = " << ymin << std::endl;
        outfile << "IV%yMax = " << ymax << std::endl;

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

    readFitness();
}

Mesh* Optimisation::mesh(int genIndex, int agentIndex)
{
    bool success = true;

    Mesh* mesh = nullptr;

    QString base_path = outputDataDirectory();
    base_path += "/Geometry";
    base_path += QString::number(genIndex + 1);
    base_path += "_";
    base_path += QString::number(agentIndex + 1);

    //Mesh file
    QString  meshpath = base_path + ".dat";
    //Results file
    QString results = base_path += ".unk";

    //Load Mesh + Results
    QFileInfo checkMesh(meshpath);
    QFileInfo checkResults(results);

    success &= checkMesh.exists();
    success &= checkMesh.isFile();
    success &= checkResults.exists();
    success &= checkResults.isFile();

    if (success)
    {
        mesh = new Mesh();
        success &= mesh->loadMesh(meshpath);
        success &= mesh->loadResults(results);

        if(success) {
            return mesh;
        } else{
            delete mesh;
        }
    }

    return nullptr;
}

bool Optimisation::readFitness()
{
    bool r = true;
    bool isok = true;

    std::string line = "";
    QString output = "";
    QString fitnessString = "";

    QString path = simulationDirectoryPath() + "/FitnessAll.txt";
    path = QDir::toNativeSeparators(path);

    QFileInfo checkFit(path);

    r &= checkFit.exists();
    r &= checkFit.isFile();

    if(!r) {
        return r;
    }

    std::ifstream infile(path.toStdString(), std::ifstream::in);
    r &= infile.is_open();

    if (r)
    {
        uint genNo = 0;
        std::vector<double> nests;
        bool testok;

        mFitness.clear();

        while(std::getline(infile, line)) {
            output = QString::fromStdString(line);
            output = output.trimmed();
            fitnessString = output;
            output.prepend( "Generation: > " );
            qDebug() << output;

            //Extract from fitnessString gen no. and then nest values
            QStringList fitnessList = fitnessString.split(" ", QString::SkipEmptyParts);

            //Read Generation Number from first row
            genNo = fitnessList.at(0).toInt(&testok);
            uint genIndex = genNo - 1;
            isok = isok && testok;

            //Extend fitness vector with new generation(s) if required
            while(mFitness.size() <= genIndex) {
                mFitness.emplace_back();
            }

            //Read Nest Fitnesses from subsequent rows
            uint agentIndex;
            for (int index = 1; index < fitnessList.size(); ++index) {
                // get agent fitness
                double agentFitness = fitnessList.at(index).toDouble(&testok);
                isok = isok && testok;

                // set agent fitness
                agentIndex = index - 1;
                if(mFitness.at(genIndex).size() <= agentIndex) {
                    mFitness.at(genIndex).emplace_back(agentFitness);
                } else {
                    mFitness.at(genIndex).at(agentIndex) = agentFitness;
                }
            }

        }

        infile.close();
    }

    return isok & r;
}

std::vector<std::vector<double> > Optimisation::allfitness() {
    return mFitness;
}

double Optimisation::fitness(int generationIndex, int agentIndex) {
    if(mFitness.size() > generationIndex && mFitness.at(generationIndex).size() > agentIndex) {
        return mFitness.at(generationIndex).at(agentIndex);
    } else {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

std::pair<double,double> Optimisation::fitnessRange() {
    std::pair<double,double> minMax;
    minMax.first = std::numeric_limits<double>::infinity();
    minMax.second = -std::numeric_limits<double>::infinity();

    for(auto& gen : mFitness)
        for(auto& agentFitness : gen) {
            if(agentFitness < minMax.first)
                minMax.first = agentFitness;
            if(agentFitness > minMax.second)
                minMax.second = agentFitness;

        }
    return minMax;
}

std::pair<double,double> Optimisation::fitnessRangeBestAgents() {
    std::pair<double,double> minMax;
    minMax.first = std::numeric_limits<double>::infinity();
    minMax.second = -std::numeric_limits<double>::infinity();

    for(auto& gen : mFitness){
        if(gen.at(0) < minMax.first)
            minMax.first = gen.at(0);
        if(gen.at(0) > minMax.second)
            minMax.second = gen.at(0);
    }
    return minMax;
}

QString Optimisation::outputText() {
    return mOutputLog;
}

bool Optimisation::readAerOptSettings(QString filePath) {
    bool success = true;
    std::ifstream infile(filePath.toStdString(), std::ifstream::in);
    success &= infile.is_open();

    QString qline, value;
    QStringList strList;
    std::string variable;

    if (success)
    {
        std::string line("");

        while (std::getline(infile, line))
        {
            qline = QString::fromStdString(line);
            strList = qline.split(QRegExp("\\s*=\\s*"));

            if(strList.length()==2) {
                variable = strList.at(0).toStdString();
                value = strList.at(1);

                if(variable == "IV%Ma") {
                    setMachNo(value.toFloat());
                }
                else if(variable == "IV%Tamb") {
                    setFreeTemp(value.toFloat());
                }
                else if(variable == "IV%Pamb") {
                    setFreePress(value.toFloat());
                }
                else if(variable == "IV%Re") {
                    setReNo(value.toFloat());
                }
                else if(variable == "IV%AlphaInflowDirection") {
                    setFreeAlpha(value.toFloat());
                }
                else if(variable == "IV%Low2Top") {
                    setNoTop(int(value.toFloat()*100));
                }
                else if(variable == "IV%NoSnap") {
                    setNoAgents(value.toInt());
                }
                else if(variable == "IV%NoNests") {
                    setNoAgents(value.toInt());
                }
                else if(variable == "IV%NoG") {
                    setNoGens(value.toInt());
                }
                else if(variable == "IV%ObjectiveFunction") {
                    int ivalue = value.toInt();
                    if(ivalue == 1)
                        setObjFunc(Enum::ObjFunc::LIFTDRAG);
                    else if(ivalue == 2)
                        setObjFunc(Enum::ObjFunc::DISTORTION);
                    else if(ivalue == 3)
                        setObjFunc(Enum::ObjFunc::MAXLIFT);
                    else if(ivalue == 4)
                        setObjFunc(Enum::ObjFunc::MINDRAG);
                    else if(ivalue == 5)
                        setObjFunc(Enum::ObjFunc::MAXDOWNFORCE);
                    else if(ivalue == 6)
                        setObjFunc(Enum::ObjFunc::MINLIFT);
                    else
                        setObjFunc(Enum::ObjFunc::FUNCNOTSET);
                }
                else if(variable == "IV%Optimiser") {
                    int ivalue = value.toInt();
                    if(ivalue == 1)
                        setOptimisationMethod(Enum::OptMethod::MCS);
                    else if(ivalue == 2)
                        setOptimisationMethod(Enum::OptMethod::DE);
                    else if(ivalue == 3)
                        setOptimisationMethod(Enum::OptMethod::PSO);
                    else
                        setOptimisationMethod(Enum::OptMethod::METHODNOTSET);
                }
                else if(variable == "IV%SimulationName") {
                    if(strList.size() > 0) {
                        strList = strList.at(1).split("run_");
                        QString label =  strList.at(1);
                        // remove trailing quotation mark
                        label.chop(1);

                        setLabel(label);
                    }
                }
                else if(variable == "IV%xrange") {
                    QStringList xrange = strList.at(1).split(QRegExp("\\s+"), QString::SkipEmptyParts);
                }
                else if(variable == "IV%yrange") {
                    QStringList yrange = strList.at(1).split(QRegExp("\\s+"), QString::SkipEmptyParts);
                }
                else if(variable == "IV%smoothfactor") {
                    QStringList smoothfactor = strList.at(1).split(QRegExp("\\s+"), QString::SkipEmptyParts);
                }
                else if(variable == "IV%runOnCluster"){
                    QString yesValue = QString("'Y'");
                    if( value.compare( yesValue ) ){
                        runOnCluster = false;
                    } else {
                        runOnCluster = true;
                    }
                }
            }
        }

    }

    return success;
}

bool Optimisation::load(QString aerOptInputFilePath) {
    bool success = true;
    success &= readAerOptSettings(aerOptInputFilePath);
    success &= readFitness();
    success &= readProfilePointsFromSimulationDir();
    success &= readInitialBoundaryPoints();
    success &= readLogFromFile();

    return success;
}

QString Optimisation::logCacheFileName() {
    QString fileName = outputDataDirectory();
    fileName += "/output.log";
    fileName = QDir::toNativeSeparators(fileName);
    return fileName;
}

void Optimisation::addToOutputLog(const QString line) {

    mOutputLog += line;

    QString fileName = logCacheFileName();

    if( !runOnCluster ){
        // If running of the cluster, the outputfile is downloaded and
        // already contains all the output
        // write output to output file
        QFile logFile(fileName);
        if(logFile.open(QIODevice::WriteOnly | QIODevice::Append)) {
            QTextStream outputStream(&logFile);
            outputStream << line;
            logFile.close();
        } else {
            qWarning() << "Write output to log file " << fileName << " failed.";
        }
    }
}

bool Optimisation::readLogFromFile() {
    QString fileName = logCacheFileName();
    QFile logFile(fileName);
    QString str;

    // read output from output file
    if(logFile.open(QIODevice::ReadOnly)) {
        QTextStream inputStream(&logFile);
        mOutputLog = inputStream.readAll();
        logFile.close();
    } else {
        qWarning() << "Read output from log file " << fileName << " failed.";
        return false;
    }

    return true;
}

bool Optimisation::readInitialBoundaryPoints() {

    QString filePath = simulationDirectoryPath() + "/Boundary_Points.txt";
    bool success = true;
    std::ifstream infile(filePath.toStdString(), std::ifstream::in);
    success &= infile.is_open();

    QString qline, value;
    QStringList strList;
    std::string variable;

    if (success)
    {
        std::string line("");

        QStringList x;
        QStringList y;
        QStringList isCP;
        QStringList sf;
        QStringList xMin;
        QStringList xMax;
        QStringList yMin;
        QStringList yMax;

        while (std::getline(infile, line))
        {
            qline = QString::fromStdString(line);
            strList = qline.split(QRegExp("\\s*=\\s*"));

            if(strList.length()==2) {
                variable = strList.at(0).toStdString();
                value = strList.at(1);

                if(variable == "IV%x") {
                    x = strList.at(1).split("\t", QString::SkipEmptyParts);

                } else if(variable == "IV%y") {
                    y = strList.at(1).split("\t", QString::SkipEmptyParts);

                } else if(variable == "IV%isCP") {
                    isCP = strList.at(1).split("\t", QString::SkipEmptyParts);

                } else if(variable == "IV%smoothingFactor") {
                    sf = strList.at(1).split("\t", QString::SkipEmptyParts);

                } else if(variable == "IV%xMin") {
                    xMin = strList.at(1).split("\t", QString::SkipEmptyParts);

                } else if(variable == "IV%xMax") {
                    xMax = strList.at(1).split("\t", QString::SkipEmptyParts);

                } else if(variable == "IV%yMin") {
                    yMin = strList.at(1).split("\t", QString::SkipEmptyParts);

                } else if(variable == "IV%yMax") {
                    yMax = strList.at(1).split("\t", QString::SkipEmptyParts);

                }
            }
        }

        std::vector<BoundaryPoint*> boundaryPs;
        std::vector<BoundaryPoint*> controlPs;

        for (int i = 0; i< x.size(); i++) {
            QString xstr = x.at(i);
            QString ystr = y.at(i);
            QString isCPstr = isCP.at(i);
            QString sfstr = sf.at(i);
            QString xMinstr = xMin.at(i);
            QString xMaxstr = xMax.at(i);
            QString yMinstr = yMin.at(i);
            QString yMaxstr = yMax.at(i);
            double xVal = xstr.toDouble();
            double yVal = ystr.toDouble();
            bool isCPVal = false;
                if (isCPstr.toInt() == 1) isCPVal = true;
            double sfVal = sfstr.toDouble();
            double xMinVal = xMinstr.toDouble();
            double xMaxVal = xMaxstr.toDouble();
            double yMinVal = yMinstr.toDouble();
            double yMaxVal = yMaxstr.toDouble();

            BoundaryPoint* bp = new BoundaryPoint(xVal,yVal);
            bp->setControlPoint(isCPVal);
            bp->setSmoothing(sfVal);

            QRectF boundingBox = *new QRectF(xMinVal, yMinVal, xMaxVal - xMinVal, yMaxVal - yMinVal);

            bp->setControlPointRect(boundingBox);
            boundaryPs.push_back(bp);
            if (isCPVal){
                controlPs.push_back(bp);
            }

        }

        this->setBoundaryPoints(boundaryPs);
        this->setControlPoints(controlPs);

        // Add Bounding Boxes to Control Points
        QString inputFilePath = simulationDirectoryPath() + "/AerOpt_InputParameters.txt";
        std::ifstream bbInfile(inputFilePath.toStdString(), std::ifstream::in);
        success &= bbInfile.is_open();

    }
}
