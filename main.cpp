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
#include <QSettings>
#include <QSplashScreen>
#include <QMessageBox>

#include "DebugOutput.h"
#include "MainWindow.h"
#include "OptimisationModel.h"
#include "windowsizing.h"

void firstTimeSetup(QString AerOptWorkDir) {
    QSettings settings;
    settings.clear();

    if(AerOptWorkDir.indexOf(' ') >= 0) {
        QMessageBox::critical(nullptr, "Error", "Invalid Path: AerOpt executable path must not contains spaces.", QMessageBox::Ok);
        exit(-1);
    }

    AerOptWorkDir = QDir::toNativeSeparators(AerOptWorkDir);
    settings.setValue("AerOpt/workingDirectory", AerOptWorkDir);
    QDir().mkdir(AerOptWorkDir);

    settings.setValue("AerOpt/startLoadPath", AerOptWorkDir);

    // setup executables
    QString exeDir = QDir::toNativeSeparators(AerOptWorkDir + "Executables/");
    QDir().mkdir(exeDir);

    QString ext = ".exe";
    QString targetext = ".exe";
#ifdef Q_OS_UNIX
        ext = ".unix";
        targetext = "";
#endif
#ifdef Q_OS_MACOS
        ext = ".mac";
        targetext = "";
#endif

    bool copySuccess;

    QString mesherExe = QDir::toNativeSeparators(exeDir + "MeshGenerator" + targetext);
    copySuccess = QFile::copy(":/executables/MeshGenerator" + ext, mesherExe);
    if(copySuccess) {
        QFile(mesherExe).setPermissions(QFileDevice::ExeUser);
        settings.setValue("mesher/exe", mesherExe);
    } else {
        qCritical() << "Mesher file copy failed";
    }

    QString aeroptExe = QDir::toNativeSeparators(exeDir + "AerOpt" + targetext);
    copySuccess = QFile::copy(":/executables/AerOpt" + ext, aeroptExe);
    if(copySuccess) {
        QFile(aeroptExe).setPermissions(QFileDevice::ExeUser);
        settings.setValue("AerOpt/executable", aeroptExe);
    } else {
        qCritical() << "AerOpt executable file copy failed";
    }

    QString exePath = QDir::toNativeSeparators(exeDir + "Delaunay_2D" + targetext);
    QFile::copy(":/executables/Delaunay_2D" + ext, exePath);
    QFile(exePath).setPermissions(QFileDevice::ExeUser);

    exePath = QDir::toNativeSeparators(exeDir + "PrePro_2D" + targetext);
    QFile::copy(":/executables/PrePro_2D" + ext, exePath);
    QFile(exePath).setPermissions(QFileDevice::ExeUser);

    exePath = QDir::toNativeSeparators(exeDir + "Solver_2D" + targetext);
    QFile::copy(":/executables/Solver_2D" + ext, exePath);
    QFile(exePath).setPermissions(QFileDevice::ExeUser);

    qDebug() << "Working Directory:" << AerOptWorkDir;
    qDebug() << "Mesher Exe Path" << mesherExe;
    qDebug() << "AeroOpt Exe Path" << aeroptExe;

    bool c = true;
    QFileInfo checkMesher(mesherExe);
    QFileInfo checkAerOpt(aeroptExe);
    c &= checkMesher.exists();
    c &= checkMesher.isFile();
    c &= checkAerOpt.exists();
    c &= checkAerOpt.isFile();

    if (!c)
    {
        qCritical() << "Application executables not available!";
    }

    // setup default profiles
    QString profileDir = QDir::toNativeSeparators(AerOptWorkDir + "profiles/");
    QDir().mkdir(profileDir);
    settings.setValue("AerOpt/profilesDefaultPath", profileDir);

    settings.beginWriteArray("profiles");
    settings.setArrayIndex(0);
    QString profileFile = QDir::toNativeSeparators(profileDir + "NACA0024.prf");
    copySuccess = QFile::copy(":/profiles/NACA0024.prf", profileFile);
    if(copySuccess)
        settings.setValue("filepath", profileFile);

    settings.setArrayIndex(1);
    profileFile = QDir::toNativeSeparators(profileDir + "NACA21120.prf");
    copySuccess = QFile::copy(":/profiles/NACA21120.prf", profileFile);
    if(copySuccess)
        settings.setValue("filepath", profileFile);

    settings.endArray();

    // Input folder is read by AeroOpt executable
    QString inFolder = AerOptWorkDir + "Input_Data/";
    inFolder = QDir::toNativeSeparators(inFolder);
    settings.setValue("AerOpt/inFolder", inFolder);

    QString AerOptInFileName = "AerOpt_InputParameters.txt";
    QString AerOptInFile = inFolder + AerOptInFileName;
    AerOptInFile = QDir::toNativeSeparators(AerOptInFile);
    settings.setValue("AerOpt/inputFile", AerOptInFile);
    settings.setValue("AerOpt/inputFileName", AerOptInFileName);

    QString AerOptNodeFile = inFolder + "Control_Nodes.txt";
    AerOptNodeFile = QDir::toNativeSeparators(AerOptNodeFile);
    settings.setValue("AerOpt/nodeFile", AerOptNodeFile);

    QString AerOptBoundaryFile = inFolder + "Boundary_Points.txt";
    AerOptBoundaryFile = QDir::toNativeSeparators(AerOptBoundaryFile);
    settings.setValue("AerOpt/boundaryFile", AerOptBoundaryFile);

    // mesh is created in a scratch directory
    QString scratchPath = QDir::toNativeSeparators(AerOptWorkDir + "Scratch/");
    settings.setValue("mesher/scratchDir", scratchPath);

    QString initialMeshFile;
    initialMeshFile = QDir(scratchPath + "Mesh.dat").absolutePath();
    initialMeshFile = QDir::toNativeSeparators(initialMeshFile);
    settings.setValue("mesher/initMeshFile", initialMeshFile);

    settings.setValue("Cluster/Username", "");
    settings.setValue("Cluster/Account", "scw1022");
    settings.setValue("Cluster/Address", "sunbird.swansea.ac.uk");
}

void checkSettings()
{
    //Setup up default settings
    QCoreApplication::setOrganizationName("Swansea University (Engineering)");
    QCoreApplication::setOrganizationDomain("engineering.swansea.ac.uk");
    QCoreApplication::setApplicationName("AerOpt");

    // get working directory
    QString appPath = QDir::fromNativeSeparators(QCoreApplication::applicationFilePath());
    QFileInfo appPathInfo(appPath);

    QString AerOptWorkDir = appPathInfo.dir().absolutePath();

    // macOS: remove part of path inside application bundle
    AerOptWorkDir = AerOptWorkDir.remove(QRegExp("/AerOptGui.app/Contents.*"));

    AerOptWorkDir = QDir::fromNativeSeparators(AerOptWorkDir + "/AerOptFiles/");

    // if working directory exists, then no need for setup
    if(!QDir(AerOptWorkDir).exists()) {
        firstTimeSetup(AerOptWorkDir);
    }
}

//Programme entry point.
int main(int argc, char *argv[])
{
    //Initialises the QT library application event loops.
    QApplication app(argc, argv);
    app.setWindowIcon(QIcon(":/images/AerOpt.png"));

    checkSettings();

    //Application main interaction classes.
    //DebugOutput::Instance();
    qInfo() << " *** Welcome to AerOpt ***";
    qInfo() << " ***     Have a nice day     ***";

    // Setup optimisation model and selection model
    OptimisationModel* optimisationModel = new OptimisationModel();

    //Main window setup and show.
    MainWindow w;
    centerAndResizeWindow(&w);
    w.setOptimisationModel(optimisationModel);
    w.show();
    w.setSplitterSizes();

    return app.exec();
}
