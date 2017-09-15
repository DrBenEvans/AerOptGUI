#include "MeshDialogModel.h"
#include <QDebug>

MeshDialogModel::MeshDialogModel(QObject *parent) : QObject(parent),
  mMesh(new Mesh(this)),
  mBoundaryPointModel(new BoundaryPointModel(this))
{
    mMeshProcess.setParent(this);
}

void MeshDialogModel::runMesher() {

    QDir workDir("/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/Project");

    bool r = true;

    mMeshPath = workDir;
    mBoundaryPointModel->clearPoints();

    QString meshInFile;
    QString meshBacFile;
    QString meshGeoFile;

    QString mesherPath = "/Volumes/HardDrive/Users/mark/AerOpt/AerOpt/FLITE/Mesher/MeshGenerator";

    meshInFile = QDir(mMeshPath.absolutePath() + QDir::separator() + mMeshPath.dirName() + ".in").absolutePath();
    meshInFile = QDir::toNativeSeparators(meshInFile);

    meshBacFile = QDir(mMeshPath.absolutePath() + QDir::separator() + mMeshPath.dirName() + ".bac").absolutePath();
    meshBacFile = QDir::toNativeSeparators(meshBacFile);

    meshGeoFile = QDir(mMeshPath.absolutePath() + QDir::separator() + mMeshPath.dirName() + ".geo").absolutePath();
    meshGeoFile = QDir::toNativeSeparators(meshGeoFile);

    mMeshDatFile = QDir(mMeshPath.absolutePath() + QDir::separator() + mMeshPath.dirName() + ".dat").absolutePath();
    mMeshDatFile = QDir::toNativeSeparators(mMeshDatFile);

    r &= mMesh->createFiles(meshInFile, meshBacFile, meshGeoFile, mMeshDatFile);

    if (r)
    {
        mMeshProcess.setWorkingDirectory( mMeshPath.absolutePath() );
        mMeshProcess.setStandardInputFile( meshInFile );
        mMeshProcess.start( mesherPath );
        mMeshProcess.waitForFinished();
        meshingFinished(mMeshProcess.exitCode(),mMeshProcess.exitStatus());
    }
}

void MeshDialogModel::meshingFinished(int exitCode, QProcess::ExitStatus exitStatus)
{
    if (exitStatus == QProcess::NormalExit && exitCode == 0)
    {
        bool r = true;

        qInfo() << "Process finished normally";

        //Set and/or check existance of output data file
        r &= mMeshProcess.exitCode() == 0;
        r &= QFile::exists(mMeshDatFile);

        //Now load points into data object
        if (r)
        {
            r &= mMesh->loadMesh(mMeshDatFile);
        }

        //When done set menu to 'Yes'
        if (r)
        {
            qInfo() << "Mesh successfully created.";
        }
        else
        {
            qWarning() << "Mesh not created!";
        }

        mBoundaryPointModel->setPoints(mMesh->boundaryPoints());
        emit meshChanged();
    }
    else
    {
        qCritical() << "Process finished abnormally with exit code: " << exitCode;
        qCritical() << "Standard Ouput";
        writeStdOutToLog();
        qCritical() << "Standard Error";
        writeStdErrToLog();
    }
}

void MeshDialogModel::stopMesher()
{
    mMeshProcess.kill();
    qInfo() << "Any running mesh jobs have been stopped!";
}

void MeshDialogModel::writeStdOutToLog()
{
    QByteArray byteArray = mMeshProcess.readAllStandardOutput();
    QStringList strLines = QString(byteArray).split("\n");

    foreach (QString line, strLines){
        qInfo() << line;
    }
}

void MeshDialogModel::writeStdErrToLog()
{
    QByteArray byteArray = mMeshProcess.readAllStandardError();
    QStringList strLines = QString(byteArray).split("\n");

    foreach (QString line, strLines){
        qCritical() << line;
    }
}

Mesh* MeshDialogModel::currentMesh() {
    return mMesh;
}

BoundaryPointModel* MeshDialogModel::boundaryPointModel() {
    return mBoundaryPointModel;
}
