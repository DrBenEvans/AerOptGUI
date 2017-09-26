#include "MeshDialogModel.h"
#include <QDebug>
#include <QSettings>

MeshDialogModel::MeshDialogModel(QObject *parent) : MeshModel(parent)
{
    mMeshProcess.setParent(this);
    mMesh = new Mesh(this);
}

MeshDialogModel::~MeshDialogModel() {
    mMeshProcess.kill();
}

void MeshDialogModel::runMesher() {
    QSettings settings;
    QDir scratchDir = QDir(settings.value("mesher/scratchDir").toString());
    QString meshDatFile = settings.value("mesher/datFile").toString();

    bool r = true;

    mBoundaryPointModel->clearPoints();
    mMesh->clear();

    QString meshInFile;
    QString meshBacFile;
    QString meshGeoFile;

    QString mesherPath = settings.value("mesher/exe").toString();

    meshInFile = QDir(scratchDir.absolutePath() + "/scratch.in").absolutePath();
    meshInFile = QDir::toNativeSeparators(meshInFile);

    meshBacFile = QDir(scratchDir.absolutePath() + "/scratch.bac").absolutePath();
    meshBacFile = QDir::toNativeSeparators(meshBacFile);

    meshGeoFile = QDir(scratchDir.absolutePath() + "/scratch.geo").absolutePath();
    meshGeoFile = QDir::toNativeSeparators(meshGeoFile);

    r &= mMesh->createFiles(meshInFile, meshBacFile, meshGeoFile, meshDatFile);

    if (r)
    {
        mMeshProcess.setWorkingDirectory( mMeshPath.absolutePath() );
        mMeshProcess.setStandardInputFile( meshInFile );
        mMeshProcess.start( mesherPath );
        mMeshProcess.waitForFinished();
        meshingFinished(mMeshProcess.exitCode(), mMeshProcess.exitStatus(), meshDatFile);
    }
}

void MeshDialogModel::meshingFinished(int exitCode, QProcess::ExitStatus exitStatus, QString meshDatFile)
{
    if (exitStatus == QProcess::NormalExit && exitCode == 0)
    {
        bool r = true;

        qInfo() << "Process finished normally";

        //Set and/or check existance of output data file
        r &= mMeshProcess.exitCode() == 0;
        r &= QFile::exists(meshDatFile);

        //Now load points into data object
        if (r)
        {
            r &= mMesh->loadMesh(meshDatFile);
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
