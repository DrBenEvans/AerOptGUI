#include "MeshDialogModel.h"
#include <QDebug>
#include <QSettings>
#include <QString>
#include "FileManipulation.h"

MeshDialogModel::MeshDialogModel(QObject *parent) : MeshModel(parent)
{
    mMeshProcess.setParent(this);
}

MeshDialogModel::~MeshDialogModel() {
    mMeshProcess.kill();
}

/**
 * @brief MeshDialogModel::runMesher Sets up the necessary files to beginning generating a mesh, then makes a method call to generate the mesh.
 * @return true iff mesh generated successfully, false if data files not created, or meshing process did not complete
 */
bool MeshDialogModel::runMesher() {
    QSettings settings;
    QString scratchDir = settings.value("mesher/scratchDir").toString();
    QString meshDatFile = settings.value("mesher/initMeshFile").toString();
    QString mesherPath = settings.value("mesher/exe").toString();

    FileManipulation::emptyFolder(scratchDir);

    mBoundaryPointModel->clearPoints();
    mMesh->clear();

    QString meshInFile;
    QString meshBacFile;
    QString meshGeoFile;

    meshInFile = QDir(scratchDir + "/scratch.in").absolutePath();
    meshInFile = QDir::toNativeSeparators(meshInFile);

    meshBacFile = QDir(scratchDir + "/scratch.bac").absolutePath();
    meshBacFile = QDir::toNativeSeparators(meshBacFile);

    meshGeoFile = QDir(scratchDir + "/scratch.geo").absolutePath();
    meshGeoFile = QDir::toNativeSeparators(meshGeoFile);

    //Try to create data files
    bool r = mMesh->createFiles(meshInFile, meshBacFile, meshGeoFile, meshDatFile);

    //If files were created successfully, then generate the mesh / start the mesher
    if (r)
    {
        mMeshProcess.setWorkingDirectory(scratchDir);
        mMeshProcess.setStandardInputFile(meshInFile);

        //Generate Mesh
        mMeshProcess.start(mesherPath);
        mMeshProcess.waitForFinished();
        r &= meshingFinished(mMeshProcess.exitCode(), mMeshProcess.exitStatus(), meshDatFile);
    }

    return r;
}

bool MeshDialogModel::meshingFinished(int exitCode, QProcess::ExitStatus exitStatus, QString meshDatFile)
{
    bool r = true;

    if (exitStatus == QProcess::NormalExit && exitCode == 0)
    {
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
            qWarning() << "Mesh not created! Check mesh file " << meshDatFile;
        }

        mBoundaryPointModel->setPoints(mMesh->boundaryPoints());
        emit meshChanged();
    }

    // If process did not complete normally then output error log
    else
    {
        qCritical() << "Process finished abnormally with exit code: " << exitCode;
        qCritical() << "Standard Ouput";
        writeStdOutToLog();
        qCritical() << "Standard Error";
        writeStdErrToLog();
        r = false;
    }

    return r;
}

/**
 * @brief MeshDialogModel::stopMesher
 * Stop meshing process immediately
 */
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
