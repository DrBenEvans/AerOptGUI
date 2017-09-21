#ifndef MESHDIALOGMODEL_H
#define MESHDIALOGMODEL_H

#include <QObject>
#include "Mesh.h"
#include "ProfileModel.h"

class MeshDialogModel : public QObject
{
    Q_OBJECT
public:
    explicit MeshDialogModel(QObject *parent = nullptr);
    ~MeshDialogModel();

    void runMesher();
    void stopMesher();
    Mesh* currentMesh();
    BoundaryPointModel* boundaryPointModel();

signals:
    void meshChanged();

private slots:
    void meshingFinished(int exitCode, QProcess::ExitStatus exitStatus, QString meshDatFile);

private:
    void writeStdOutToLog();
    void writeStdErrToLog();

    Mesh* mMesh;
    BoundaryPointModel* mBoundaryPointModel;
    QDir mMeshPath;
    QProcess mMeshProcess;
};

#endif // MESHDIALOGMODEL_H
