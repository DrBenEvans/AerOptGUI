#ifndef MESHDIALOGMODEL_H
#define MESHDIALOGMODEL_H

#include <QObject>
#include "Mesh.h"

class MeshDialogModel : public QObject
{
    Q_OBJECT
public:
    explicit MeshDialogModel(QObject *parent = nullptr);

    void runMesher();
    void stopMesher();
    Mesh* currentMesh();
    BoundaryPointModel* boundaryPointModel();

signals:
    void meshChanged();

private slots:
    void meshingFinished(int exitCode, QProcess::ExitStatus exitStatus);

private:
    void writeStdOutToLog();
    void writeStdErrToLog();

    Mesh* mMesh;
    BoundaryPointModel* mBoundaryPointModel;
    QDir mMeshPath;
    QString mMeshDatFile;
    QProcess mMeshProcess;
};

#endif // MESHDIALOGMODEL_H
