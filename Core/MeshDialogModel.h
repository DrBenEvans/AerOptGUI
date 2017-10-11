#ifndef MESHDIALOGMODEL_H
#define MESHDIALOGMODEL_H

#include <QObject>
#include "Mesh.h"
#include "MeshModel.h"
#include "ProfileModel.h"

class MeshDialogModel : public MeshModel
{
    Q_OBJECT
public:
    explicit MeshDialogModel(QObject *parent = nullptr);
    ~MeshDialogModel();

    bool runMesher();
    void stopMesher();

private slots:
    bool meshingFinished(int exitCode, QProcess::ExitStatus exitStatus, QString meshDatFile);

private:
    void writeStdOutToLog();
    void writeStdErrToLog();

    QProcess mMeshProcess;
};

#endif // MESHDIALOGMODEL_H
