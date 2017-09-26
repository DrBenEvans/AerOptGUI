#ifndef MESHMODEL_H
#define MESHMODEL_H

#include <QObject>
#include "Mesh.h"
#include "BoundaryPointModel.h"

class MeshModel : public QObject
{
    Q_OBJECT
public:
    MeshModel(QObject* parent);
    Mesh* currentMesh();
    BoundaryPointModel* boundaryPointModel();

signals:
    void meshChanged();

protected:
    Mesh* mMesh;
    BoundaryPointModel* mBoundaryPointModel;
};

#endif // MESHMODEL_H
