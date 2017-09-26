#include "MeshModel.h"

MeshModel::MeshModel(QObject *parent) :
    QObject(parent),
    mMesh(nullptr),
    mBoundaryPointModel(new BoundaryPointModel)
{

}

Mesh* MeshModel::currentMesh() {
    return mMesh;
}

BoundaryPointModel* MeshModel::boundaryPointModel() {
    return mBoundaryPointModel;
}
