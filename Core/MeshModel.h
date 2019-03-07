#ifndef MESHMODEL_H
#define MESHMODEL_H

#include <QObject>
#include "Mesh.h"
#include "BoundaryPointModel.h"

/**
 * @brief The MeshModel class allows you to interact with a mesh and a set of boundary points where the mesh meets the profile.
 */
class MeshModel : public QObject
{
    Q_OBJECT
public:
    /**
     * @brief MeshModel Constructor Method.
     * Initialises with an empty boundary point model and no mesh.
     * @param parent
     */
    MeshModel(QObject* parent);

    /**
     * @brief currentMesh Returns the current mesh.
     * @return the current mesh.
     */
    Mesh* currentMesh();

    /**
     * @brief setCurrentMesh Replaces the current mesh with a new one. The previous mesh is discarded.
     * Emits meshChanged().
     * @param mesh The new mesh.
     */
    void setCurrentMesh(Mesh* mesh);

    /**
     * @brief boundaryPointModel Returns the Boundary Points Model for intertacting
     * with the points where the mesh meets the profile.
     * @return Boundary Point Model
     */
    BoundaryPointModel* boundaryPointModel();

signals:
    /**
     * @brief meshChanged Emits if a new mesh is present.
     */
    void meshChanged();

protected:
    /**
     * @brief mMesh The mesh.
     */
    Mesh* mMesh;
    /**
     * @brief mBoundaryPointModel The model for interacting with Boundary Points between the mesh nad profile.
     */
    BoundaryPointModel* mBoundaryPointModel;
};

#endif // MESHMODEL_H
