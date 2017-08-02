#ifndef MESH_H
#define MESH_H

#include <QWidget>
#include <QProcess>
#include <QCoreApplication>
#include <QDir>
#include "Enumerations.h"
#include "CustomTypes.h"
#include "BoundaryPoint.h"

class Mesh : public QObject
{
    Q_OBJECT
public:
    Mesh(ProfilePoints profilePoints);

    // Getters
    bool getBoundaryLayerFlag();
    int getNumBoundaryLayers();
    qreal getBoundaryLayerThickness();
    /**
     * @brief getMeshConnectivities
     * @return A reference to the list of mesh connectivities.
     * Gets a reference to the mesh connectivity indices.
     */
    const MeshConnectivities& getMeshConnectivities() const;
    /**
     * @brief getBConnects
     * @return A reference to the current boundary connectivities.
     * Returns a reference to a list of pairs of indices representing the
     * boundary connectivities.
     */
    const BConnectivities& getBConnects() const;
    /**
     * @brief getMeshBoundary
     * @return A reference to the current mesh boundary.
     * Gets a reference to the current bounday data.
     */
    Boundaries &getMeshBoundary();
    /**
     * @brief getControlPoints
     * @return A reference to the list of control point indices.
     * Gets a reference to the current list of selected control points.
     */
    const std::list<unsigned int>& getControlPoints() const;

    // Setters
    void setBoundaryLayerFlag(bool hasBoundaryLayer);
    void setNumBoundaryLayers(int num_layers);
    void setBoundaryLayerThickness(qreal layer_thickness);
    void setMeshDensity(const Enum::Mesh& meshDensity);
    void stopMesher();

    // Other
    /**
     * @brief runMesher
     * run the mesher
     */
    void runMesher(QDir workDir);
    /**
     * @brief loadResults
     * loads the results
     */
    bool loadResults(const std::string& filePath);
    /**
     * @brief loadMesh
     * loads the mesh
     */
    bool loadMesh(const QString &filePath);
    /**
     * @brief loadMeshProfile
     * Loads the mesh profile
     */
    bool loadMeshProfile(const QString &filePath);
    /**
     * @brief clearMeshConnectivities
     * Clears the list of mesh connectivites.
     */
    void clearMeshConnectivities();
    /**
     * @brief selectControlPoint
     * @param index Index of the control point to be selected.
     * Sets the control point of interest as selected.
     */
    void selectControlPoint(const unsigned int& index);
    /**
     * @brief checkControlPointIntegrity
     * @return True if the control point is OK.
     * Checks the integrity of the control points.
     */
    bool checkControlPointIntegrity();
    //mesh points
    /**
     * @brief addMeshPoint
     * @param pair The point location pair as x and y.
     * Adds a point to the list of mesh points
     * where the pair is x and y locations.
     */
    void addMeshPoint(const std::pair<float,float>& pair);
    /**
     * @brief clearMeshPoints
     * Clears the current mesh points.
     */
    void clearMeshPoints();
    /**
     * @brief getMeshPoints
     * @return A reference to the list of x y pairs.
     * Gets a reference to the current list of mesh points.
     */
    const MeshPoints& getMeshPoints() const;
    //mesh connectivity
    /**
     * @brief addMeshConnectivity
     * @param tuple A triple of indices representing a triangular cell.
     * Adds mesh cell connectivity indeces.
     */
    void addMeshConnectivity(const std::tuple<uint, uint, uint>& tuple);
    //mesh results
    /**
     * @brief addMeshData
     * @param tuple A quintuple of results data (rho, u, v, e, p).
     * Sets the results data for each node or point in the mesh.
     */
    void addMeshData(const std::tuple<float, float, float, float, float>& tuple);
    /**
     * @brief clearMeshData
     * Clears the mesh results data.
     */
    void clearMeshData();
    /**
     * @brief getMeshData
     * @return A reference to the results data.
     * Gets a reference to the list of results data.
     */
    const MeshResults& getMeshData() const;
    /**
     * @brief getControlPoint
     * @param index Index of the boundary point
     * of interest.
     * @return A point and a rectangle representing
     * the control point of interest.
     * This function returns a reference to the control
     * point of interest.
     */
    BoundaryPoint &getControlPoint(const uint& index);
    /**
     * @brief checkBoundaryIntegrity
     * @return True if the boundary integrity is OK.
     * Checks the integrity of the boundary.
     */
    bool checkBoundaryIntegrity();
    /**
     * @brief hasMeshBoundary
     * @return Check if mesh boundary has been set
     */
    bool hasMeshBoundary() const;
    /**
     * @brief resetBoundary
     * Resets the boundary back to the initial boundary.
     */
    void resetBoundary();
    /**
     * @brief clear
     * clears all mesh data
     */
    void clear();
    //mesh boundary points
    /**
     * @brief addBoundaryPoint
     * @param gen The current generation number.
     * @param x Location x of point.
     * @param y Location x of point.
     * Adds a mesh boundary point.
     */
    void addBoundaryPoint(const float& x, const float& y);
    /**
     * @brief addBConnectivity
     * @param gen The current generation number.
     * @param a Boundary point index.
     * @param b Boundary point index.
     * Adds boundary point connectivity from a to b.
     */
    void addBConnectivity(const uint& a, const uint& b);

    Enum::Mesh getMeshDensity() const;

signals:
    void meshUpdate();

private slots:
    void meshingFinished(int exitCode, QProcess::ExitStatus exitStatus);

private:
    bool createInputFile(const std::string& meshInFile,
                         const std::string& meshBacFile,
                         const std::string& meshGeoFile,
                         const std::string& mMeshDatFile);
    bool createBacFile(const std::string& meshBacFile);
    bool createGeoFile(const std::string& meshGeoFile);
    bool loadMeshProfileType1(const QString& filePath);
    bool loadMeshProfileType2(const QString& filePath);
    void writeStdOutToLog();
    void writeStdErrToLog();
    /**
     * @brief loadMeshType1
     * @param filePath
     * @param data
     * @return
     */
    bool loadMeshType1(const QString& filePath);
    /**
     * @brief loadMeshType2
     * @param filePath
     * @param data
     * @return
     */
    bool loadMeshType2(const QString& filePath);

    // Mesher attributes
    QProcess mMeshProcess;
    QDir mMeshPath;
    Enum::Mesh mMeshDensity;
    QString mMeshDatFile;
    ProfilePoints mProfilePoints;

    //Mesh Attributes
    Boundaries mMeshProfile;
    BConnectivities mBoundConnects;
    std::list<uint> mControlPoints;
    MeshPoints mMeshPoints;
    MeshConnectivities mMeshConnectivities;
    MeshResults mMeshResults;

    // boundary layer parameters
    bool mHasBoundaryLayer;
    int mNumBoundaryLayers;
    qreal mBoundaryLayerThickness;

};

#endif // MESH_H





