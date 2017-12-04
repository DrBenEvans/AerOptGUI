#ifndef MESH_H
#define MESH_H

#include <QWidget>
#include <QProcess>
#include <QCoreApplication>
#include <QDir>
#include "Enumerations.h"
#include "CustomTypes.h"
#include "BoundaryPointModel.h"

class Mesh : public QObject
{
    Q_OBJECT
public:
    Mesh(QObject* parent = 0);

    // getters

    int numberBoundaryLayers() const;
    Enum::Mesh meshDensity() const;
    qreal boundaryLayerThickness() const;
    double growthFactor() const;
    ProfilePoints boundaryPoints();

    /**
     * @brief profilePoints
     * @return list of profile points from which the mesh has been or will be constructed
     */
    ProfilePoints profilePoints() const;

    /**
     * @brief getMeshPoints
     * @return A reference to the list of x y pairs.
     * Gets a reference to the current list of mesh points.
     */
    const MeshPoints& getMeshPoints() const;

    /**
     * @brief getBConnects
     * @return A reference to the current boundary connectivities.
     * Returns a reference to a list of pairs of indices representing the
     * boundary connectivities.
     */
    const BConnectivities& getBConnects() const;

    /**
     * @brief getMeshConnectivities
     * @return A reference to the list of mesh connectivities.
     * Gets a reference to the mesh connectivity indices.
     */
    const MeshConnectivities& getMeshConnectivities() const;

    /**
     * @brief getMeshData
     * @return A reference to the results data.
     * Gets a reference to the list of results data.
     */
    const MeshResults& getMeshData() const;

    // setters

    void setNumberBoundaryLayers(int num_layers);
    void setBoundaryLayerThickness(qreal layer_thickness);
    void setGrowthFactor(qreal factor);
    void setMeshDensity(const Enum::Mesh& meshDensity);
    void setProfilePoints(ProfilePoints profilePoints);

    // actions
    bool createFiles(QString meshInFile, QString meshBacFile, QString meshGeoFile, QString meshDatFile);

    /**
     * @brief loadMesh
     * loads the mesh
     */
    bool loadMesh(const QString &filePath);
    /**
     * @brief loadResults
     * loads the results
     */
    bool loadResults(const QString& filePath);

    /**
     * @brief clear
     * clears all mesh data
     */
    void clear();
private:

    // Other
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
     * @param index Index of the control node to be selected.
     * Sets the control node of interest as selected.
     */
    void selectControlPoint(const unsigned int& index);
    /**
     * @brief checkControlPointIntegrity
     * @return True if the control node is OK.
     * Checks the integrity of the control nodes.
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

private:
    bool createInputFile(const std::string& meshInFile,
                         const std::string& meshBacFile,
                         const std::string& meshGeoFile,
                         const std::string& meshDatFile);
    bool createBacFile(const std::string& meshBacFile);
    bool createGeoFile(const std::string& meshGeoFile);
    bool loadMeshProfileType1(const QString& filePath);
    bool loadMeshProfileType2(const QString& filePath);
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
    //mesh connectivity
    /**
     * @brief addMeshConnectivity
     * @param tuple A triple of indices representing a triangular cell.
     * Adds mesh cell connectivity indeces.
     */
    void addMeshConnectivity(const std::tuple<uint, uint, uint>& tuple);

    // Mesher attributes
    Enum::Mesh mMeshDensity;
    ProfilePoints mProfilePoints;

    //Mesh Attributes
    BConnectivities mBoundConnects;
    void resetBConnectivity();
    std::list<uint> mControlPoints;
    MeshPoints mMeshPoints;
    MeshConnectivities mMeshConnectivities;
    MeshResults mMeshResults;

    ProfilePoints mBoundaryPoints;

    // boundary layer parameters
    int mNumBoundaryLayers;
    qreal mGrowthFactor = 2.7182818285;
    qreal mBoundaryLayerThickness;

};

#endif // MESH_H





