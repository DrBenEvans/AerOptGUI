/*********************************************
**
**	Created on: 	09/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		ProjectData.h
**
**********************************************/

#ifndef PROJECTDATA_H
#define PROJECTDATA_H

#include <vector>
#include <list>
#include <utility>
#include <QRectF>

#include "Enumerations.h"
#include "BoundaryPoint.h"

//Abbreviate long type names
typedef std::vector<std::vector<BoundaryPoint*>> Boundaries;
typedef std::vector<std::vector<std::pair<uint,uint>>> BConnectivities;
typedef std::vector<std::pair<float,float>> MeshPoints;
typedef std::vector<std::tuple<uint,uint,uint>> MeshConnectivities;
typedef std::vector<std::tuple<float,float,float,float,float>> MeshResults;

/**
 * @brief The ProjectData class
 * This class stores all aspects of the
 * project data accessed from the treeview
 * and canvas classes.
 */
class ProjectData
{
public:
	/**
	 * @brief ProjectData
	 * Constructor and distructor for
	 * the ProjectData class.
	 */
	ProjectData();
	~ProjectData();

	/**
	 * @brief clearProject
	 * Clears all project data.
	 */
	void clearProject();

	//Profile Attributes
	/**
	 * @brief clearProfile
	 * Clears the current profile data.
	 */
	void clearProfile();
	/**
	 * @brief addPoint
	 * @param x Dimension x location of point.
	 * @param y Dimension y location of point.
	 * Adds points to the current profile.
	 */
	void addPoint(const float& x, const float& y);
	/**
	 * @brief checkProfileIntegrity
	 * @return True if integrity of profile is OK.
	 * Checks the integrity of the current profile.
	 */
	bool checkProfileIntegrity();
	/**
	 * @brief getProfile
	 * @return A reference to the profile as a series of x y pairs.
	 */
	const std::list<std::pair<float,float>>& getProfile() const;


	//Mesh Attributes
	/**
	 * @brief clearMesh
	 * Clears the mesh data.
	 */
	void clearMesh();
	//mesh boundary points
	/**
	 * @brief addBoundaryPoint
	 * @param gen The current generation number.
	 * @param x Location x of point.
	 * @param y Location x of point.
	 * Adds a mesh boundary point.
	 */
	void addBoundaryPoint(const uint& gen, const float& x, const float& y);
	/**
	 * @brief addBConnectivity
	 * @param gen The current generation number.
	 * @param a Boundary point index.
	 * @param b Boundary point index.
	 * Adds boundary point connectivity from a to b.
	 */
	void addBConnectivity(const uint& gen, const uint& a, const uint& b);
	/**
	 * @brief getControlPoint
	 * @param genNo The generation of interest.
	 * @param index Index of the boundary point
	 * of interest.
	 * @return A point and a rectangle representing
	 * the control point of interest.
	 * This function returns a reference to the control
	 * point of interest.
	 */
    BoundaryPoint *getControlPoint(const uint& genNo, const uint& index);
	/**
	 * @brief checkBoundaryIntegrity
	 * @return True if the boundary integrity is OK.
	 * Checks the integrity of the boundary.
	 */
	bool checkBoundaryIntegrity();
	/**
	 * @brief getBoundary
	 * @return A reference to the current boundary.
	 * Gets a reference to the current bounday data.
	 */
	const Boundaries& getBoundary() const;
	/**
	 * @brief getBConnects
	 * @return A reference to the current boundary connectivities.
	 * Returns a reference to a list of pairs of indices representing the
	 * boundary connectivities.
	 */
	const BConnectivities& getBConnects() const;
	/**
	 * @brief resetBoundary
	 * Resets the boundary back to the initial boundary.
	 */
	void resetBoundary();
	//control points
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
	/**
	 * @brief getControlPoints
	 * @return A reference to the list of control point indices.
	 * Gets a reference to the current list of selected control points.
	 */
	const std::list<unsigned int>& getControlPoints() const;
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
	/**
	 * @brief clearMeshConnectivities
	 * Clears the list of mesh connectivites.
	 */
	void clearMeshConnectivities();
	/**
	 * @brief getMeshConnectivities
	 * @return A reference to the list of mesh connectivities.
	 * Gets a reference to the mesh connectivity indices.
	 */
	const MeshConnectivities& getMeshConnectivities() const;
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

	//Getters and Setters for the class variables
	/**
	 * @brief profile
	 * @return True if profile is set.
	 */
	bool profile() const;
	/**
	 * @brief setProfile
	 * @param profile Sets profile.
	 */
	void setProfile(bool profile);
	/**
	 * @brief mesh
	 * @return True if profile is set.
	 */
	bool mesh() const;
	/**
	 * @brief setMesh
	 * @param mesh Sets mesh.
	 */
	void setMesh(bool mesh);
    /**
	 * @brief setFunction
	 * @param function Sets the function.
	 */
	void setFunction(bool function);
	/**
	 * @brief boundary
	 * @return True if boundary is set.
	 */
	bool boundary() const;
	/**
	 * @brief setBoundary
	 * @param boundary Sets the boundary.
	 */
	void setBoundary(bool boundary);
	/**
	 * @brief optimiser
	 * @return True if optimiser is set.
	 */
	bool optimiser() const;
	/**
	 * @brief setOptimiser
	 * @param optimiser Sets optimiser.
	 */
    void setOptimisationMethod(int method_index);
    /**
     * @brief setOptimisation method
     * @param method_index is the method index as defined by ordering in OptimiserDialog.ui.
     */
	void setOptimiser(bool optimiser);
	/**
	 * @brief runTime
	 * @return True if runtime is set.
	 */
	bool runTime() const;
	/**
	 * @brief setRunTime
	 * @param runTime
	 */
	void setRunTime(bool runTime);
	/**
	 * @brief renderProfile
	 * @return True if render profile is set.
	 */
	bool renderProfile() const;
	/**
	 * @brief setRenderProfile
	 * @param renderProfile Sets render profile.
	 */
	void setRenderProfile(bool renderProfile);
	/**
	 * @brief renderMesh
	 * @return True if render mesh is set.
	 */
	bool renderMesh() const;
	/**
	 * @brief setRenderMesh
	 * @param renderMesh Sets render mesh.
	 */
	void setRenderMesh(bool renderMesh);
	/**
	 * @brief objFunc
	 * @return The current objective function.
	 */
	Enum::ObjFunc objFunc() const;
	/**
	 * @brief setObjFunc
	 * @param objFunc sets the current objective function.
	 */
	void setObjFunc(const Enum::ObjFunc& objFunc);
	/**
	 * @brief machNo
	 * @return Gets the current mach number.
	 */
	float machNo() const;
	/**
	 * @brief setMachNo
	 * @param machNo Sets the current mach number.
	 */
	void setMachNo(float machNo);
	/**
	 * @brief reNo
	 * @return Gets the current Reynalds number.
	 */
	float reNo() const;
	/**
	 * @brief setReNo
	 * @param reNo Sets the current reynalds number.
	 */
	void setReNo(float reNo);
	/**
	 * @brief freeAlpha
	 * @return Gets the current angle of attack.
	 */
	float freeAlpha() const;
	/**
	 * @brief setFreeAlpha
	 * @param freeAlpha Sets the current angle of attack.
	 */
	void setFreeAlpha(float freeAlpha);
	/**
	 * @brief freePress
	 * @return Gets the current free pressure.
	 */
	float freePress() const;
	/**
	 * @brief setFreePress
	 * @param freePress Sets the current pree pressure.
	 */
	void setFreePress(float freePress);
	/**
	 * @brief freeTemp
	 * @return Gets the current free temperature.
	 */
	float freeTemp() const;
	/**
	 * @brief setFreeTemp
	 * @param freeTemp Sets the current free temperature.
	 */
	void setFreeTemp(float freeTemp);
	/**
	 * @brief noAgents
	 * @return Gets the number of agents.
	 */
	int noAgents() const;
	/**
	 * @brief setNoAgents
	 * @param noAgents Sets the number of agents.
	 */
	void setNoAgents(int noAgents);
	/**
	 * @brief noGens
	 * @return Gets the number of generations.
	 */
	int noGens() const;
	/**
	 * @brief setNoGens
	 * @param noGens Sets the number of generations.
	 */
	void setNoGens(int noGens);
	/**
	 * @brief meshDensity
	 * @return Gets the current mesh density requirement.
	 */
	Enum::Mesh meshDensity() const;
	/**
	 * @brief setMeshDensity
	 * @param meshDensity Sets the current mesh density requirement.
	 */
	void setMeshDensity(const Enum::Mesh& meshDensity);
	/**
	 * @brief projectPathSet
	 * @return True if the project path is set.
	 */
	bool projectPathSet() const;
	/**
	 * @brief setProjectPathSet
	 * @param projectPathSet Sets wether the project path is set.
	 */
	void setProjectPathSet(bool projectPathSet);

    int getOptimisationMethod() const;
    /**
     * @brief getOptimisation method
     * @param returns the method index as defined by ordering in OptimiserDialog.ui.
     */
    int getNoTop() const;
    void setNoTop(int noTop);

    bool getBoundaryLayerFlag();
    int getNumBoundaryLayers();
    qreal getBoundaryLayerThickness();
    void setBoundaryLayerFlag(bool hasBoundaryLayer);
    void setNumBoundaryLayers(int num_layers);
    void setBoundaryLayerThickness(qreal layer_thickness);

private:
    bool mProfileSet;
    bool mProjectPathSet;
	bool mMeshSet;
	bool mBoundarySet;
	bool mOptimiserSet;
	bool mRunTimeSet;
	bool mRenderProfile;
	bool mRenderMesh;

	//Profile Attributes
	bool checkDuplicates();
	bool checkClockwise();
	bool checkNormalised();
	std::list<std::pair<float,float>> mProfile;

	//Mesh Attributes
	Boundaries mMeshProfile;
	BConnectivities mBoundConnects;
	std::list<uint> mControlPoints;
	MeshPoints mMeshPoints;
	MeshConnectivities mMeshConnectivities;
	MeshResults mMeshResults;
	Enum::Mesh mMeshDensity;

	//Objective function attributes
	Enum::ObjFunc mObjFunc;

	//Boundary condition attributes
	float mMachNo;
	float mReNo;
	float mFreeAlpha;
	float mFreePress;
	float mFreeTemp;

	//Optimiser parameters
    int mOptimisationMethod;
	int mNoAgents;
	int mNoGens;
    int mNoTop;

    // boundary layer parameters
    bool mHasBoundaryLayer;
    int mNumBoundaryLayers;
    qreal mBoundaryLayerThickness;

};

#endif // PROJECTDATA_H
