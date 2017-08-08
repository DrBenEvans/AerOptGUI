/*********************************************
**
**	Created on: 	09/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		OptimisationRun.h
**
**********************************************/

#ifndef OptimisationRun_H
#define OptimisationRun_H

#include <vector>
#include <list>
#include <utility>
#include <QRectF>

#include "Enumerations.h"
#include "Profile.h"
#include "Mesh.h"
#include "BoundaryPoint.h"
#include <QString>

/**
 * @brief The OptimisationRun class
 * This class stores project data
 */
class OptimisationRun
{
public:
	/**
     * @brief OptimisationRun
     * Constructor and destructor for
     * the OptimisationRun class.
	 */
    OptimisationRun();
    ~OptimisationRun();

	/**
	 * @brief clearProject
	 * Clears all project data.
	 */
	void clearProject();

	//Mesh Attributes
	/**
	 * @brief clearMesh
	 * Clears the mesh data.
	 */
	void clearMesh();

    //Getters and Setters for the class variables
    /**
     * @brief setProfile
     * @param Profile object to use
     */
    void setProfile(QSharedPointer<Profile> profile);
    /**
	 * @brief mesh
	 * @return True if profile is set.
	 */
	bool mesh() const;
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
    void setOptimisationMethod(Enum::OptMethod method);
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
     * @brief getOptimisationMethod
     * @param returns the method index as defined by ordering in OptimiserDialog.ui.
     */
    Enum::OptMethod getOptimisationMethod() const;
    QSharedPointer<Mesh> newMesh();
    int getNoTop() const;
    void setNoTop(int noTop);

    QString getLabel();
    void setLabel(QString label);

    const std::list<std::pair<float,float>> getProfile() const;
    QSharedPointer<Profile> getProfileObj() const;

    QSharedPointer<Mesh> getMesh();
    void finishConfigure();

private:
    QSharedPointer<Profile> mProfile;

    QString mLabel;

	//Objective function attributes
	Enum::ObjFunc mObjFunc;

	//Boundary condition attributes
	float mMachNo;
	float mReNo;
	float mFreeAlpha;
	float mFreePress;
	float mFreeTemp;

	//Optimiser parameters
    Enum::OptMethod mOptimisationMethod;
	int mNoAgents;
	int mNoGens;
    int mNoTop;

    QSharedPointer<Mesh> mMesh;

};

#endif // OptimisationRun_H
