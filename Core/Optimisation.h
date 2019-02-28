#ifndef OptimisationRun_H
#define OptimisationRun_H

#include <vector>
#include <list>
#include <utility>
#include <QRectF>
#include <memory>

#include "Enumerations.h"
#include "Profile.h"
#include "Mesh.h"
#include "ProcessManager.h"
#include "BoundaryPoint.h"
#include <QString>

class OptimisationModel;

/**
 * @brief The OptimisationRun class
 * This class stores project data
 */
class Optimisation
{
public:
	/**
     * @brief OptimisationRun
     * Constructor and destructor for
     * the OptimisationRun class.
	 */
    Optimisation();
    ~Optimisation();

    //Getters and Setters for the class variables
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
    /**
     * @brief readFitness
     * @param read the fitness for this optimisation
     */
    bool readFitness();
    /**
     * @brief mesh
     * @param getter method for meshes
     */
    Mesh *mesh(int genIndex, int agentIndex);
    /**
     * @brief load
     * @param load the optimisation from disk (based on label)
     */
    bool load(QString aerOptInputFilePath);

    /**
     * @brief getNoTop Return the percentage of lower-fitness agents that should be discarded during an iteration of the optimisation algorithm
     * @return percentage of lower-fitness agents that should be discarded
     */
    int getNoTop() const;

    /**
     * @brief setNoTop Set the percentage of lower-fitness agents that should be discarded during an iteration of the optimisation algorithm
     * @param noTop The percentage of lower-fitness agents that should be discarded
     */
    void setNoTop(int noTop);

    /**
     * @brief label Returns the name / label for the optimisation
     * @return Optimisation name
     */
    QString label() const;

    /**
     * @brief setLabel Set the label / name for the optimisation
     * @param label The label
     */
    void setLabel(QString label);


    /**
     * @brief simulationDirectoryName Returns a modified version of the label, suitable for directory naming.
     * Any non-alphanumeric character in label is replaced with '_'
     * @return Modified label
     */
    QString simulationDirectoryName();

    /**
     * @brief simulationDirectoryPath Returns the OS-specific directory path for the optimisation
     * @return the directory path for the optimisation
     */
    QString simulationDirectoryPath();

    /**
     * @brief initMesh Return the initial mesh.
     * @return the initial mesh
     */
    Mesh *initMesh();


    // control nodes
    std::vector<BoundaryPoint*> controlPoints();
    void setControlPoints(std::vector<BoundaryPoint*> controlPoints);
    int controlPointCount();

    /**
     * @brief Optimisation::run Run AerOpt on these optimisation parameters
     * @return true iff optimisation completed successfully, false if there was a file-handling error.
     */
    bool run();


    /**
     * @brief setModel Set the optimisation Model Handler for interacting with the UI.
     * @param model
     */
    void setModel(OptimisationModel* model);

    /**
     * @brief allfitness
     * @return a vector containing the fitness of each agent
     */
    std::vector<std::vector<double>> allfitness();

    /**
     * @brief fitness Returns the fitness for a given generation and agent
     * @param generationIndex Generation Index
     * @param agentIndex Agent Index
     * @return If agent and generation indices are valid, return the specified fitness value, otherwise return quiet_NaN
     */
    double fitness(int generationIndex, int agentIndex);

    /**
     * @brief outputText Return the output log
     * @return output log
     */
    QString outputText();

    /**
     * @brief fitnessRange Return the minimum and maximum fitness values present in an optimisation (for all generations)
     * @return First element of pair is the minimum fitness, second is the maximum
     */
    std::pair<double,double> fitnessRange();

    /**
     * @brief initProfilePoints Returns the list of points for the initial profile
     * @return The list of points for the initial profile
     */
    ProfilePoints initProfilePoints();

private:

    void optimiserFinished(int exitCode, QProcess::ExitStatus exitStatus);
    bool createAerOptInFile(const QString &filePath);
    bool createAerOptNodeFile(const QString &filePath);
    bool saveCurrentProfile(const QString& path);
    QString outputDataDirectory();
    bool readAerOptSettings(QString filePath);
    void addToOutputLog(const QString line);

    void writeProfilePointsToSimulationDir();
    bool readProfilePointsFromSimulationDir();
    void setInitProfilePoints(ProfilePoints profilePoints);

    bool readLogFromFile();
    QString logCacheFileName();

    // copy files
    void copyFileToSimulationDir(QString source);
    QString aerOptNodeFileCopyPath();
    QString aerOptInputFileCopyPath();

    QString mLabel = "";

	//Objective function attributes
    /**
     * @brief mObjFunc The function to be optimised.
     */
    Enum::ObjFunc mObjFunc;


	//Boundary condition attributes
    /**
     * @brief mMachNo Mach Number
     */
	float mMachNo;
    /**
     * @brief mReNo Reynolds Number
     */
	float mReNo;
    /**
     * @brief mFreeAlpha Angle of Attack
     */
	float mFreeAlpha;
    /**
     * @brief mFreePress Pressure Absolute
     */
	float mFreePress;
    /**
     * @brief mFreeTemp Temperature Absolute
     */
	float mFreeTemp;


	//Optimiser parameters
    /**
     * @brief mOptimisationMethod The optimisation algorithm
     */
    Enum::OptMethod mOptimisationMethod;
    /**
     * @brief mNoAgents The number of agents (solutions) per generation
     */
	int mNoAgents;
    /**
     * @brief mNoGens The number of generations for which the optimisation algorithm should iterate
     */
	int mNoGens;

    /**
     * @brief mNoTop The percentage of low-fitness agents that should be discarded
     */
    int mNoTop;

    Mesh* mInitMesh;
    std::vector<BoundaryPoint*> mControlPoints;

    /**
     * @brief mFitness Vector containing the fitness of each agent
     */
    std::vector<std::vector<double>> mFitness;

    ProcessManager* mProcess = nullptr;

    /**
     * @brief mOptimisationModel Optimisation Model Handler for interacting with the UI.
     */
    OptimisationModel* mOptimisationModel;

    /**
     * @brief mOutputLog Complete output log for optimisation process.
     */
    QString mOutputLog = "";

    /**
     * @brief mProfilePoints A list of the points of the initial profile.
     */
    ProfilePoints mProfilePoints;
};

#endif // OptimisationRun_H
