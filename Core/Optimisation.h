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

#include <clusterManager.h>


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
     * @brief setBoundary Set list of mesh boundary points
	 * @param boundary Sets the boundary.
	 */
    void setBoundaryPoints(std::vector<BoundaryPoint*> boundaryPoints);

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


    bool readFitness();

    /**
     * @brief getProfile Returns the initial optimisation Profile
     * @return
     */
    ProfilePoints getProfile();
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

    int getNoTop() const;
    void setNoTop(int noTop);

    /**
     * @brief label Returns the identification label / optimisation name.
     * @return Optimisation identification label.
     */
    QString label() const;
    QString simulationDirectoryName();
    QString simulationDirectoryPath();

    /**
     * @brief setLabel Set the identification label
     * @param label Optimisation label.
     */
    void setLabel(QString label);
    Mesh *initMesh();

    // control nodes
    /**
     * @brief controlPoints Returns a vector of all Boundary Points that are control nodes.
     * @return Vector of Boundary Points that are control nodes
     */
    std::vector<BoundaryPoint*> controlPoints();

    /**
     * @brief boundaryPoints Returns a vector of all boundary points of the initial mesh.
     * @return
     */
    std::vector<BoundaryPoint*> initialBoundaryPoints();

    /**
     * @brief setControlPoints Sets the vector of control nodes.
     * @param controlPoints A vector of Boundary Points for which isControlPoint == true
     */
    void setControlPoints(std::vector<BoundaryPoint*> controlPoints);

    /**
     * @brief controlPointCount Returns the number of control points to optimise.
     * @return The number of control points
     */
    int controlPointCount();

    bool run();

    /**
     * @brief setModel Set the model that this optimisation belongs to.
     * @param model
     */
    void setModel(OptimisationModel* model);

    /**
     * @brief allfitness Returns a 2D vector of the fitnesses for all agents in each generation.
     * @return 2D Vector of fitness values
     */
    std::vector<std::vector<double>> allfitness();

    /**
     * @brief fitness Returns the fitness of a given agent within a specified generation.
     * @param generationIndex Generation number
     * @param agentIndex Agent index
     * @return The fitness of the given agent in the specified generation
     */
    double fitness(int generationIndex, int agentIndex);

    /**
     * @brief outputText Returns the output log.
     * @return The output log
     */
    QString outputText();

    /**
     * @brief fitnessRange Returns the minimum and maximum fitness values for this optimisation, across all generations.
     * @return A pair of double values of the format <minimum, maximum>.
     */
    std::pair<double,double> fitnessRange();

    /**
     * @brief fitnessRange Returns the minimum and maximum fitness values for the
     * best agents of each generation in this optimisation.
     * @return A pair of double values of the format <minimum, maximum>.
     */
    std::pair<double,double> fitnessRangeBestAgents();

    /**
     * @brief initProfilePoints Returns the initial profile.
     * @return mProfilePoints
     */
    ProfilePoints initProfilePoints();

    bool runOnCluster = false;
    QString mClusterPassword = "";

private:
    void optimiserFinished(int exitCode, QProcess::ExitStatus exitStatus);

    /**
     * @brief createAerOptInFile Writes optimisation data to an external file.
     * @param filePath File location
     * @return true iff file written successfully
     */
    bool createAerOptInFile(const QString &filePath);

    /**
     * @brief createAerOptNodeFile Writes control node coordinates to an external file.
     * @param filePath File location
     * @return true iff file written successfully
     */
    bool createAerOptNodeFile(const QString &filePath);

    /**
     * @brief createAerOptBoundaryPointFile Writes boundary point information for initial mesh.
     * @param filePath File location
     * @return true iff file written successfully
     */
    bool createAerOptBoundaryPointFile(const QString &filePath);
    bool saveCurrentProfile(const QString& path);
    QString outputDataDirectory();
    bool readAerOptSettings(QString filePath);

    /**
     * @brief addToOutputLog Adds a line to mOutputLog and the output log file.
     * @param line Text to add to output log
     */
    void addToOutputLog(const QString line);

    void writeProfilePointsToSimulationDir();
    bool readProfilePointsFromSimulationDir();

    /**
     * @brief setInitProfilePoints Set the initial profile.
     * @param profilePoints The initial profile.
     */
    void setInitProfilePoints(ProfilePoints profilePoints);

    /**
     * @brief readLogFromFile Load the output log into mOutputLog from the log file as defined in logCacheFileName().
     * @return true iff file loaded successfully
     */
    bool readLogFromFile();

    /**
     * @brief readInitialBoundaryPoints Load the initial boundary points.
     * @return true iff file loaded successfully.
     */
    bool readInitialBoundaryPoints();
    /**
     * @brief logCacheFileName Returns the file path for the output log file.
     * @return Output log file directory path.
     */
    QString logCacheFileName();

    // copy files
    void copyFileToSimulationDir(QString source);
    QString aerOptNodeFileCopyPath();
    QString aerOptInputFileCopyPath();

    /**
     * Identification label / optimisation name
    */
    QString mLabel = "";

	//Objective function attributes
    /**
     * @brief mObjFunc Fitness function that agents are being optimised for.
     */
    Enum::ObjFunc mObjFunc;

	//Boundary condition attributes
    /**
     * @brief mMachNo Mach number
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
     * @brief mOptimisationMethod Optimisation Algorithm
     */
    Enum::OptMethod mOptimisationMethod;

    /**
     * @brief mNoAgents Number of agents
     */
	int mNoAgents;

    /**
     * @brief mNoGens Number of generations
     */
	int mNoGens;

    /**
     * @brief mNoTop
     */
    int mNoTop;

    Mesh* mInitMesh;

    /**
     * @brief mControlPoints Vector of Boundary Points that are also optimisable control nodes.
     */
    std::vector<BoundaryPoint*> mControlPoints;

    /**
     * @brief mBoundaryPoints Vector of all Boundary Points.
     */
    std::vector<BoundaryPoint*> mBoundaryPoints;

    /**
     * @brief mFitness 2D vector of fitnesses for all agents across all generations.
     * Indexed by <generation number<agent index>>.
     */
    std::vector<std::vector<double>> mFitness;

    ProcessManager* mProcess = nullptr;

    /**
     * @brief mOptimisationModel Model that this optimisation belongs to.
     */
    OptimisationModel* mOptimisationModel;

    /**
     * @brief mOutputLog Process output log
     */
    QString mOutputLog = "";

    /**
     * @brief mProfilePoints The original profile before any optimisation.
     * E.g. NACA0024 or NACA21120
     */
    ProfilePoints mProfilePoints;

    clusterManager* clusterChecker = nullptr;
};



#endif // OptimisationRun_H
