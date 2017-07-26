/*********************************************
**
**	Created on: 	08/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		TreeView.h
**
**********************************************/

#ifndef TREEVIEW_H
#define TREEVIEW_H

#include "ui_TreeView.h"
#include <QProcess>
#include <QFileSystemWatcher>
#include <QDialog>

//Forward declarations
class Menu;
class OptimisationRun;
class QTreeWidgetItem;
class Canvas;
class PlotterDialog;

/**
 * @brief The AppController class
 * This class ised used for interecting with the
 * user via a tree type menu structure.
 */
class AppController : public QDialog
{
    Q_OBJECT
	
public:
	/**
	 * @brief TreeView
     * @param data A reference to the OptimisationRun class.
	 * @param canvas A reference to the Canvas class.
     * This class depends on the OptimisationRun and Canvas classes.
	 */
    AppController(OptimisationRun& data, Canvas& canvas, QWidget *parent = 0);
    ~AppController();

public slots:
	void clearProject();
	/**
	 * @brief loadProject
	 * Sets the project working directory from user input.
	 */
	void loadProject();

	//Sub menus
	/**
	 * @brief runMesher
	 * Runs the mesher tool to generate an initial mesh
	 * from the imported profile stored in the project
	 * working directory, and input from the mesh dialogue.
	 */
	void runMesher();//2
    /**
	 * @brief setBoundary
     * Sets the flow conditions from user input
	 * from boundary condition dialogue.
	 */
	void setBoundary();//4
	/**
	 * @brief setOptimiser
	 * Sets the optimiser parameters from the user
	 * input optimiser dialogue.
	 */
	void setOptimiser();//5
	/**
	 * @brief runAerOpt
	 * Runs the aeropt programme in a new thread.
	 */
	void runAerOpt();//6
	/**
	 * @brief stopMesher
	 * Stops any running mesher threads.
	 */
	void stopMesher();//2
	/**
	 * @brief stopAerOpt
	 * Stops any running aeropt threads.
	 */
	void stopAerOpt();//6
	/**
	 * @brief showGraph
	 * Displays the fitness graph.
	 */
	void showGraph();
    /**
     * @brief createNewSimulationObject
     * @return
     */
    void new_simulation();
    /**
     * @brief configureCurrentSimulationObject
     * @return
     */
    void configure_current_simulation();

private slots:
	/**
	 * @brief meshOutput
	 * Outputs mesher std out.
	 */
	void meshOutput();
	/**
	 * @brief meshError
	 * Outputs mesher std err.
	 */
	void meshError();
	/**
	 * @brief processOutput
	 * Outputs aeropt std out.
	 */
	void processOutput();
	/**
	 * @brief processError
	 * Outputs aeropt std err.
	 */
	void processError();
	/**
	 * @brief meshingStarted
	 * Called at mesher start.
	 */
	void meshingStarted();
	/**
	 * @brief meshingFinished
	 * @param exitCode mesher exit code.
	 * @param exitStatus mesher exit status.
	 * Called at mesher finished.
	 */
	void meshingFinished(int exitCode, QProcess::ExitStatus exitStatus);
	/**
	 * @brief optimiserStarted
	 * Caller at optimiser start.
	 */
	void optimiserStarted();
	/**
	 * @brief readDirectory
	 * @param path The directory to read.
	 * Reads the directory specified directory
	 * checking for aeropt output files.
	 */
	void readDirectory(const QString& path);
	/**
	 * @brief readFitness
	 * @param path The fitness file to read.
	 * Reads the contents of a fitness file.
	 */
	void readFitness(const QString& path);
	/**
	 * @brief optimiserFinished
	 * @param exitCode Optimiser exit code.
	 * @param exitStatus Optimiser exit status.
	 */
    void optimiserFinished(int exitCode, QProcess::ExitStatus exitStatus);

private:
	/**
	 * @brief loadProfile
	 * @param filePath
	 * @param data
	 * @return True if file loaded sucsessfully.
	 * Loads the profile data from file specified by path,
     * and stores the data in OptimisationRun.
	 */
    bool loadProfile(const std::string& filePath, OptimisationRun& data);
	/**
	 * @brief createInputFile
	 * @param meshInFile
	 * @param meshBacFile
	 * @param meshGeoFile
	 * @param meshDatFile
	 * @return
	 */
	bool createInputFile(const std::string& meshInFile,
						 const std::string& meshBacFile,
						 const std::string& meshGeoFile,
						 const std::string& meshDatFile);
	/**
	 * @brief createBacFile
	 * @param meshBacFile
	 * @return
	 */
	bool createBacFile(const std::string& meshBacFile);
	/**
	 * @brief createGeoFile
	 * @param meshGeoFile
	 * @param data
	 * @return
	 */
    bool createGeoFile(const std::string& meshGeoFile, OptimisationRun& data);
	/**
	 * @brief loadMeshProfile
	 * @param genNo
	 * @param filePath
	 * @param data
	 * @return
	 */
    bool loadMeshProfile(const uint genNo, const std::string& filePath, OptimisationRun& data);
	/**
	 * @brief loadMeshProfileType1
	 * @param genNo
	 * @param filePath
	 * @param data
	 * @return
	 */
    bool loadMeshProfileType1(const uint genNo, const std::string& filePath, OptimisationRun& data);
	/**
	 * @brief loadMeshProfileType2
	 * @param genNo
	 * @param filePath
	 * @param data
	 * @return
	 */
    bool loadMeshProfileType2(const uint genNo, const std::string& filePath, OptimisationRun& data);
	/**
	 * @brief loadMesh
	 * @param filePath
	 * @param data
	 * @return
	 */
    bool loadMesh(const std::string& filePath, OptimisationRun& data);
	/**
	 * @brief loadMeshType1
	 * @param filePath
	 * @param data
	 * @return
	 */
    bool loadMeshType1(const std::string& filePath, OptimisationRun& data);
	/**
	 * @brief loadMeshType2
	 * @param filePath
	 * @param data
	 * @return
	 */
    bool loadMeshType2(const std::string& filePath, OptimisationRun& data);
	/**
	 * @brief loadResults
	 * @param filePath
	 * @param data
	 * @return
	 */
    bool loadResults(const std::string& filePath, OptimisationRun& data);
	/**
	 * @brief createAerOptInFile
	 * @param filePath
	 * @param data
	 * @return
	 */
    bool createAerOptInFile(const std::string& filePath, OptimisationRun& data);
	/**
	 * @brief createAerOptNodeFile
	 * @param filePath
	 * @param data
	 * @return
	 */
    bool createAerOptNodeFile(const std::string& filePath, OptimisationRun& data);

	/**
	 * @brief emptyFolder
	 * @param path
	 * @return
	 */
	bool emptyFolder(const QString& path);
	/**
	 * @brief copyFolder
	 * @param source
	 * @param dest
	 * @return
	 */
	bool copyFolder(const QString& source, const QString& dest);
	/**
	 * @brief copyFile
	 * @param source
	 * @param dest
	 * @return
	 */
	bool copyFile(const QString& source, const QString& dest);
	/**
	 * @brief removeFolder
	 * @param path
	 * @return
	 */
	bool removeFolder(const QString& path);

	/**
	 * @brief saveCurrentProfile
	 * @param path
	 * @param data
	 * @return
	 */
    bool saveCurrentProfile(const QString& path, const OptimisationRun& data);


	QTreeWidgetItem* mRoot;

	uint sGenNo;

    QObject* mParent;

	QProcess myMeshProcess;
	QProcess myOptProcess;
	QFileSystemWatcher myDirWatcher;

	QString mAppPath;
	QString mMesherPath;
	QString mAerOptPath;
	QString mProjectDirectory;
	QTreeWidgetItem* mCurrentNode;

	PlotterDialog* mPlotter;

	Canvas& mCanvas;
    OptimisationRun& mData;
};

#endif // TREEVIEW_H
