/*********************************************
**
**	Created on: 	08/04/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:		TreeView.h
**
**********************************************/

#ifndef TREEVIEW_H
#define TREEVIEW_H

#include <QFileSystemWatcher>
#include <QDialog>
#include "Mesh.h"

//Forward declarations
class Menu;
class OptimisationRun;
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
	 * @brief runAerOpt
	 * Runs the aeropt programme in a new thread.
	 */
	void runAerOpt();//6
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

private slots:
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
	 * @brief loadMesh
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
    bool saveCurrentProfile(const QString& path, OptimisationRun &data);


	uint sGenNo;

    QObject* mParent;

	QProcess myOptProcess;
	QFileSystemWatcher myDirWatcher;

	QString mAppPath;
	QString mMesherPath;
	QString mAerOptPath;
	QString mProjectDirectory;

	PlotterDialog* mPlotter;

	Canvas& mCanvas;
    OptimisationRun& mData;
};

#endif // TREEVIEW_H
