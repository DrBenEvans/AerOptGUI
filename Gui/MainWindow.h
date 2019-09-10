#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "Optimisation.h"
#include "ProfileModel.h"
#include "MeshModel.h"
#include "OptimisationModel.h"
#include "ProfileView.h"
#include "ParallelCoordinatesWindow.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void setOptimisationModel(OptimisationModel* optimisationModel);
    void setSplitterSizes();

public slots:
    void newOptimisation();
    void setCurrentOptimisationIndex(int index);
    void optimisationFitnessChanged(int index);
    void optimisationOutputChanged(int index);
    void onSelectedPointChange(int iGen, int agent);

private slots:
    void on_actionShowLog_triggered();
    void revealFilesCurrentOptimisation();
    void clearCurrentSelection();
    void setPointSelected(bool pointSelected);
    void setupGraphicsView();

    void on_loadOptimisationButton_clicked();
    void on_actionCurrent_Optimisation_Settings_triggered();

    void on_actionVisualSettings_triggered();

    void on_actionShow_triggered();

private:
    void setMeshViewSimulation(int iGen, int agent);
    void setLogText();

    /**
     * @brief openParallelCoordinatesWindow Opens a parallel coordinate in a new window
     */
    void openParallelCoordinatesWindow();


    Optimisation* currentOptimisation();
    Ui::MainWindow *ui;
    ProfileModel mProfileModel;
    OptimisationModel* mOptimisationModel;
    int mCurrentOptimisationIndex = -1;
    MeshModel* mCurrentMeshViewModel;
    ProfileView* mProfileView = nullptr;
    QString mClusterPassword="";

    ParallelCoordinatesWindow* pcWindow;
};

#endif // MAINWINDOW_H
